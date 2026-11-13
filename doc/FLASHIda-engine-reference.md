# FLASHIda real-time engine — reference

The full text behind `CLAUDE.md`'s one-line summaries, moved here on 2026-11-04 so it is read
when needed rather than loaded into every session. `CLAUDE.md` carries the rules; this file carries
the mechanisms, the measurements and the history behind them. ADR implementation notes for the
whole workspace are in the parent repo's `docs/adr/README.md`.

## Where the code lives

Headers under `src/openms/include/OpenMS/ANALYSIS/TOPDOWN/`, sources under
`src/openms/source/ANALYSIS/TOPDOWN/` — note the source tree has **no** `OpenMS/` path segment
(`source/OpenMS/…` does not exist). `sources.cmake` is **not** a reliable header inventory:
`CandidateAdmission.h`, `Ms2Params.h`, `NotchSelection.h` and `ProteoformTracker.h` are live but
unregistered (verified 2026-09-04 — the list is longer than it used to be). Glob the directory
instead. It costs nothing to register a new header there and it keeps the drift from growing;
`ScanRole.h`, `ScanCommandJoin.h`, `KineticModel.h`, `Refine.h`, `RefineSettings.h` and `SearchPool.h` are.
`FLASHIda/SearchPool.{h,cpp}` (ADR-0062) is the fixed pool of search workers with a bounded backlog that
matches MS2s beside the analysis thread; it knows nothing about spectra, and `FLASHIda` builds its jobs.
`FLASHIda/ProteomeSearch.{h,cpp}` (ADR-0061 d.3) is FLASHTnT's own search of one MS2: the Protein
reference loaded once (decoy-prefixed entries dropped and counted), `FLASHTnTAlgorithm::prepareDatabase`
built once and shared, and a `Worker` per search thread that owns its `FLASHTnTAlgorithm` and calls `run()`
on a one-spectrum input (the MS2 converted with `DeconvolvedSpectrum::toSpectrum`, its raw spectrum handed
to the annotation through a loader). With a target-only reference, every run-level step of `run()` (protein
FDR, cascade, q-values) has nothing to decide, so a spectrum's PrSMs are those a whole run of MS2s gives
it. `FLASHIda_ProteomeSearch_test` pins it on a private E. coli slice (`test-data-private/ecoli_parity/`,
parent repo only). Nothing in the engine calls it yet: push 4b wires it into the search pool.

### The untouchable boundary

The directory split is a near-perfect proxy for what you may modify.

| | Files | Status |
|---|---|---|
| **FLASHDeconv** (flat `TOPDOWN/`) | `DeconvolvedSpectrum`, `SpectralDeconvolution`, `FLASHDeconvAlgorithm`, `FLASHHelperClasses`, `MassFeatureTrace`, `PeakGroup`, `PeakGroupScoring`, `Qvalue`, `TopDownIsobaricQuantification` | **OFF-LIMITS** |
| **FLASHTnT** (flat `TOPDOWN/`) | `FLASHTnTAlgorithm`, `FLASHTaggerAlgorithm`, `FLASHGappedTaggerAlgorithm`, `FLASHExtenderAlgorithm` | **OFF-LIMITS** |
| **FLASHIda engine** | `FLASHIda.{h,cpp}`, `FLASHIdaBridgeFunctions.{h,cpp}` (flat) + everything in `TOPDOWN/FLASHIda/` | Fair game |

**"Fair game" means "fair game to propose".** The table says what is off-limits outright. Every
change under `src/` (engine, bridge, tests, `.cmake`, inside the boundary or outside it) needs the
owner's approval before it is written. An approved plan that names the files and describes the
change **is** that approval; do not re-ask per file. Anything the plan did not describe (a new file,
a changed mechanism, an extra test) goes back to the owner first. Ambiguous cases count as source
changes. See `../CLAUDE.md`, Working rules.

Two deliberate exceptions to the clean split, both worth knowing before you conclude the boundary
is airtight:

- **`OptimizationMetadata` is IDA code living in FLASHDeconv territory.** Its header sits flat
  under `TOPDOWN/` and it is stored as `std::optional<OptimizationMetadata> opt_metadata_` on
  `DeconvolvedSpectrum`, populated by `Exploration.cpp` and serialized through `toSpectrum()`.
  It is a sanctioned IDA hook with its own CI test (`DeconvolvedSpectrum_OptimizationMetadata_test`).
- **Only two IDA files reach into FLASHTnT**: `FragmentAnalysis.cpp` (constructs
  `FLASHTaggerAlgorithm` :430, `runMatching` :556, `FLASHExtenderAlgorithm` :566) and
  `PrecursorSelection.cpp` (constructs `FLASHTaggerAlgorithm` :934, calls the static
  `fillMatchedPositionsAndFlankingMassDiffs` :972). The IDA path drives Tagger + Extender
  **directly** and never uses `FLASHTnTAlgorithm` — that orchestrator belongs to the TOPP tool.
  Tuning knobs arrive via the `identification.flashtnt` config block
  (`config_.identification().flashtnt`, a `FlashTnTSettings`; ADR-0061 moved it out of the top level
  and out of `TargetingConfig`), not through `FLASHTnTAlgorithm`'s own defaults.
  **Every Param handed to either algorithm is built in one of two places** —
  `FragmentAnalysis::buildTaggerParam` / `buildExtenderParam` (:197 / :213, which construct a
  throwaway instance purely for `getDefaults()`). All three call sites route through them, so a new
  knob is added once, not per site. `fixed_mod` is set **unconditionally** there: an empty config
  list is passed verbatim as an empty Param list, *not* left at the algorithms' declared `{""}`
  placeholder — pinned by `FragmentAnalysis_test`'s `build*Param_sets_fixed_mod_when_empty`. Since ADR-0067 both algorithms
  **read** it (`FixedModifications::fromNames`), and since its push 2 `Config` accepts an authored list and builds it
  into `IdentificationConfig::fixed_mods`: the one table every FLASHIda ladder (`computePTMAdjustedFragmentMasses`,
  `runTagBasedFragmentMatching_`'s flanking masses and `vec_pro`, `matchFromPrSM`), MS3 match
  (`MS3FragmentMatcher::computeTheoreticalMasses`/`computeProteinPrefixMasses`/`calibrateAndScore`), Refine ladder and
  logged ProForma reads. A fixed modification is never a PrSM modification, so a ladder that skipped the table would be
  off by its mass on every fragment with the residue. The static readers default the table to `none()` for the tests;
  every production caller passes `config_.identification().fixed_mods`. The logged `proteoform` columns
  (`identification.tsv`, `pooled_identification.tsv`, `scan_commands.tsv`'s `ms3_proteoform`) are prefixed with its
  ProForma global modifications (empty: none). Under `quantification.scan` ms3/both a fixed C is refused: the ADR-0064
  label already fixes every C.

> **Name trap.** `TopDownIsobaricQuantification` (flat, off-limits, used only by
> `FLASHDeconvAlgorithm`) and `FLASHIda/Quantification` (IDA, fair game) are unrelated — but the
> IDA one now **borrows its vocabulary**: `quantification.labelling` takes the same seven strings
> that tool's `setValidStrings("type", …)` pins, minus `"none"`, because here
> `quantification.enabled` is the switch. All seven schemes are live (ADR-0038); the
> `TMTSixPlexQuantitationMethod` hardcode and its `// TODO: Variable channel extractors` are gone.
>
> ⚠️ Correction is applied now (`IsobaricQuantifier`, stock matrix, `normalization` left at its
> `false` default), and **an all-zero `correction_matrix` is handled by FLASHIda rather than passed
> through**: `IsobaricIsotopeCorrector` *throws* on an identity matrix, so delegating would abort
> every scan instead of skipping the correction.
>
> ⚠️ **`channel_normalization` (ADR-0057) is applied in `Quantification::evaluate()` and nowhere
> else** — subtracted in log2 *after* the empty check, so a factor never decides whether a channel
> was measured, and never written into `Result::channels`, which is what `scan_results.tsv` logs.
> `decide()` needs no code of its own: it already hands its aggregate to `evaluate()`. Moving the
> division to the top of `evaluate()` (or into `measure()`) would revalue the logged channels and
> make the offline pipeline normalise twice; `Quantification_test`'s
> `a_factor_shifts_the_ratio_and_leaves_the_reported_channels_raw` and
> `emptiness_is_judged_on_the_raw_value` pin both.

⚠️ **Five owner-approved edits sit inside the boundary.** The fifth (ADR-0074 d.10, 2026-11-11) is one appended field,
`FLASHAnnotatorAlgorithm::Result::explained_intensity`: the intensity of the peaks `annotateOne_` already flags in
`terminal_peaks` (every peak of an accepted terminal envelope), summed once, so envelopes that share peaks, as z, z+1 and
z+2 do, are not counted twice. No existing value changes. The first is `PeakGroup::updateSNR_` (`fix(TOPDOWN)`
`9721595ca4`, 2026-10-26, on top of the FDmerge merge `d0954cf37b`). FDmerge's `updateChargeRange_` passes
over up to two weak charges, so a charge with **no peak** can now sit inside a PeakGroup's range. With
no signal and no noise, the two `1e-6` regularisers cancel and that charge scores SNR **exactly 1.0**,
above the real charges. It then became the representative charge, and `buildMS2` wrote an isolation
centre of −5.5. The fix gives a charge with no signal SNR 0. It is a separate commit so Kyowon can take
it into FDmerge or reject it; nothing was pushed to `FDmerge` itself.

The third (ADR-0061 d.3 as amended for push 4c, 2026-11-07) is one line, `FLASHTnTAlgorithm.cpp:1131`
`setLogType(CMD);` → `if (! prepared_) setLogType(CMD);`: a prepared `run()` keeps its caller's log type, so the
real-time workers (`ProgressLogger::NONE`) print nothing per MS2; an offline run never prepares and is unchanged.

The second is `FLASHTnTAlgorithm`'s **prepared database** (ADR-0061 d.3, 2026-11-06), so that the real-time
search builds `run()`'s database side once instead of per MS2 (the gapped-tag index alone is seconds for
E. coli). `prepareDatabase()` holds that side moved as it was: the entries I→L, the sequences with and
without X, the decoy factor, the full sequences and fixed masses of a PEFF database, the blind-mod map, the
protein mass sets, the residue sums and their indexes, the gapped-tag index and the mass ladders.
`setPreparedDatabase()` makes `run()` reuse it; without one, `run()` builds its own, so offline FLASHTnT runs
the same code. Three things changed for offline runs, none of them a PrSM: the database work and its log
lines come before the spectra are read; the precursor-ladder index is built whether or not a spectrum has a
precursor mass (the searches still gate on it); and a prepared `run()` binds the ladders `const` and never
batches precursor windows (`preparePrecursorWindows`/`clearPrecursorWindows`), because the batch writes the
shared index. An unlisted mass is swept, which gives the same windows, and one spectrum never reaches the
256 masses the batch needs. `FLASHIda_ProteomeSearch_test` pins `run()` after the edit to a capture taken
before it (`test-data-private/ecoli_parity/tnt_reference.tsv`, run 42290157689), per scan as a multiset, because
`run()` is not bit-reproducible: modification masses differ in their last digits between two runs of one binary, a
scan's PrSM order differs between builds, and RawScore's chance rate is pooled per spectrum over its PrSMs' decoy
fragments (`FLASHAnnotatorAlgorithm::annotate`), so one decoy match at a tolerance edge moves every RawScore of that
spectrum (scan 1005, 0.4–1.0 %, in one of four runs). PrSM identity is compared exactly, RawScore to rel 2e-2, and the
score only by 8 × its RawScore change.

The fourth is ADR-0067 (2026-11-09): **fixed modifications are one residue-mass table**, `FixedModifications`
(`TOPDOWN/FixedModifications.h`, outside the boundary). Every residue-mass read in the Tagger, GappedTagger,
Extender, `FLASHTnTAlgorithm`, `FLASHAnnotatorAlgorithm`, `MassLadderAlgorithm`, `ProteoformDatabase` and
`FLASHTnTFile` goes through `mass(c)`: a residue that is not fixed is read live from `ResidueDB`, as before, so an
empty table changes no double; a fixed one is a `Residue` copy with `setMonoWeight(full + diff)`, the mass the
FLASHTnT tool's former process-wide overwrite gave it. `fixed_mod` is a parameter of `FLASHTnTAlgorithm` (handed to
the tagger, gapped tagger, extender and annotator, like `ion_type`), of the Annotator and GappedTagger, and the
Tagger's and Extender's old declared-but-unread `fixed_mod` is now read. Static readers take the table as a trailing
argument defaulting to `FixedModifications::none()`. A `PreparedDatabase` records the table it was built with, and
`run()` throws on another. The FLASHTnT tool passes `-fixed_mod` through and no longer touches `ResidueDB`; the
FLASHAnnotator tool reads a FLASHTnT file's table from its ProForma prefix (rows must agree) and has no option for it.
`FLASHIda_ProteomeSearch_test`'s last section pins the parameter route to the old overwrite on real data. Every other
FLASHDeconv/FLASHTnT file is exactly as FDmerge has it.

### Ownership tree — declaration order is load-bearing

`FLASHIda` is the single orchestrator and owns all ten subsystems **by value**. C++ initializes
members in declaration order and several bind references to earlier siblings, so **reordering the
member declarations is undefined behavior**, not a style choice:

```
config_ (Config) → logger_ (IdaLogger) → queue_ (ScanCommandQueue) → deconv_ (Deconvolution)
  → fragments_ (FragmentAnalysis) → selection_ (PrecursorSelection) → quant_ (Quantification)
  → [analysis_mutex_, exploration_active_, current_faims_cv_] → faims_ (FAIMS)
  → exploration_ (Exploration) → tracker_ (ProteoformTracker)
```

Reference-holding: `PrecursorSelection{const Config&, Deconvolution&}`,
`Exploration{const Config&, FragmentAnalysis&}`, `ProteoformTracker{const Config&, IdaLogger&}`;
config-reference-only: `IdaLogger`, `ScanCommandQueue`, `FragmentAnalysis`, `Quantification`.
`FAIMS` and `Deconvolution` take `const Config&` but **copy** what they need. The ctor init list
mirrors the declaration order exactly (`FLASHIda.cpp:54-64`).

The search pool is not one of the ten: `search_pool_` is a `std::unique_ptr<SearchPool>`, declared last,
created in the ctor body only when an `identification` section is present (ADR-0062), and **stopped and
joined first in `~FLASHIda`'s body**, before any member a running apply touches (tracker, queue, logger) is
destroyed. Waiting jobs are discarded; a running one finishes, its apply included.

### Enqueue discipline: return, don't push

`FLASHIda.cpp` is the **sole enqueuer**, and since ADR-0060 its sites come in three kinds, pinned as a
table by `FLASHIda_Replan_test::push_site_inventory` (a new or moved site fails it until the table says
what the site is): a survey's MS2s enter through **one** `queue_.replacePlan(` in `runSurveyMS1_`; with
`characterization.per_proteoform` on, a governed MS3 plan enters through **one** `queue_.replaceMs3Plan(` in
`pushGovernedMs3_` (ADR-0066, below; the `stampAndPush_` MS3 sites are then its off branch); five
**permanent** sites go through `stampAndPush_` (`FLASHIda.h`, the one `queue_.push` that stamps a
precursor id first: refine re-dispatch, exploration MS2 return, the bought `R`, MS3 targeting,
exploration MS3 return; the tag follow-up `C` went with ADR-0061; since ADR-0062 "MS3 targeting" sits in
`finishRegularMS2_` and "exploration MS2 return" in `finishVariantMS2_`, which a search worker's apply
reaches as well as `processScan`); three MS1-family `queue_.push(` sites (CV transition, cycle time,
monitor) plus the idle survey's `pushAndDequeue`; and one `registerAdjustable` (the quant verdict). No
`.push(` on a `ScanCommandQueue` exists anywhere else under `TOPDOWN/`. The `ScanCommandQueue&` that
`Exploration::initiate/feedResult/applyVariant/initiateNextLevel/applyNextLevel` receive is a **command
builder, not an enqueue target**: `buildMS2`/`buildMS3`/`buildFollowUp` allocate a tracking id and return a
`ScanCommand` by value. Any new command-producing component must follow the same contract.

They also do **not** register in the pending-scan map — that happens only in `dequeue()` and in
explicit `registerPending()` calls, which `FLASHIda.cpp` makes itself on the drain path.
⚠️ The doc comment at `ScanCommandQueue.h:59-60` claiming builders register is **stale**.

### A survey's MS2s are a plan until the drain commits them (ADR-0060)

```
runSurveyMS1_ (analysis thread, analysis_mutex_)          getNextScanCommand (event thread, no analysis lock)
  applyStamps(takeCommitted())  fold what the drain committed   dequeue_ -- commit the unit (queue_mutex_)
  filterAndRank -> journal(); rollbackAcquisitionMemory()       its stamps --> inbox_ --> takeCommitted()
  units per nominal mass (+ remainder re-plan, bought-R adjust)
  replacePlan(cv, units, adjust, replace_remainder_of)  -- one queue_mutex_ hold
  then, OUTSIDE the lock: closeGroup / dropGroup / forgetTracking_ for dropped AND reconciled units
```

- **`memory_` holds committed state only.** `filterAndRank` still writes acquisition memory as it admits
  (the admission logic was not touched); `rollbackAcquisitionMemory()` then restores the snapshot taken on
  the MS1 path, and the writes survive only as the journal's **stamps** (`PlanUnit.h`). The rollback is
  guarded by `acq_entry_valid_`, so it is a no-op after a non-MS1 call and can never restore an older
  survey's state. Survey-level stamps (refused candidates' `AllMassRt`) are held to the next survey's fold,
  which is equivalent because nothing reads them in between.
- **The commit rule.** A species commits when any of its commands is dequeued, and a committed species is
  never dropped. `dequeue_` commits every charge the command isolates (`isolatedCharges`: anchor + stage-0
  notches), so a multiplexed authored set cannot leak re-selectable charges. Only a non-quant `separate`
  remainder (re-derived by `rechooseCharges`, sent charges excluded) and the bought `R` stay adjustable.
- **`replacePlan` reconciles at the swap.** A new unit whose species committed after the survey read memory
  is refused by `wouldBar`, which must mirror `admitCandidate`'s exclusion gate and qscore bar — a
  bidirectional drift-guard test pins the two together. The last `MassQscore` is taken **by vector
  position**, so `committed_since_fold_` stays in commit order. A remainder replacement inherits the old
  unit's `sent_charges` *at swap time* and drops any charge sent in between; if the old unit emptied, it is
  refused. Reconciled units are reported in full (`reconciled_scan_ids` / `_quant_groups` /
  `_sweep_groups`) because their `addMember` / `initiate` side effects already happened and must be undone
  like a dropped unit's.
- **The bought `R`.** Registered with `ScanCommandQueue::isolatedCharges` of the group's measured `Q`s
  (`GroupDecision::q_covered_charges`) and the species' nominal mass **at registration** — the adjust keys on
  that, never on the `R`'s current `mono_mass`, which `reisolate` rewrites. An empty covered set registers
  nothing. `reisolate` shares `buildMS2`'s geometry byte for byte (`applyIsolation_`) and rebuilds the
  description's `@charge`, because `replacePlan`'s adjust does not touch the description.
- **Registry hygiene.** `cleanupExpired` and `cancelByScanIds` both call `forgetQueued_`, so an expired or
  cancelled member leaves its unit and an emptied staged unit is discarded with its stamps (never acquired,
  never excluded). An unregistered push whose parent sits in a still-staged unit prints `[REPLAN-INVARIANT]`.
- **`decision_info_` is consumed at the drain** (`takeDecisionInfo` erases), so a test must read
  `decision_survey_id` from `scan_commands.tsv`, not from `FLASHIdaTestAccess::decisionInfo` after a drain.
- ⚠️ **`Exploration::initiate` can return empty before it mints a group id**, so `lastGroupId()` then names
  the PREVIOUS group; `runSurveyMS1_` records a sweep group only after a non-empty return.
- ⚠️ **Depth-1 drives never re-plan.** `runInterleaved` and the C# harness return every MS2 before the next
  survey, and no committed config sets `cycle_time`, so outside FAIMS cycling the goldens cannot see this.
  `FLASHIda_Replan_test` drives stale surveys directly and forces the race through the private
  `beforeReplan_` seam (`FLASHIdaTestAccess::onBeforeReplan`), which commits in the gap every round.

### A proteoform's MS3s are one live plan (ADR-0066)

With `identification.characterization.per_proteoform` on (default off), each proteoform has one live MS3 plan: a
newer identification replaces its unsent MS3s, a sent MS3 commits its fragment, a sent fragment is acquired again
only at 2× its sent intensity, and a proteoform gets at most `max_targets` distinct cleavages per run.

```
dispatchFromModel (analysis thread, analysis_mutex_)            getNextScanCommand (event thread, no analysis lock)
  identified model + key on: key = accession, region, ms_ctx->mono_mass   dequeue_ -- commitMs3_ (queue_mutex_)
  snap = queue_.ms3Snapshot(key)                                            sent[cleavage][z] = max(intensity)
  planNextScans(pid, &gate) -> NextLevelResult::governed (items)       cleanupExpired -- expired uncommitted item
pushGovernedMs3_                                                            --> ms3_expired_inbox_
  setPrecursorForTracking_ every new id; beforeMs3Replan_ (test seam)
  replaceMs3Plan(key, items, gen, max_targets)  -- one queue_mutex_ hold
  then, OUTSIDE the lock: forgetTracking_, ms2_context_cache_, undoMs3Items_, [MS3-REPLAN], operator line
```

- **Where it lives.** `FLASHIda/Ms3Registry.h` is header-only and pure, like `PlanUnit.h`: `Ms3ProteoformKey`
  (accession, region with −1 normalised to 0/P, commanded mass), `sameProteoform` (same accession and region;
  masses within the MS1 tolerance or `sameSpeciesMass`), `Ms3PlanItem` (one fixed-CE command, one MS3 exploration
  group, or one (fragment, charge) under `separate`), `ms3Admits` (the 2× rule), `ms3NewCleavageRoom`,
  `Ms3ReplanResult`, `Ms3Expired`. The registry is in `ScanCommandQueue` (`ms3_entries_`, `ms3_item_of_`,
  `ms3_expired_inbox_`), the gate in `ProteoformTracker` (`Ms3PlanGate`, `gateAdmitsCharge`, `CleavageRoom`,
  `withdrawMs3`), the governed branch in `Exploration::dispatchFromModel` (`GovernedMs3`), the swap and the undo in
  `FLASHIda` (`pushGovernedMs3_`, `undoMs3Items_`, `undoExpiredMs3_`).
- **The lock is `queue_mutex_`.** `ms3Snapshot`, `replaceMs3Plan`, `takeMs3Expired`, the commit and the expiry all
  run under it; the undo runs on the analysis thread after the swap, outside it. The lock order stays analysis →
  queue, and the drain takes no new lock (ADR-0025).
- **The gate, at planning.** Per charge: m/z > 0 and |z| ≥ the charge floor (the two skips `dispatchFromModel`
  applied after planning, moved in so they no longer burn a slot), a window SNR > 1 over the dispatching MS2 in
  `buildMS3`'s window (`max(iso_width, 2.0)` Th, so it equals the later `window_snr`), and `ms3Admits`. A fragment
  with no admitted charge never enters a pool; the anchor is the strongest admitted charge and the notches come
  from the admitted charges. The cap: a cleavage never sent takes one of `max_targets − |sent|` slots, and once they
  are used a new cleavage is **skipped** (`add_target` returns `true`, so callers do not `break`); a re-acquisition
  takes none. `planExhaustive_` and `planDiscriminating_` apply both before any stamp (ADR-0023 D-d). A null gate is
  byte-identical to the ungated plan.
- **The commit at dequeue.** `dequeue_` calls `commitMs3_` after the plan-unit block. An item's first dequeued
  command commits it: `committed`, `++generation`, `sent[cleavage][z] = max` for every charge it acquires (anchor +
  stage-1 notches), a commit-log entry and its (owner, nominal). A committed item is never withdrawn; an MS3
  exploration group commits at its first variant and runs to completion.
- **The swap.** `replaceMs3Plan` matches the entry (closest anchor mass, then lowest id; created if none), withdraws
  every uncommitted item (queued commands only), then inserts the new items in plan order, each re-checked with
  `ms3Admits` and the cap against `sent` **now**; a refused item is reconciled and its `ms3_cmd_proteoform_` stash
  erased. It reports withdrawn, reconciled and inserted ids, `committed_in_gap` (commit-log entries after the
  snapshot's generation) and the entry's committed (owner, nominal) pairs. An empty governed plan still swaps, so it
  withdraws. A slot refused at the swap is not refilled until the next identification.
- **The undo.** `undoMs3Items_` groups items by owner and calls `tracker_.withdrawMs3`: the nominal leaves
  `dispatched_nominal_masses`, and an undecided discriminating round's record (nearest within the level-2
  tolerance) is un-dispatched, not marked returned-empty. An item keeps both when a sibling of the same owner and
  nominal is committed or was inserted by the same swap. `RoundComplete` (records still dispatched, all returned)
  makes the caller run `decideRefine`: decided, not re-dispatched. Every withdrawn sweep group is dropped
  (`exploration_.dropGroup`; a failure is an `[MS3-REPLAN-INVARIANT]` Warn event), then `publishExplorationState_`.
- **The expiry inbox.** `cleanupExpired` moves an expired registered, uncommitted item (its other variants taken off
  the queue) into `ms3_expired_inbox_`; a committed item's expired variant only leaves the registry.
  `undoExpiredMs3_` runs at the top of `processScan`, `applyRegularSearch_` and `applyVariantSearch_`, before
  planning, and undoes the items with their entries' committed nominals as of the take (`Ms3Expired`).
- **Inert with the key off.** No governed branch runs, the registry stays empty, `commitMs3_` returns on the empty
  map, and `undoExpiredMs3_` returns at once. The `cleanupExpired` leak of an expired MS3's MS2 context (and a
  stranded round or group) is unchanged with the key off.
- **Events, not columns.** `[MS3-REPLAN] key= entry= owner= withdrawn= inserted= reconciled= committed_in_gap=` is an
  Info event (ADR-0065); `[MS3-DISPATCH] commands=` counts what was planned, before the swap. The governed batch's
  operator MS3 line is written at the swap from the inserted ids (ions, the first command's activation and CE, the
  scan count); a batch withdrawn or refused whole has none. `commands_pushed`/`child_ids` count inserted ids only.
- ⚠️ **Ungoverned with the key on**: an unidentified exhaustive model (no key, ADR-0055) and the production MS3 of
  a committed MS3 exploration group, both through `stampAndPush_`.
- ⚠️ **Not covered**: `cancelByScanIds` (the stop path) and `closeOut` leave the MS3 registry alone; entries are
  never evicted (the cap is per run); a refine flip that changes the region makes a new key, and the old entry's
  queued MS3s are not withdrawn.
- ⚠️ **Depth-1 drives drain an MS3 (priority 1) before the next MS2**, so the goldens see the gate, the 2× rule and
  the cap. `FLASHIda_Replan_test` forces the gap commit through `beforeMs3Replan_`
  (`FLASHIdaTestAccess::onBeforeMs3Replan`).

### An MS3 is admitted by its predicted bond yield (ADR-0069)

With `identification.characterization.yield_gate` present (and, for the quantification `Q` MS3,
`quantification.ms3_yield_gate`), a target is admitted only at a charge whose predicted bond yield
E = p × (L − 1) reaches `threshold` (default 1). Absent, every plan is exactly as before.

- **Where it lives.** `FLASHIda/YieldGateSettings.h` holds the config (threshold + six coefficients, core-5 defaults),
  `FLASHIda/Ms3YieldModel.h` the pure model: `predictBondYield` (clip, logistic, × (L − 1); L ≤ 0 untested),
  `yieldAdmits` (≥), `yieldInputsOf(ms_ctx, scores, ion)` for the planner and `yieldInputsOf(ms3)` for the row,
  `ms3IonOf` (the one ion parser), `logYieldVerdict` (`[MS3-YIELD]`).
- **The planner.** `Ms3PlanGate::yield` / `ms_ctx`; `yieldWalk` runs after `chargesOf` / `gatedRawCharges` and before
  `CleavageRoom::admit`, in `add_target` and in the exhaustive and discriminating emit loops, so a refused target
  takes no place and no slot. `dispatchFromModel` builds a snapshot-less gate for an ungoverned plan when the yield
  gate is on (it then gains the isolatable and charge-floor tests in the walk). `dispatchQuant` walks
  `labelledFragmentCharges` (SNR > 1, then yield).
- **The ranking** (amended 2026-11-11). `YieldGateConfig::rank_by` is `Yield` by default. Mapped plans sort through
  `rankCandidates` (intensity, then every verdict made, then a stable sort by the anchor's E, refused last), used
  for each candidate list and for the final dispatch order; raw-mass plans evaluate every pool entry and call
  `rankPoolByYield` (tested by E, then untested, then refused, stable); `dispatchQuant` collects the buyable
  fragments and takes the first n by E. `Intensity` runs the code path the gate had. `yieldWalk` returns the
  anchor's E and whether the test ran through two optional out-parameters.
- **`characterization.order`.** `Ms3SendOrder::Random` permutes the vector `planNextScans` returned, in
  `dispatchFromModel`, with `sweepPermutation(n, kSweepSeed ^ ms_ctx->scan_id)`. The planner never sees it.
  `Ms3SendOrder::Interleaved` (ADR-0073) is made in the planner, because the groups need the proteoform model:
  after `rankCandidates(targets)`, `planNextScans` calls `interleaveTargets_`, which builds each target's groups
  and hands them to the pure, static `ProteoformTracker::interleavedOrder(groups_of, served, counted)`. Each turn
  goes to the target whose least-served group has been served least; a tie goes to the higher-ranked target; a
  target in no group follows every grouped one. Groups: under `coverage` and `fragments` the two termini (first
  letter of `ion_type`: a/b/c or x/y/z); under `ambiguity` every modification with `candidate_start <
  candidate_end`, and a target belongs to each one it contains (`fragmentContains`, the planner's own test), so a
  container of two serves both. Under a governed plan (`per_proteoform`) `served` starts from the keys of
  `snap->sent`, read from the key and not from `m.fragments` (a sent cleavage need not be a fragment of the current
  identification; a raw-mass key `u…` belongs to no group), and a target already among them is `counted`: it takes
  its turn and serves nobody a second time. Ungoverned, every plan starts at zero. The emit loop is unchanged, so
  under `fragment_charges: separate` a fragment's charges stay together and count once. `Config` refuses the value
  under `exhaustive` and `discriminating`, the latter for the whole run, including the ambiguity plans made after
  the round is decided.
- **The row.** `writeScanCommandRow` recomputes `predicted_bond_yield` from the command under the gate of its role
  (`Q` → quantification, any other MS3 → characterization); −1 otherwise. No map carries it: the planner's inputs are
  what `buildMS3` copies into the command (`FLASHIda_Ms3YieldModel_test::row_inputs_equal_plan_inputs`).
- ⚠️ **A discriminating round whose every record is yield-refused dispatches nothing and stays undecided**, as one the
  cleavage cap fills: a refused record stays in the pool (ranked last under `rank_by: yield`) and is skipped in the
  emit loop, so `dead_end` does not fire.

## `processScan` — five admission gates

`FLASHIda.cpp:79-169`. Before any MS-level branching, a scan must clear all of these; each returns
`0`. This is the always-on **engine-id-echo** contract: a spectrum is analysed only if the engine
itself minted and dispatched its id.

| # | Gate | Trace |
|---|---|---|
| 1 | `scan_description` (or null) shorter than 3 chars | silent |
| 2 | `!roleObserves(echo_role)` — AGC calibration scan; resolves pending, returns. Reads the ECHO, the one gate entitled to (ADR-0042); `Agc` is the only `observes == false` row | silent |
| 3 | `peekPending(decode(desc[0..2]))` empty — id never emitted | `[TRACK-RESOLVE] status=not_found` |
| 4 | `resolvePending` empty after a successful peek (upstream race) | `status=context_lost_race` |
| 5 | context support: `required_stages` = 0/1/2 for ms_level 1/2/3; rejects if ms_level ∉ {1,2,3}, or `parent_ctx.msn_level != ms_level`, or `parent_ctx.num_stages < required_stages` (a **less-than**, not an equality) | `status=context_unsupported` |

A sixth early return follows the gates: MS1 with `config_.level(1).selection == None` returns 0
before selection. None of these paths reaches the TSV writers, so an un-commanded scan produces no
log row — but gates 3–5 do leave an `events.log` trace.

**Gate 3 also refuses a result that arrived too late** (FDmerge `96f8ab3732`). Only an answer
(`resolvePending`) or an exploration abort (`cancelByScanIds`) empties `pending_scan_map_`, so a scan the
host or instrument dropped, or one an MS1-only replay never feeds, used to stay for the life of the engine.
Every registration is now numbered (`registerPending_`, called from `dequeue_` ahead of ADR-0060's commit
block). Every `IN_FLIGHT_SWEEP_EVERY` (1024) registrations, entries more than `MAX_IN_FLIGHT` (10 000)
registrations old are dropped with a `[TRACK-LOST] id=…` event. A result arriving after that is
`not_found`. The age is counted in commands **sent**, never in time, so an instrument pause cannot age an
entry. It is deliberately not tied to `scheduling.scan_timeout`, which only drops **queued** commands
(`cleanupExpired`).

**`processScan` takes retention time in SECONDS** (ADR-0051). The C export converts the host's
minutes, and it is the only production caller; the C++ tests call `processScan` directly with
seconds. Everything behind it decides in seconds — `rt_window`, the expiry of acquisition memory,
inclusion rows, target logs. The two log streams that record it (`scan_results` `rt`,
`ida.log` `RT`) write `rt / 60`, because their readers read minutes; the event markers
(`[TRACK-RESOLVE] rt=`, `[CHARGE-SET] rt=`) print seconds.

### A matched MS2 is searched beside the analysis thread (ADR-0062)

With an `identification` section the engine owns a `SearchPool` of `identification.search_threads` workers
(unset = the core count minus two, at least 1; each worker sets OpenMP's team size on start, because a new
thread does not inherit it). Two job kinds exist, exactly the matcher's two callers on an MS2:

| Scan | Under `analysis_mutex_` in `processScan` | On the worker, no lock | The apply, under `analysis_mutex_` |
|---|---|---|---|
| Regular MS2 (not `Q`) | deconvolve, then copy the stored MS2, `parent_ctx` and the half-filled row into a `RegularSearch` | `Exploration::searchNextLevel` | `applyRegularSearch_` → `finishRegularMS2_` (`applyNextLevel`: tracker seeding, MS3 planning and pushes), then the `scan_results` + `identification` rows |
| MS2 exploration variant the matcher reads (`variantNeedsSearch`) | `Exploration::prepareVariant` (route, deconvolve, share) into a `VariantSearch` | `Exploration::searchVariant` | `applyVariantSearch_` → `applyVariant` (score, pool, group completion) → `finishVariantMS2_`, then the rows |

`initiateNextLevel == applyNextLevel(searchNextLevel(...))` and `feedResult == applyVariant(prepareVariant(...))`,
split at their existing seams with no logic change; without a pool, or for a `Q` scan, an MS3 variant, a trap
or empty variant, the inline path runs both halves as before. The cascade re-match of a completed MS2 sweep's
winner (`initiateNextLevel` inside `applyVariant`) still runs under the lock.

- **`analysis_mutex_` has two acquisition sites**: `processScan` and the worker's apply. The drain takes
  neither (ADR-0025). A worker never holds the pool mutex while it runs a search or an apply, so the
  hierarchy `analysis_mutex_` → { `queue_mutex_`, `precursor_map_mutex_`, the `IdaLogger` stream mutexes,
  the `SearchPool` mutex } stays acyclic.
- **Engine output goes through `EventLog`** (ADR-0065, `FLASHIda/EventLog.h`): process-wide, one mutex, one
  whole line per call. `EventLine() << …` composes an event and emits it at scope exit; `EventLog::say` writes an
  operator line (console + file). `IdaLogger` opens `events.log` with the run folder and closes it in its destructor;
  lines before that (Config's warnings) are kept, at most 1,000, and written first. Hazard: two engines alive at once
  share the one file; the second `openFile` replaces the first.
- **A quantification MS3 is bought on the apply** (ADR-0064): under `quantification.scan` `ms3`/`both`, an MS2 whose
  identification passes `reachesThreshold_` (and is neither charge-named nor already closed out) gets
  `Exploration::dispatchQuant` before its characterization MS3s are pushed: `ProteoformTracker::labelledFragments`
  ('C', any ion type, not full length), window SNR > 1 on the MS2, the most intense `ms3_targets`, `buildMS3` with
  `ms_settings.ms3_quant` and marker `'Q'`. It returns on the regular MS3 branch, is read by
  `Quantification::measure(..., 2, "ms3_quant")` (level 2: the extractor walks a one-spectrum experiment), and is never cached for annotation. The search fixes iodoTMT6plex on
  every C as an entry of the fixed-modification table: `Config` appends `iodoTMT6plex (C)` to `flashtnt.fixed_mod` under
  `readsAtMs3()` (ADR-0067 d.7), so FLASHTnT, every ladder, the MS3 matching, Refine and the logged ProForma carry it.
- **The operator view** (ADR-0065 push 2): `PrecursorBoard board_` is fed at the drain's send, an MS2's return,
  the apply and the survey, under its own mutex (lock order `analysis_mutex_` -> board; the drain and the render
  thread hold no `analysis_mutex_`). `ConsoleView::start` returns null unless stdout is a VT-capable Windows
  console; its render thread snapshots the board (`boardSnapshot_`) and paints through `EventLog::paint`.
  `~FLASHIda` stops the view, then the pool, then prints the run summary through `EventLog::summary`.
  A row describes its **best identification** in every value it shows (`onResult`, highest
  `search_identification_score`, ties to the later): that MS2's commanded mass and anchor charge, the score, protein,
  ProForma, and `cov`, the Pooled coverage of that identification's precursor id. Coverage reaches the board through
  `ProteoformTracker::setPooledObserver`, called right after each `pooled_identification.tsv` row is written
  (`emitPooledIDRow`) with that row's `precursor_id`/`nominal_mass`/`trigger`/`trigger_scan_id`/`coverage_pct`;
  `PrecursorBoard::onPooled` keeps the latest per precursor id and returns an `[OP-MS3]` line when an MS3 fold (an
  ion trigger, not `MS2` or `MS3-DECIDE`) raised it. The status line's `searching w/c` is `SearchPool::waiting()`
  against `backlog()` (the capacity); `identified` counts the Identified state only, closed out apart. The summary's
  search line gives `search_ms` and `search_wait_ms` (nearest-rank median, p95).
- **A close-out runs on the apply** (ADR-0061 push 5, `closeOutIfIdentified_`): an MS2 whose
  `search_identification_score` reaches `identification.threshold` arms `CloseOutMassRt`/`CloseOutMzRt` (max-applied)
  through `ScanCommandQueue::closeOut`, which also withdraws the species' uncommitted units at every CV under
  `queue_mutex_` (the lock order holds: `analysis_mutex_` → `queue_mutex_`). `closed_out_` then skips another
  Precursor's later scan of the species inside `rt_window` (`closed_out_skip`, finished as a dropped search); the
  closing Precursor's own scans are still searched. A charge-named inclusion species is never closed out.
- **A rescue set is bought on the apply** (ADR-0068, ADR-0070, `finishRegularMS2_`, right after the close-out
  verdict): with `ms_settings.ms2_rescue` (`IdentificationConfig::rescue_scans`, `hasRescue()`), a searched, non-Q
  MS2 that is short of the close-out threshold (`!reachesThreshold_`: no proteoform, or a score below it) on a species
  not closed out (`speciesClosedOut_`, no closing-Precursor exemption: a sibling charge may have closed it while this
  scan was searched) buys `buildFollowUp(parent_ctx, entry, 'R')` for every entry, in order, all at priority 0, at
  most once per **precursor id**. Hazards:
  - **The cap is per precursor id** (ADR-0070 d.3; ADR-0068 capped per detection). Under `separate` each charge is its
    own precursor id and buys its own set. `rescuePrecursorFree_` claims `rescued_precursors_`; every rescue carries
    its trigger's precursor id, so its own return finds the id taken. Depth is one: no rescue buys anything.
  - **Register before the push**, every scan of the set. `queue_.registerRescue(id, nominal)` goes into `rescue_`,
    which `closeOut` withdraws while queued (`withdrawn_rescue_ids`, one `[RESCUE-WITHDRAWN]` per scan). `dequeue_`
    and `forgetQueued_` erase it, so a sent rescue is never withdrawn. `replacePlan` never reads it. Both calls are
    rows of `push_site_inventory`.
  - **Every entry's activation has its own prepared search**, the ion-set block in the constructor; an ETD rescue
    searched on b/y would find nothing.
  - **Config:** a bare object (ADR-0068's shape), `[]` and two identical raw entries are refused; more than one entry
    with MS3 on and `per_proteoform` off loads with a `[CONFIG-WARN]` (each rescue's identification plans its own MS3).
  - **Events:** one per set, `[RESCUE] precursor_id=… nominal=… score=… reason=unidentified|below_threshold n=… ids=…
    trigger=…` (`ids` space-separated: the id alphabet contains `;`).
- **A Window scan is bought on the apply** (ADR-0071, `buyWindowScans_`, its one push site): with
  `identification.window_scan`, a searched, non-Q MS2 whose search found a proteoform (the identification row's own
  test; `finishRegularMS2_`, after the rescue block) buys, for each stage-0 window (the anchor, then
  `notchesForStage(ctx, 0)`), an MS1 command whose `first_mass`/`last_mass` are that window's commanded bounds widened
  by `window_scan_margin`, at priority 0, at the MS2's CV, with the MS2 as `parent_scan_id`, role `'W'`, pushed with
  `stampAndPush_` under the MS2's precursor id. The command is `makeMS1(*window_scan_config)` when
  `ms_settings.window_scan` is authored, else `makeMS1()`. Its description is `{id}W{mass token}@{charge}`
  (`writeWindowDescription`): the MS2's `mono_mass` and the charge of the window it reads. Hazards:
  - **A completed MS2 sweep buys too, once** (`finishVariantMS2_`). `Exploration` sets `FeedResultInfo::window_ctx`
    on the row that completes an MS2 group: the command of the first pre-scan with `found_proteoform`, unless the
    group emitted a level-2 command (its production scan, an `'R'` that buys for itself). The parent is the completing
    pre-scan, and the bought ids join its `child_ids` and `commands_pushed`.
  - **It decides nothing because its role does not.** `'W'` has `decides == false`, so ADR-0042's MS1 fork sends its
    return to `observeMonitoringMS1_`, and the drain does not reset the survey clock for it. `monitor_deconv_` is built
    when `window_scan` is on, as well as for monitor scans.
  - **It is registered nowhere.** No survey replaces it and no close-out withdraws it; it is uncapped.
  - **Every reader that treats "an MS1 that is not `'M'`" as a survey must also skip `'W'`:** both test harnesses (they
    answer it with the last survey cut to the command's range) and the offline `FLASHIdaWrapper` MS1 path.
  - **Events:** `[WINDOW-SCAN] precursor_id=… trigger=… windows=N`.
- **A sweep completes when its searches return**, not its scans. **Rows of searched MS2s are written in
  completion order** in production.
- **A bounded backlog drops the oldest** (d.5): 4 waiting jobs per worker (`kSearchBacklogPerThread`,
  provisional until push 4 logs `search_ms`/`search_wait_ms`). The dropped job's `drop()` runs in the
  submitting `processScan`, under its lock: `[SEARCH-DROPPED] id=…`, the row written with the identification
  columns at the not-run sentinels (`tag_count` -1), nothing matched, planned or identified. A dropped sweep
  variant is received unscored under `fragment_count`, scored as usual under the other metrics, and never
  pooled, so its group still completes.
- **`waitForSearchIdle()`** (bridge export `WaitForSearchIdle`, the sixth) blocks until nothing waits or runs
  and returns the commands the applies pushed since the previous call; `0` at once with no pool. Harnesses
  only: production never calls it, and calling it with `analysis_mutex_` held deadlocks (an apply would wait
  for the caller).

### Scan-description wire format

`{3-char base-94 tracking id}{1-char type marker}{payload}`. The marker is decoded in exactly ONE
place — `FLASHIda/ScanRole.h`'s `constexpr` traits table (ADR-0042), which
`IdaLogger::scanTypeFromDescription_` now delegates to:

| Marker | `scan_type` | observes | decides | |
|---|---|---|---|---|
| `S` | `survey` | ✓ | ✓ | survey MS1 |
| `A` | `agc` | **✗** | ✗ | AGC calibration; peaks discarded before the spectrum is read |
| `M` | `monitor` | ✓ | **✗** | monitor scan — observed, never acted upon (ADR-0042) |
| `R` | `recording` | ✓ | ✓ | any data-acquiring MS2 *or* MS3, **and the identification scan a quantification verdict buys** |
| `Q` | `quantification` | ✓ | ✓ | the quantification scan (ADR-0038) |
| `E` | `exploration` | ✓ | ✓ | exploration variant |
| `W` | `window` | ✓ | **✗** | Window scan: a full MS1 whose range is one window of an identified MS2; observed, never acted upon (ADR-0071) |

The `scan_type` strings are a **golden column** and the table is now their only source — `"recording"`
is a poor name for an identification MS2 and must stay. Pinned by
`FLASHIda_exploration_test::roleName_matches_the_legacy_scan_type_strings`, which also asserts both
predicates as EXACT sets. Adding an enumerator without its table row **fails to compile**
(`static_assert` on totality), which is deliberate: `ENABLE_GCC_WERROR` is off and MSVC's C4062 is
off-by-default, so a `switch` without `default` would only warn.

⚠️ **An UNDECODABLE description — too short, or a marker no row claims — resolves as a SURVEY**:
`observes` and `decides` both true. That is what keeps the migrated read sites byte-identical with
the `char` comparisons they replaced. Do not "fail closed" here; it would silently stop processing
every scan the engine did not mint.

⚠️ **Which copy of the role is authoritative: the ENGINE's** (ADR-0042, settling a question ADR-0008
and ADR-0035 left open). Two copies exist and `processScan` names them `echo_role` and `cmd_role`.
The **observe** gate reads the instrument's echo, because it must run before `resolvePending` and
there is no engine record yet. Every **decision** gate reads `parent_ctx`, the engine's own queued
command, which cannot be truncated or lost by the trailer round trip.

⚠️ **`F` is retired** (ADR-0038). It marked the scan quantification *bought* — which the engine
never measured — while the scan it *did* measure was the base MS2, whose activation could not
release the reporter ion. `Q` marks the scan quantification **measures**: rostered once per
selected precursor, and the only scan measured. What it buys is an ordinary `R`, which is what
keeps `R` meaning "identification MS2" in every mode; only *when* it fires changes.

⚠️ **WHICH verdicts buy is authored, not hardcoded** (ADR-0039). `quantification.identify` is a cut
point on the verdict ladder — `differential` (default; ADR-0038's exact behaviour) | `quantified` |
`all` | `none` — and `quantification.enriched_in` names the **condition** a differential species
must be enriched in (or `"either"`). Never a direction: `fold_change = gmean(conditions[0]) /
gmean(conditions[1])` (geometric means, ADR-0056), so `"up"` would mean *enriched in whichever condition is listed first*.
The predicate is `quantBuysIdentification` in `FLASHIda.cpp`'s anonymous namespace, and it reads
**`Result::condition_means`, never `fold_change`** — a wholly-absent condition is `Differential`
with `fold_change == -1.0`, a sentinel, so a `fold_change`-based direction test drops exactly the
species the experiment exists to find. Two fixtures pin it, one per absence direction, because with
`conditions[0]` absent a bug spelled `fold_change < 1.0` coincidentally agrees with the right
answer. Since identification and MS3 are suppressed on a `Q` and ride the bought `R`, `identify` also
governs what gets characterized — with no MS3 code involved.

`processScan` branches on the marker twice: `desc[3]=='A'` (gate 2) and
`is_quant_scan = desc[3]=='Q'` (`is_follow_up_scan`, for the retired `C`, went with ADR-0061). The
screen gate keys on
`'Q'` rather than on "is not a follow-up", which is what makes it safe for the bought scan to be an
indistinguishable `'R'` — an `R` is never re-screened, so depth-one holds without the marker having
to separate a follow-up from a primary.

> The ABI field is `char scan_description[256]`, but **every writer caps at `snprintf(dst, 16, …)`
> — 15 chars + NUL.** That is why MS2/MS3 mass tokens go through
> `ScanCommandQueue::formatMassToken`, which computes a token budget and degrades decimal precision
> from 6 down to 0 so the trailing `@{charge}{ion_type}{index}` is never truncated.

## `getNextScanCommand` — never returns 0

Six steps: (1) scheduled AGC prescan → (2) cycle-time MS1 → (2b) monitor scan (ADR-0042) → (3) cleanup → (4) priority dequeue →
(5) idle fallback. **Every path returns 1.** When all four queues are empty, Step 5 mints a fresh
priority-3 survey MS1, pushes it, and re-enters `dequeue()` to return it — the instrument is never
observably starved. The only `return 0` in the whole file belongs to `processScan`.
`FLASHIda.h:99`'s "0 if error" comment is imprecise.

Callers must therefore terminate on `msn_level == 1 && priority == 3`, or bound their iterations.
See `../CLAUDE.md`.

⚠️ **Step 2b is NOT an `else` of Step 2**, although the two are behaviourally complementary (Step 2
fires only when nothing is sweeping, Step 2b only when something is). `queue_.takeMonitorDue()` must
be reached on EVERY drain, because its `sweeping == false` path is what RE-ARMS a level for its next
episode. Put it in an `else` and the arm survives the end of a sweep, so the next sweep silently
loses its anchor scan — visible only as a missing first data point. Both levels are always asked,
even once one has answered yes, so each consumes its own clock: it is one MS1 and it serves both.
`takeMonitorDue` is a CONSUME, not a peek — it tests and updates under one `queue_mutex_`
acquisition, so two concurrent drains cannot both be told the same monitor scan is due. Its cadence
is one anchor per episode, then one per `interval_ms`; the anchor is deliberately independent of the
interval, so a sweep shorter than one interval still yields the reading that says what the source
looked like while it ran.

⚠️ **The drain tail's `recordMS1Time()` carries a `roleDecides` conjunct** (ADR-0042). A monitor scan
satisfies `msn_level == 1 && is_agc == 0` exactly as a survey does, so without it the scan stands in
for one — which is character-for-character ADR-0031's bug, with `scheduling.cycle_time` as the
victim instead of `agc_interval_seconds`: monitoring a long sweep would starve the very survey the
cycle time exists to force.

⚠️ **Step 5 emits NO prescan** (ADR-0031). It used to fabricate one as filler *and* call
`recordAGCTime()`, which reset the timer Step 1 reads — so Step 1 could only fire in a run whose
queue stayed busy for a whole interval, and `agc_interval_seconds` never governed the real cadence.
Step 1 is now the only `makeAGC()` caller. Two details of Step 5's shape are load-bearing: it
re-enters `dequeue()` rather than returning its local copy, which is what makes the survey inherit
`recordMS1Time()`, the timestamps, pending-map registration and the rich log row from the Step 4
tail; and that re-entry is also the correct read under a concurrent `processScan`, since work pushed
between the two dequeues wins on priority.

Two subtleties:

- **Step 2 does not return.** The cycle-time MS1 is skipped while `exploration_active_` is true;
  when it fires it stamps priority 0, pushes, and deliberately falls through to Steps 3–4.
- **The five steps are not atomic** — no lock is held across the body (explicit comment at :619).
  Each `ScanCommandQueue` call takes `queue_mutex_` independently, so one drain is many short
  critical sections. The drain takes **no** `analysis_mutex_` anywhere, the `writeScanCommandRow` calls
  included (ADR-0025; `IdaLogger` serialises each stream itself), and ADR-0060's commit happens inside
  `dequeue_`'s own `queue_mutex_` hold. Lock order wherever both are taken: `analysis_mutex_` → `queue_mutex_`.
  `analysis_mutex_` itself has two acquisition sites, `processScan` and a search worker's apply (ADR-0062).
  The only lock-free `processScan` → `getNextScanCommand` channel is two atomics
  (`exploration_active_`, `current_faims_cv_`), written release / read acquire.

This matters because the two bridge calls genuinely run on **different threads**: `processScan` on
the C# TPL Dataflow `ActionBlock` thread (serialized with itself, `MaxDegreeOfParallelism = 1`),
`getNextScanCommand` on the instrument's event thread, and since ADR-0062 a matched MS2's apply on a
search-pool worker.

### Priority ladder (0 highest → 3 lowest)

| P | Commands |
|---|---|
| 0 | FAIMS CV-transition MS1; cycle-time MS1; the **bought scan** — the identification scan a quantification verdict buys (`R`, ADR-0038; the tagging conditional follow-up `C` was retired with ADR-0061); the **monitor scan** (`M`, ADR-0042) — 3 is the idle sentinel and 1/2 are the sweep's own lanes, so 0 is the only priority at which its authored interval can govern |
| 1 | MS3, including exploration MS3 variants |
| 2 | MS2 from MS1 selection; exploration MS2 variants |
| 3 | idle survey MS1 — **and nothing else** |

The **AGC prescan is not in this table**: `makeAGC()` sets `priority = 0`, but Step 1 returns the
command directly and never pushes it, so that field is never read. Its ordering comes from where the
check sits — ahead of the dequeue — not from its value.

Priority 3 is a contract, not a detail: three C# drain loops treat a priority-3 MS1 as "the queue is
drained" (ADR-0031). Emitting anything else at 3 truncates them silently. Pinned by
`FLASHIda_ProcessScan_test::only_the_idle_survey_is_emitted_at_priority_3` ∥ C# `IdleSurveySentinelTests`.

Assigned at build time, never re-ranked. Within a priority `dequeue()` is strict FIFO, which is
what makes `scan_commands.tsv` row order and `child_ids` order deterministic. `push()` clamps
out-of-range priorities into [0,3].

## Where deconvolution actually runs — 4 call sites, 3 engines

- **MS1** is *not* deconvolved from `processScan`: it happens inside
  `PrecursorSelection::filterAndRank` → `deconv_.deconvolveMS1(…)`, just before mass-exclusion.
- **Regular MS2/MS3** — `deconv_.deconvolveMSn(…)`, called directly from `processScan`.
- **Exploration variants use a different engine.** `Exploration` owns a *second* `Deconvolution`
  instance built with `config.explorationToleranceList()` rather than `config.toleranceList()`.
  So `deconv_.storedMS2()` is empty/stale during an exploration group, and exploration results can
  use a different ppm tolerance than production scans.
- **Monitor scans use a THIRD engine** (ADR-0042). `FLASHIda::monitor_deconv_` is a
  `unique_ptr<Deconvolution>`, null unless a level enables `monitor_ms1`, built with the SAME
  `config.toleranceList()` as `deconv_` so a monitor reading is numerically comparable with a
  survey's. It is private on purpose, not incidentally: `deconv_` holds ONE `SpectralDeconvolution`
  shared between MS1 and MSn, and `deconvolveMS1` restores the global charge window on it (flagged
  `// LOAD-BEARING` in `Deconvolution.cpp`), so routing a monitor scan through `deconv_` mid-sweep
  would reset that window under the running sweep. It self-corrects — `deconvolveMSn` re-narrows on
  every call — but that is an inference across two call sites, and a private instance makes it
  impossible instead. **Do not "simplify" it away.** Which of the three engines carries a fed rt is
  also the only spectrum-independent way to tell which MS1 arm ran, which is what
  `the_engine_record_decides_the_role_not_the_echo` asserts on.

`target_mode` semantics (and the fact that 2/3 are swapped in every doc comment) are documented in
`../CLAUDE.md` — mode **2 is deep/in-depth**, mode **3 is exclusion**.

## Characterization / MS3

**`ProteoformTracker::planNextScans` is the only producer of MS3 targets.** `initiateNextLevel`
builds commands only through `dispatchFromModel`, which needs a model that is identified **or**
under `exhaustive` (ADR-0055); every other path returns zero commands. No tracker, no finalized
model, or an unidentified precursor in any other mode ⇒ **zero MS3 by design**, with no intensity
or legacy fallback.

- **`identification` gates MS2 matching; its `characterization` block gates MS3** (ADR-0061).
  `Config::identification()` returns an `IdentificationConfig`: `sequence` (from
  `identification.reference.sequence`, refused empty), `method`, `refine`, `has_refine_block`,
  `flashtnt` (`FlashTnTSettings`: `min_length`, `max_length`, `allow_gap`, `max_aa_in_gap`,
  `fixed_mod`, `max_blind_mod_count`, `max_mod_mass`), `fixed_mods` (ADR-0067: the table of `flashtnt.fixed_mod`,
  built at load) and `enabled()` (= non-empty sequence). Read
  sites: `config_.identification().sequence` (was `characterization().protein_sequence`),
  `.method`/`.refine`, `.flashtnt.*` (was `targeting().min_tag_length` etc.); the three MS2-matching
  gates in `FLASHIda.cpp` read `identification().enabled()`. `CharacterizationConfig` no longer
  carries the sequence, method or refine block. Root `characterization`/`flashtnt` and
  `identification.characterization.protein_sequence`/`.method`/`.refine` throw with a migration
  message.
- **`identification.characterization.mode`** (`ambiguity | coverage | fragments | exhaustive |
  discriminating`) picks the targets: each value *is* an objective, one for one (ADR-0013, extended
  by ADR-0023, ADR-0049, ADR-0052). It is required when the block is present, and `"off"` throws
  ("delete the block"); `CharacterizationMode::Off` remains as the internal state for an absent
  block, which means no MS3 at all.
  **No objective is a set cover any more** (ADR-0052): `coverage` takes every fragment whose interior
  still holds an MS2-unwitnessed bond, `fragments` every mapped fragment with a bond inside it, both
  strongest `best_ms2` first until `max_targets`. The uncovered set is computed once and never
  shrinks as targets are picked — the deleted `uncovered.erase` loop assumed an MS3 witnesses every
  bond inside its parent and planned 2 of 67 targets on ubiquitin. One `eligible()` lambda in
  `planNextScans` drops a **full-length ion** (`ion_index >= L`: b_L/y_L, the precursor itself) from
  the ambiguity, coverage and fragments pools; it sits where the pools are BUILT, not in
  `add_target`, whose `false` means "budget full" to callers that `break` on it. `planExhaustive_`
  keeps full-length masses. The objective dispatch is now an explicit `if / else if / else` — the old
  bare `else` meant any new objective silently ran coverage.
  **Unknown values throw** — the old `objective` parse was `if (== "coverage") … else Ambiguity`, so
  `"Coverage"` or any typo silently meant ambiguity, and with `mode` now carrying the on/off bit a
  typo'd `"Off"` would have silently *enabled* MS3. The struct still exposes
  `characterization().objective`, derived from `mode`, so every existing read site is unchanged.
  ⚠ **`mode` itself has zero read sites outside `Config.cpp`** — it is parsed, projected and quoted
  in error text there and nowhere else; the engine branches on `.objective`
  (`ProteoformTracker::planNextScans`, plus one `== Discriminating` test in `FLASHIda.cpp`'s MS3
  return path). A new mode that does not also assign an objective therefore
  ships byte-identical to `ambiguity`: accepted, green, and inert (ADR-0023 D-a).
- **`exhaustive`** (ADR-0023) targets **every deconvolved mass of the winner MS2 scan**, not only the
  masses that matched the winning proteoform — in the reference cytC MS2 that is 117 masses rather
  than 44. A mass that maps keeps its real `ion_type`/`ion_index` and is matched as today; a mass
  that does not is labelled `'u'`/`0`, acquired identically, and **logged rather than matched**.
  `'u'` is an in-engine sentinel only: it never reaches the wire, because `ion_index == 0` sends it
  down `buildMS3`'s no-ion branch. `MS3FragmentMatcher::isKnownIonClass` is the positive guard that
  makes every projection site refuse an unknown class **on the class, never on the index** — the two
  fields travel independently, and the old `default:`/`else` fallback was the *suffix* branch, so an
  index-only guard would let a `'u'` with a plausible index cut a real suffix and match against it.
  Its pool floor is the new `identification.characterization.min_target_mass` (Da, `0` = off), and an unassigned
  target's CE sweep is forced onto `RemainingPrecursor` whatever the config asked for — gated
  `msn_level >= 3`, because `initiate`'s `ion_type` defaults to `'\0'` and an ungated test would
  drag every MS2 exploration group with it (ADR-0023 decisions 5 and 11, D-f, D-g).
  **An unidentified Precursor is characterized too** (ADR-0055, amending ADR-0023 decision 10). Its
  winner scan is the scan the dispatch comes from — the regular MS2, the sweep metric's winner
  (`variants[best_idx]`, passed by the cascade) or the production follow-up — captured by
  `initiateNextLevel`'s pooled-baseline block through `finalizeMS2(pid, dispatch_scan_id)`, a
  fallback read only when nothing identified. The model gains `winner_peaks`/`winner_params` and
  nothing else: no sequence, no fragments, no pooled row, and `foldMs3` refuses it. Every target is
  therefore `'u'`/`0`. `planNextScans`' `unidentified_precursor` guard exempts `Exhaustive` alone, and
  `[MS3-PLAN] … identified=0|1` says which case ran. Under `fragment_count` the sweep winner of an
  unidentified Precursor is `ce_min` (every variant scores 0), accepted.
- **`Config::applyCharacterizationMode_`** projects `mode` onto `levels_[2].selection` and
  `levels_[3].selection` after the whole document is parsed. MS3 needs BOTH non-`None`
  (`FLASHIda.cpp:366`, `Exploration.cpp:728` and `:730`), so driving them from one enum makes
  "MS3 on with MS2 off" unrepresentable. It also assigns **level 1** from
  `precursor_selection.rank_by` — do not remove that: `MSLevelConfig::selection` defaults to `None`,
  so an unassigned level 1 makes `FLASHIda.cpp:168` short-circuit every MS1 and the instrument
  acquires *nothing*, silently.
- **`selection_strategy.ms3.selection` is DELETED** (ADR-0014). It used to be the on/off gate *and*
  the MS2-matcher chooser (`intensity`/`qscore` → `getTopFragmentMatches`, `terminal_fragments` →
  `getTerminalFragmentIons`, `ambiguity_resolution` → `getAmbiguityEnclosingIons`). The gate is now
  `mode`; the matcher is hardcoded to `getTopFragmentMatches`, which is what every MS3-enabled
  config selected anyway, so the other two matchers survive in the engine but are unreachable from
  config.
- **Budget and charge floor come from level 2**, not 3: `config_.level(2).max_targets` (hardcoded
  inside `planNextScans`) and level 2's `min_charge`. Both are now AUTHORED as
  `identification.characterization.max_targets` / `.min_fragment_charge` and projected onto level 2
  by `Config::applyCharacterizationMode_`, so the read sites are unchanged while the keys finally sit
  where the feature does. Quirk: level 2's `max_targets` also caps the matched fragments behind the
  MS2 `tic_coverage` column, so under the old `mode: off` with a sequence an authored value still
  counted there; an `identification` with no `characterization` block uses the default 3. `selection_strategy.ms3.max_targets` / `.min_charge` were parsed and
  **never read** — four committed configs set 200 and ran 3 — and are deleted.

### `fragmentContains` vs `fragmentBrackets` (ADR-0005)

Two similar-looking predicates in `ProteoformTracker`'s anonymous namespace, used for opposite
purposes — do not swap them:

- **`fragmentContains(f, rs, re)`** = `cover_start <= rs && cover_end >= re`. **Targeting only.**
- **`fragmentBrackets(f, rs, re)`** = cleavage strictly inside the range. **Narrowing only**
  (`narrowModifications_`).

Targeting on `brackets` is self-defeating — a bracketing ion would already have localized the mod —
and dispatched zero MS3 in every mode.

Ambiguity target ordering: keep mods with `candidate_start < candidate_end`; sort **widest range
first**; per mod collect containing fragments with a `best_ms2`, sorted by best-MS2 intensity
descending; **round-robin by rank** (every mod gets its strongest container before any mod gets a
second), deduped by `(ion_type, ion_index)`, until the budget fills; then re-sort by intensity.

> `selectNextLevelTargets` fills `masses/charges/wstarts/wends/ion_types/frag_indices/frag_scores`
> and **none of them is read again**. Only `frag_result` and `qscores[i]` (summed into
> `tic_coverage`) survive.

A returning regular MS3 is scored against the **live tracker winner**
(`tracker_.buildWinnerProteoformContext(precursor_id)`), not its own MS2 context — the cached MS2
context supplies only `fragment_ion_type`/`fragment_ion_index`. No winner ⇒ empty context ⇒ the
MS3 matches nothing. The cache entry is erased after use, so a duplicate/late MS3 yields no row.

⚠️ **An MS3 that arrives with zero peaks clears the stored deconvolution** (OpenMS `3ea17e49cb`). The
regular MS3 branch calls `deconvolveMSn` only when `mzs`/`ints`/`length > 0`, and `deconvolveMSn` is
what clears `storedMS2()` on entry — so an empty MS3 used to read back the **previous** MSn scan's
deconvolution and was identified, folded and logged with it, in every mode, on the instrument too.
The `else deconv_.clearStoredMS2();` is the fix; it is the regular-path twin of the exploration
engine's `last_result_deconvolved_` (ADR-0045 below). **No unit test pins it** (`FLASHIda_ProcessScan_test`
was outside the approved set of the change that found it); the `refine_discriminating_cytc` golden's
`scan_results` MS3 rows at `mass_count 0` are the only witness. Do not remove the `else` as dead code.

### An MS3 sweep calibrates each ion without itself (ADR-0075)

`MS3FragmentMatcher::calibrateAndScore` pools the loose-pass (500 ppm) errors of every variant, each kept with the
theoretical entry it fell on (`MatchDetail::theoretical_index`). Each entry's offset is the median of that pool with
every error on the entry, in every variant, left out; an entry with no error takes the whole pool's median, and an entry
whose errors are the whole pool takes 0 (no other ion to learn from). The tight pass then tests each entry against the
observed mass corrected by its own offset (`matchSpectrum`'s `ppm_offsets`; absent, `matchSpectrum` is unchanged).
Why: a sweep's variants repeat its masses, and a pooled median that holds a mass's own copies anchors on them; on the
2026-11-12 decoy ladders the whole-pool median gave 9.7–26 × λ chance matches per MS3 scan, leaving out the judged
error alone 5.6–22 ×, leaving out the ion 0.4–3.7 ×. The row's `ppm_offset` / `correction_factor` and
`[MS3-SCORE] calibration_ppm` log the whole pool's median. `Refine::pieceTest` calls it with loose 0 (nothing pooled,
every offset 0); its own shared-piece median (`Refine.cpp`) is a separate, single-spectrum calibration and unchanged.

### Any MS2 activation feeds an MS3; only a b/y precursor's MS3 is annotated (ADR-0054)

No activation refuses an MS3: `winnerIsMs3Capable` (`"" | HCD | CID`) is **deleted**, so `exhaustive`
and `discriminating` plan behind an ETD/EThcD/ETciD winner like the other modes, whose MS3s already
went out as `ETD;HCD`. What ETD changes is the targets' **ion type** (c/z), and two predicates in
`MS3FragmentMatcher` keep the two questions apart — **do not merge them**:

- **`isKnownIonClass`** (a b c x y z) — *has a frame*: a span to target and to project through.
  Unchanged. It also keys `Exploration`'s ADR-0023 forced `remaining_precursor` sweep, so narrowing it
  would force every c/z sweep too, which ADR-0054 d.5 rejects (a c/z sweep keeps its metric; under
  `fragment_count` every variant scores 0 and `ce_min` wins, knowingly).
- **`isAnnotatableMs3Precursor`** (b y) — *its sub-ions are predicted correctly*. `getMS3IonTypes`
  reads c/a as b and z/x as y, so a c precursor's yb/ya pieces sit ~17 Da low. Gates the two places
  sub-ions are predicted: `calibrateAndScore` (the regular MS3 return and the `fragment_count`
  re-score both go through it) and `Refine::pieceTest`'s per-side `predict`. A gated MS3 is acquired
  and logged and matches nothing, so the existing `if` skips the identification row, `feedScan` and
  `foldMs3`; the regular return prints `[MS3-UNANNOTATED] precursor_id= ion= scan=` (an event,
  frame-having classes only — `'u'`/`'\0'` stay silent).

`decideRefine` keeps path a **without a recount** when no dispatched record has a b/y side
(`[DISCRIMINATING] … decided=kept reason=no_annotatable_pieces`, no `[MS3-DECIDE]`, no pooled row):
with no MS3 confirmation the recount degenerates to an MS2-only ranking that can flip the winner. The
whole-protein matcher an MS3 *exploration variant* also runs (`computeFragmentMatch_`) never reads the
precursor class and is left alone — its full-protein ions are right whatever the precursor was. Its match
feeds the variant's provisional `fragment_count` score only: at MS3 it is never published as an
identification (`feedResultImpl_` assigns `identification.result` for an MS2 group only), because it carries
no equivalent ion and no adjusted mass. An MS3 variant's `identification.tsv` row is the calibrated
re-score's; a measuring-metric pre-scan and a variant of a dropped sweep have none.
`ETciD` (one spelling, PSI-MS MS:1003182) joins both `needsCollisionEnergy` and `needsReactionTime`.
**No ETD→MS3 has run on hardware**; whether the instrument takes `ReactionTime "7;0"` is open.

### `identification.method` picks the MS2 annotator; Refine is the second one (ADR-0049)

(`characterization.method` and `.refine` until ADR-0061; read as `config_.identification().method` /
`.refine`.)

`FragmentAnalysis::getTopFragmentMatches` is **the one switch**: `AnnotationMethod::Refine` →
`runRefineMatching_`, else `runTagBasedFragmentMatching_` (FLASHTagger + FLASHExtender, the default
`flashtnt`). Every MS2 route reaches it — the regular MS2 through `initiateNextLevel` →
`ProteoformTracker::selectNextLevelTargets` and every exploration variant through
`computeFragmentMatch_` → `ProteoformTracker::identifyExplorationFragments` — so there is no second
place to switch. The two config-unreachable matchers
(`getAmbiguityEnclosingIons`, `getTerminalFragmentIons`) still call the tag path directly.
`runRefineMatching_` touches nothing in FLASHTnT, so the two-files-reach-into-FLASHTnT rule above is
unchanged.

- **`Refine` (`FLASHIda/Refine.{h,cpp}`) is pure static functions** — `ladder`, `detect`, `assemble`,
  `discriminating`, `ptmSites`, `pieceTest`, `recount` — with no engine state, the same shape as
  `KineticModel`. `RefineSettings` lives in its own `RefineSettings.h` **because of an include cycle**,
  not for tidiness: `Refine.h` needs `FragmentAnalysis.h`, which includes `Config.h`, which holds the
  settings by value. Do not fold it back. `RefineSettings::tolerance_ppm` is not authored — `Config`
  copies level 2's tolerance into it after parsing the tolerances, so the two can never disagree.
  `Detect::levels` keeps **every** level with its `Level::excluded` reason (`None | NeutralLoss |
  NeutronSatellite | Capped`, first match wins); `Detect::n_kept` is what is built from and what
  `identification.tsv` logs as `n_levels`.
- **The contract is the tag path's** (decision 6): score = S, `ptm_sites` = the best path's windows,
  `fragments` = its explained ions, `tag_count` −1 (since ADR-0061 push 4c the searched MS2's own tag count is set
  after it), plus one extra member,
  `ProteoformMatch::refine` (`shared_ptr<const RefineCandidates>`: Detect, every candidate path, the
  precursor mass difference). `getTopFragmentMatches` **resets it before the switch**, so a reused
  result never carries a stale record; only the tracker reads it.
- **The tracker owns the round.** At finalize the winner's record becomes `ProteoformModel::refine`
  (`RefineState`); with two or more totals, the masses the two best paths explain differently become
  `DiscriminatingRecord`s, intensity-descending. Under `CharacterizationObjective::Discriminating`,
  `planNextScans` returns `planDiscriminating_` **until the round's decision fires** (budget counts
  species, nominal-mass dedupe, emitted through `emitRawMassTarget_`, which it shares with
  `planExhaustive_`), and **ambiguity** after it — so a scan with only one total runs as ambiguity from
  the start. **A round that can never close through MS3 is decided without a recount or a pooled
  row**: when the pool is exhausted with nothing awaiting its return (every record refused by the
  fragment-charge floor or dispatched in an earlier survey), `planDiscriminating_` sets `decided`,
  prints `[DISCRIMINATING] precursor_id=N decided=fallthrough reason=pool_exhausted`, and
  `planNextScans` plans ambiguity **in the same call** (an ETD winner is no longer a second reason —
  ADR-0054; see the section above); while a dispatched record is still unreturned it waits instead
  (`[MS3-PLAN] … targets=0`). A record only path b explains (`ion_a` empty) is dispatched with the
  **unassigned** label (`'u'`/0, as `planExhaustive_` marks a mass with no identity), never under
  path b's label; and a record is frozen once the decision has fired — `recordDiscriminatingScan`
  returns false and a late MS3 is folded by the ambiguity path. Under `refine`,
  `narrowModifications_` runs `shift_only`: only a fragment that carries the shift moves a boundary.
- **The MS3 return path closes the round** (`FLASHIda.cpp`, regular MS3 branch). Each returning MS3
  feeds `recordDiscriminatingScan` — whether or not the winner-context match found anything, and with
  an empty spectrum when the scan deconvolved to nothing, because the round closes only when **every
  dispatched record has returned**. Then `refineRoundClose_` calls `decideRefine` (recount over
  **every** record, dispatched or not: `Refine::recount` throws `InvalidParameter` when the two paths'
  shared counts disagree, which a partial target list produces; a flip re-seeds the modifications and
  fragments from path b; the pooled row's trigger is `MS3-DECIDE`) and re-dispatches through
  `Exploration::dispatchFromModel` — the tail of `initiateNextLevel`, extracted byte-identical —
  replaying the **winner MS2's own command** (`ProteoformModel::winner_cmd`), never `ms2_ctx`, which
  is the *first* fed scan's. The close runs **after** the context-cache erase: the re-dispatch inserts
  into `ms2_context_cache_` and may rehash.
- **`Config::validate` throws three ways**: a `refine` block under `flashtnt`; `discriminating` under
  `flashtnt`; and `discriminating` with level-3 exploration enabled — an exploration MS3 returns on the
  exploration branch, never reaches `recordDiscriminatingScan`, and the round could never complete.
- **Markers, events only** (`events.log`; none in `ida.log`, ADR-0032): `[REFINE]` (one attribution line naming the
  precursor and scan, because `[LEVELS]`/`[CANDIDATES]` print `scan=0` — the stored MS2 carries no
  scan number), `[LEVELS]`, `[CANDIDATES]` (one per total), `[DISCRIMINATING]`, `[MS3-DECIDE]`,
  `[MS3-REDISPATCH]`.

Accepted limits, each deliberate for now: MS3 evidence folded during a round is **dropped on a flip**
(the fragments are rebuilt from the winner scan as path b explains it; the round's staged scans are
gone); `Refine::assemble` returns an empty result **silently** on a prefix/suffix ladder size
mismatch; a re-matched non-winner (`getRematchedNonWinnerMatch`) is **stale after a flip**;
`[REFINE]` is **not printed on the exploration MS2 path**; and the golden's empty-fixture MS3 rows
carry rt 0 (cosmetic — the zero-peak fixture has no retention time).

### A sweep is exactly its grid, and every pre-scan carries its precursor share (ADR-0048)

**Sweep order.** `exploration.order` (per level, default `"ascending"`) decides only *when* the grid is acquired. Under
`"random"` `initiate` permutes the whole variant list with `Exploration::sweepPermutation(n, kSweepSeed ^ group_id)`
(splitmix64 Fisher–Yates, pinned to a literal in `FLASHIda_exploration_test`), so each triggered sweep has its own
reproducible order across every activation and both axes. `variant_index` stays the acquisition position;
`ExplorationVariant::grid_position` keeps the grid position, and the winner loop breaks a score tie on it, so the
order never changes which pre-scan wins. Nothing is logged that was not before.

`exploration.activations` is a **list**, and `buildVariants_` emits each activation's whole sweep as a
contiguous block before starting the next — and **nothing else**. The per-activation CE-0 baseline
of ADR-0029 is gone (superseded): no reference scan is synthesized under any metric, so
`ExplorationVariant::is_baseline`, `ExplorationGroup::baseline_intensity`, the suppression-when-the-
grid-contains-the-turn-off-point logic and the "empty reference bars its activation" arm no longer
exist, and `variant_index` is one running counter from 0. The grid itself is **index arithmetic with
the top edge appended** (ADR-0047 d.8): `n = floor((max−min)/step + 1e-9) + 1` points from `min`, then
`max` if the last point fell short — so `15..40/2` ends at 40, not 39 (the old `ce += step` loop
dropped the edge; `etc/method.json`'s MS3 grid is exactly that shape).

What every variant contributes instead is its **precursor share** (`ExplorationVariant::precursor_share`,
computed by `precursorShare_` through `KineticModel::precursorShare`): Σ raw intensity inside the
anchor's isolation window over Σ raw intensity of the **same scan**, −1 when the scan carried nothing.
It is read for every variant under every metric, on trap and Orbitrap alike, and it is what
`scan_results.precursor_share` logs. `remaining_precursor` scores `max(0, 1 − |share − target|)`
against `remaining_precursor_target` — the key kept its name, its value is now a share — and a
share of −1 scores **`-1.0`**, "not scored", so an unmeasured variant can never win while a fully
depleted one (share 0) competes like any other.

⚠️ **The reaction-time floor survives the baseline it was written for.** A reaction time of 0 is
**rejected** by the instrument, so `Config::validate` still throws on a `reaction_time_min` below
`Config.h`'s `MIN_REACTION_TIME_MS` (0.03 ms) when an ETD-family activation is swept. Gate it on
`needsReactionTime`, never unconditionally: every CE-only config leaves `reaction_time_min` at its 0
default. `ScanCommandQueue::cancelByScanIds` keeps having no production caller: a newer survey drops an
unsent sweep through `replacePlan` and forgets it with `Exploration::dropGroup` (ADR-0060).

### The Kinetic Model decides by the modeled maximum (ADR-0047)

`ExplorationMetric::KineticModel` (`"kinetic_model"`, logged as `4`) is the fourth metric and the odd
one: **no variant wins.** `computeExplorationScore_`'s arm returns `-1.0` for every pre-scan (the
share is still read and logged), and at group completion `feedResultImpl_` fits the Gompertz curve
`P(CE) = P₀·exp(−exp((CE − c)/w))` to the (CE, share) points (`KineticModel::fitGompertz`), runs the
fragment-pool recursion `B(k) = e^{−k₂}B(k−1) + s_k(1−e^{−k₂})/k₂` on the fitted rate `s = −dP/dCE`
over 1..`ce_max` with the level's authored `kinetic_k2[activation]`, and takes the integer argmax
(`KineticModel::modeledMaximum`). That energy is acquired as the follow-up **at every level, overrides
or not** — it is computed, so there is never a pre-scan to cascade from; the follow-up is rebuilt from
variant 0 (`ref_idx`), whose stage-0 context, fragment scores and `window_snr` the whole group shares.
`winner_index` stays −1 and `winner_tracking_id` empty. A fit with fewer than `kMinPoints` (3) finite
shares, or a non-converging lattice, ends the group with `[EXPL-ABORT] reason=no-fit` and **no
follow-up** — there is no fallback to the authored energy (owner decision, d.7). An edge maximum
(`ce == 1` or `ce == ce_max`) is still acquired and flagged. `[EXPL-MODELED]` prints the whole fit;
`FeedResultInfo::KineticInfo` carries it to the `scan_results` row (`fit_p0 … energy_source`).

⚠️ **The fit is the study's lattice search, ported verbatim — not Levenberg-Marquardt.** ADR-0047
d.3 said LM; implementation measured LM against the vendored pipeline on 1080 fits and got the same
CE* on 941 (within ±1 on 1035) — the Gompertz RSS has a flat valley in (P₀, w), and a 5-start LM made
it *worse* (921). The lattice port agrees on 1080/1080 with identical mean regret, so it is what ships
(`KineticModel.cpp` documents the lattice: c on 0..60/0.5, w on geomspace(0.25, 25, 41), closed-form
P₀ clipped to 1.1× the highest share, a ±0.5/×1.3 refinement). `KineticModel_test` pins that
equivalence against `FlashIDA/test-data/kinetic/kinetic_expected.tsv`, which is generated by the
vendored pipeline in `FlashIDA/test-scripts/kinetic/` — change one without the other and the test
says so.

⚠️ **`k2` is never defaulted, and the kinetic sweep is one CE-swept activation.** `Config` throws
when the metric is `kinetic_model` and the swept activation has no `kinetic_model.k2` entry, when the
object is authored under any other metric (an authored value nobody reads is the ADR-0040 stomped-key
shape), when more than one activation is swept (one k₂ feeds one pool), and when the swept activation
is ETD-family (nothing on the CE axis to fit). An entry for an unswept activation is a `[CONFIG-WARN]`.
The trap check admits both share-based metrics; the multiplexed rejection covers both too. The only
measured k₂ is HCD's (LOPO 0.0364–0.0390, all-in 0.0388); the MS3 golden's CID value is a
**placeholder**, and no kinetic sweep has yet run on hardware.

### The annotator metrics: raw_score, bond_coverage, explained_signal (ADR-0074)

Three `ExplorationMetric` values, appended (`RawScore` 5, `BondCoverage` 6, `ExplainedSignal` 7; the number is the logged
`exploration_metric`), read from `FLASHAnnotatorAlgorithm`'s **terminal** fragments of one PrSM in the pre-scan's raw
spectrum. `FLASHIda/AnnotatorMetrics` holds the pure functions; `isAnnotatorMetric()` is the predicate. They are reading
metrics.

- **MS2.** A variant is scored on its own best PrSM. `raw_score` is the hit's `RawScore` (before FLASHTnT's blind-shift
  penalty, so the decision value is the logged `raw_score`) and `bond_coverage` its `RawCoverage / 100`; both are copied in
  `FragmentAnalysis::matchSearched_`, for every searched MS2, sweep or not. `explained_signal` needs envelopes FLASHTnT does
  not expose, so under that metric alone the variant's search job (`FLASHIda.cpp`, no lock held) calls
  `AnnotatorMetrics::annotateSearched` with a per-thread annotator (`t_annotators`) and stores the value on
  `SearchedPrSM::explained_signal`. FLASHTnT annotates all PrSMs of a spectrum in one call and pools the chance rate over
  them, so a second annotation of the best PrSM alone need not reproduce `RawScore`: the engine never re-derives it at MS2.
- **MS3.** `computeExplorationScore_` returns `-1.0` at arrival. When the group completes, the `FragmentCount` batch runs
  unchanged for the variants' `identification_result` (and the extra identification rows), then every received, non-trap
  variant's raw peaks (`ExplorationVariant::result.getOriginalSpectrum()`) are annotated against
  `AnnotatorMetrics::ms3TargetPrSM`: the target's subsequence with the proteoform's modifications inside it, fragment charges
  capped at the target's charge, the ion types of the variant's own MS3 activation, tolerance `level(3).tolerance_ppm`. A y
  target is a truncated chain; a b target carries `−18.010565` Da on its last residue, which the annotator adds to every
  suffix ion and no prefix ion (the MS3 matcher's `yb`; its `ya` has no annotator counterpart and is not counted). All three
  values land on the variant's `identification_result`; the configured one becomes its score.
- **Unscored is `-1.0`, never 0**: a dropped search, a variant that identified nothing, a c/z/a/x MS3 target
  (`isAnnotatableMs3Precursor`), an MS3 group with no identified proteoform. A group with no scored variant ends
  `[EXPL-ABORT] reason=no-winner`. (`fragment_count` scores such a sweep 0 and crowns the first grid point; that is unchanged.)
- **Load**: each needs an `identification` section; an `IonTrap` override is refused, as under a counting metric.
- A variant with no matcher fragment writes no identification row, so its three values are logged only when it matched.

### Reading vs measuring exploration metrics — the MS3 evidence rule (ADR-0020)

**`FragmentCount` and the annotator metrics (above) identify their own MS3 pre-scans.** For `FragmentCount`: It scores a variant *by* matching it, so
`ProteoformTracker::scoreCalibratedVariants` — the only matcher that works on an MS3 spectrum, since
it is seeded with the proteoform context and the parent ion index — runs inside a
`metric == FragmentCount` guard in `Exploration.cpp`. `MassCount`, `RemainingPrecursor` and
`KineticModel` are **measuring**: they score from peak counts and the precursor share and never invoke
a matcher. `isMeasuringMetric()` (`Config.h`, beside the enum) is the classification; do not open-code
`== FragmentCount`.

Three things follow, and the third is the trap:

- `ExplorationVariant::identification_result` has **one write site**, inside that guard. Under a
  measuring metric it is empty for every variant.
- The inline fold at MS3 group completion is gated on that field, so under a measuring metric it
  cannot fire — which is why a completed MS3 sweep used to contribute *nothing*: no
  `identification.tsv` row, no `foldMs3`, no pooled evidence, after paying for every pre-scan.
- Therefore a **measuring metric at MS3 always emits a production scan**, regardless of `overrides`.
  The post-decision dispatch fires on
  `kinetic || !overrides.empty() || (msn_level >= 3 && isMeasuringMetric(...))` — the first disjunct
  is ADR-0047 d.5: the Kinetic Model's energy is computed, so there is never a pre-scan to cascade
  from, at **either** level. The follow-up returns on the *regular* MSn path and is matched there.

`overrides` means *"the pre-scans ran at degraded settings"* — variants are built from
`scans[0] + overrides`, the production scan from `scans[0]` alone. That is the whole reason the
original gate keyed on it. **MS2 is exempt** from the metric condition: `computeFragmentMatch_`'s
whole-protein matcher is correct at MS2, so every MS2 variant is identified under every metric and an
MS2 group cascades instead of re-acquiring — **with one exception that belongs to the analyzer, not
the metric** (ADR-0045, below).

### A trap pre-scan is measured, never identified (ADR-0045)

A variant whose **own `ScanCommand`** names the ion trap (`isTrapAnalyzer(cmd.analyzer)`, read in
`Exploration::feedResult` — the engine's **first** read of `ScanCommand::analyzer`; every other site
is a write) is neither deconvolved nor matched nor fed to the ProteoformTracker, at any level, under
any metric. The only thing read from it is the raw precursor-window sum behind the remaining-precursor
ratio. Unit-resolution trap peaks have no business in a ppm-tolerance identification — and until
ADR-0044 withdrew the ADR-0026 scan-range binding this was inert by accident, because a ~2 Th window
deconvolves to nothing. Three mechanics, none of them optional:

- `ExplorationVariant::measured_only` is set from the command, never from the returning scan
  (ADR-0042's rule for every decision gate); `feedScan` and the FragmentCount batch re-score both
  honour it, the latter with a `nullptr` rather than a scored-0 spectrum.
- **`exploration_deconv_->storedMS2()` is stale after a skipped deconvolution.** `FLASHIda.cpp` reads
  `explorationDeconvMassCount()` / `explorationDeconvSpectrum()` for the `scan_results` row, so without
  `last_result_deconvolved_` a trap variant would log the *previous* variant's masses. Pinned by
  `a_trap_pre_scan_is_measured_never_identified`.
- **Only `remaining_precursor` may score a trap sweep.** `Config::validate` refuses an `IonTrap`
  override under `mass_count`/`fragment_count`: with nothing deconvolved they score every variant 0,
  and winner selection (seeded `-1.0`, strictly greater) crowns `ce_min` every time, silently. The
  `-1.0` arm in `feedResultImpl_` for that pairing is therefore unreachable from config and exists to
  fail closed. The ADR-0023 forced sweep is compatible by construction — it forces the one metric that
  can read such a scan — but the load rule reads the *configured* metric, so a `fragment_count` +
  `IonTrap` exhaustive config is refused even though its unassigned masses would have been fine.

A trap sweep always carries overrides, so it is closed by the production follow-up through ADR-0020's
gate #1; that scan carries all of the sweep's evidence. `scan_commands.tsv` has **no analyzer column**
— the analyzer on both sides is pinned by `production_scan_is_built_from_the_unoverridden_config`
(ctest), not by any golden.

`planNextScans` emits `[MS3-PLAN] … reason=…` naming which of its zero-target causes fired
(`all_mods_localized` is the common, *correct* one: ambiguity mode with everything already localized).
an event — like every engine marker it is in the run folder's `events.log` (ADR-0065), not on the console.

## Config (`FLASHIda/Config.cpp`)

Unknown keys are hard-rejected by one free function `rejectUnknownKeys(obj, allowed, path)` called
~18 times, each with its own hand-written allowed-set, throwing `std::invalid_argument` naming the
offenders and pointing at `FlashIDA/test-data/config_schema_reference.json`. Two exemptions:
`exploration.overrides` (a dynamic string→string map) and a top-level `ms3`, which gets a dedicated
**migration** error before the root check runs.

> **Resolved:** `selection_strategy` is deleted, and with it the one `std::runtime_error` in the
> file. A *present* `selection_strategy` now throws `std::invalid_argument` with a migration message
> naming all seven destinations, so every config error in this file is now the same exception type.

⚠ **`scheduling.target_depth` is listed as allowed and then never read, and that is correct**
(ADR-0033). It sizes the *instrument's* custom-scan queue, which only the host can do — `Flash.cs`
is the sole submitter and the engine has no visibility into how many commands the instrument holds.
It appears in the allowed-set purely because the schema is strict on both sides (ADR-0007) and
`config_schema_reference.json` is generated from the C# model, so a host-only key C++ rejected would
fail `ConfigSchemaParity_test` rather than any real config. **Do not "wire it up"** — there is
nothing engine-side for it to control. It is the only key in the file with this shape; if a second
one appears, that is the point to introduce a marker rather than let the pattern go unremarked.

Traps worth internalizing before touching config code:

- **`.value(key, default)` fallbacks are dead in production.** C# `ToCppJson` emits every key
  unconditionally, so the C++ literals never fire and several disagree with the effective C#
  defaults. Don't read a C++ default and believe it.
- **`kScanKeys` is one lenient 17-key union validating every scan object, and one parser reads all
  17 at all five sites.** Validation and parsing share a single source of truth (`Config.cpp:70-73`
  records why: when the two lists drifted, nine keys passed validation and were then discarded —
  notably `reaction_time` on a follow-up scan, which made an ETD follow-up unconfigurable). The
  set is exposed as `Config::scanKeys()`.
  **The constraint has never been on this side — it is the C# emit set.** A key C# does not emit
  cannot reach C++ regardless of what `kScanKeys` admits, which is how `rf_lens`, `source_cid`,
  `source_cid_scaling` and `scan_rate` were unreachable from `method.json` for every MSn scan while
  C++ parsed, copied and sent them (ADR-0011). `ConfigSchemaParity_test::GeneratedReference_CarriesEveryScanKey`
  now compares the C#-generated reference against `scanKeys()` so the two cannot drift again.
  Level asymmetry that remains and is deliberate: `ms_settings.ms1` omits the five stage-carried
  keys (activation/collision_energy/reaction_time/reagent_*) because `makeMS1` emits
  `num_stages = 0`, so they could never reach an MS1 scan.
- **`ScanConfig.analyzer`'s default flips by parse site**: the in-class default is `"Orbitrap"`, but
  the `ms_settings.ms{1,2,3}` parsers override it with `""` (strncpy'd straight into the
  `ScanCommand`). No parse site keeps `"Orbitrap"` any more: the last one, `tagging.follow_up_scan`'s
  dead inline-object branch, went with ADR-0061.
  **The value is a closed set** (ADR-0045 decision 3): `"Orbitrap" | "IonTrap" | ""`, checked on the
  **raw JSON** in `parseScanConfig` and in the overrides loop — `str()`'s `is_string()` guard would
  otherwise coerce `"analyzer": 5` to `""` and load clean. `null` stays a default because the
  generated `config_schema_reference.json` carries nulls. `isKnownAnalyzer` / `isTrapAnalyzer` live at
  namespace scope beside `needsReactionTime` and are pinned as an exact set by
  `Config_SchemaProjection_test::analyzer_is_a_closed_set`. An **unset** analyzer is not the trap: a
  pre-scan that reaches the trap through the instrument method's default is invisible to the rule
  (recorded gap).
- **Selection is never authored per level any more** — it is projected from `identification.characterization.mode`
  and `precursor_selection.rank_by`, so there is no "level present without a selection key" case
  left. `Config::level(n)` still returns a static `default_level_` with `selection = None` for a
  level not in `levels_`, but `levels_` is now always `{1,2,3}`. `MSLevelConfig::selection`'s
  in-class default is also `None` (it was `Intensity`, which meant merely *defining* an
  `ms_settings.msN` block switched that level on, because that parse materialises the level before
  any selection is read).
- **`applyOverrides` silently ignores unknown keys** (hand-written if/else chain over 17 fields, no
  `else`), and override **values must be JSON strings** — `"collision_energy": 30` throws an
  nlohmann `type_error`; `"30"` works. The base is the level's `scans[0]`, which is why `validate()`
  requires exactly one scan config at an exploring level.
- `deconvolution.tol` must always carry ≥ 3 entries — C++ throws below `max_level`, and
  `Config` now materializes levels {1,2,3} unconditionally (it used to hold only because
  `selection_strategy` was required and every config named all three) — `toleranceList()` walks
  `levels_` positionally, so a missing level would shift every `tols_[ms_level-1]` index.

## Logging (`FLASHIda/IdaLogger.cpp`)

Five append-only streams under **one folder**, `runtime.log_dir`, joined with the fixed basenames
below. All five open together or none do — an empty `log_dir` opens nothing, which is how a fixture
with no `runtime` section (or `"runtime": {}`) logs nothing. `IdaLogger` never creates the
directory: the host does, and a missing one leaves every stream silently closed. Note the authored
`method.json` layer gives `""` the *opposite* meaning (`.` = the working directory) — ADR-0015.
`IdaLogger` is **lock-agnostic** — locking is
the caller's responsibility.

| Stream | Cols | Role |
|---|---|---|
| `ida.log` | — | free-text MS1 summary (not a TSV) — the record of acquisition *decisions*, and the only stream keyed on the **instrument** scan number. It **no longer has a reader outside FLASHIda** (ADR-0046); `PrecursorSelection`'s `target_log_files` loader is the one that remains |
| `scan_commands.tsv` | 38 | one row per **dequeued** command, so a re-planned-away command has none; the wide MS3-fragment stream (32→34 ADR-0026 `first_mass`/`last_mass`, between `faims_enabled` and the trailing `enqueue_ts`; →36 ADR-0060 `decision_survey_id`/`isolation_purity` in the same gap — the survey whose ranking last chose the window and the purity read there, both empty on follow-ups, MS3 and MS1 rows; a bought `R` carries its `Q`'s survey and an empty purity unless a newer survey adjusted it; →37 ADR-0069 `predicted_bond_yield` in the same gap, the predicted bond yield an MS3's yield gate admitted it at, −1 on every row no yield gate tested; →38 ADR-0072 `max_it` in the same gap, the injection time the command was sent with, `0.0000` when its scan config sets none). **The stream FLASHDeconv reads** (`-FD:scan_commands`, ADR-0046) — see *Offline coupling* below |
| `scan_results.tsv` | 56 | pure acquisition-**event** log per `processScan` (34→29 identification payload moved out, →28 per-charge deconv restructure, →29 `deconv_qscores`; →32 the identification-YIELD block `tag_count`/`fragment_count`/`tic_coverage` after column 12; →36 ADR-0038's quant block `quant_channels`/`quant_condition_means`/`quant_fold_change`/`quant_verdict`, **appended** after `dequeue_ts` because every `r.colIndex` pinned in `FLASHIda_LoggingFields_test` is ≤ 21 and appending is the only placement that invalidates none of them; →41 ADR-0040's `consensus_*` block, renamed `group_*` by ADR-0056 (→53 with `quant_p_value`, `group_p_value`, `group_channels`; no vote — one test on the group's share-normalised aggregate); →50 ADR-0047's kinetic block `fit_p0`/`fit_ce_max_s`/`fit_width`/`k2`/`modeled_ce`/`modeled_b`/`fit_points`/`fit_rss`/`energy_source`, written on the row that completes a `kinetic_model` group and sentinel-valued everywhere else, `energy_source` alone also `measured`/`none` on every other metric's completing row. Column 12 is **`precursor_share`** (ADR-0048, renamed from `remaining_ratio` because the quantity changed): the anchor window's share of the scan's raw signal on EVERY MS2/MS3 row — pre-scan, follow-up, regular — read in the window the command isolated last, −1 on MS1 rows; →56 ADR-0061 push 4b's `search`/`search_wait_ms`/`search_ms`, the search pool's verdict and timing, the two timings masked in the goldens) |
| `identification.tsv` | 52 | per-scan MS2/MS3 identification leaf (32→34: `tag_count` beside `flash_extender_score`, `fragment_qscores` inside the aligned fragment-mass table; →41 ADR-0049's Refine block `annotation_score`/`total_shift`/`second_total_shift`/`second_score`/`cut`/`n_levels`/`precursor_mass_difference`, **appended**, valued on the regular AND the exploration self-identification MS2 rows under `refine` and sentinel-valued everywhere else — `-1`, empty, empty, `-1`, `-1`, `-1`, empty. A real score of −1 is told from the sentinel by a non-empty `total_shift`; `flash_extender_score` and `tag_count` stay −1 under `refine`; →46 ADR-0061 push 4b's `protein`/`raw_score`/`identification_score`/`closed_out`/`inferred_precursor_mass` after that block on every row, only `protein` valued until push 4c; →47 ADR-0064's `quant_ms3`; →49 ADR-0074's `bond_coverage`/`explained_signal`, appended last, `-1.0000` where the row has no annotated PrSM; 2026-11-12 `identification_score` and `flash_extender_score` swap places, the Identification score next to `proteoform` (col 8) and the extension score after `raw_score` (col 44); →52 ADR-0075's `chance_p`/`chance_e`/`confirmation_p`, appended last, written from their log10 in scientific notation, `-1` where the row carries none — see *The chance-match p-value* below) |
| `pooled_identification.tsv` | 20 | per-precursor cumulative proteoform trajectory (→20 ADR-0061 push 4b's `protein`, the accession the model's winner chose) |

- **An MS3 row's `ms2_fragments` / `ms2_fragment_masses` are in the Proteoform frame** (`CONTEXT.md`): the
  equivalent ion is labelled and weighed over the proteoform's own residues, `[region_start, region_end)`,
  as an MS2 row of that proteoform would. `FragmentMatch::equiv_frame_index` is the logged label, set by
  `calibrateAndScore`; `equiv_index` (counted over the reference protein) is the pooling key and is never
  written. A b sub-fragment of a b precursor therefore reads the same label and mass in `ms3_fragments` and
  `ms2_fragments`, on any region. The pooled model finds an MS3 fragment's peak by its measured mass and
  pools the adjusted one.
- **Whether an MS3 fragment's equivalent ion carries an ambiguous shift it straddles is one rule**,
  `FragmentAnalysis::equivalentCarriesShift`, read by the logged mass (`coveredAmbiguousInEquivFrame`), the
  identification row's localization (`narrowFragmentPTMSites`) and the pooled localization
  (`narrowModifications_` Pass B). Same direction: the sub-ion's own verdict. Complement: the opposite, since the
  precursor holds the shift. An N-terminal net loss (negative shift anchored at residue 1: Met-excision +
  N-alpha-acetyl) stays on the prefix equivalent of a complement ion whichever way the complement matched.
  Under `identification.method: refine` the pooled range does not narrow on a net loss a complement ion only
  implies (the complement matched without it): no ion measured the shift (ADR-0049).
- **`ida.log`'s two identifiers are different quantities and neither substitutes for the other**
  (ADR-0035). `MS1 Scan#` is the **instrument's** scan number, arriving as a trailing `int` on
  `processScan` and used for nothing else; `Access ID` on the same line is the base-94 **tracking
  id**, which is the join key to the other four streams. The port had `MS1 Scan#` carrying the
  tracking id, which made FLASHDeconv's old `ida.log` join unsatisfiable rather than merely
  mis-keyed — tracking ids count 1, 2, 3… and never reach the thousands instrument scan numbers do.
  (That consumer — `precursor_map_for_ida_` / `findPrecursorPeakGroupsFormIdaLog_` — **is gone**:
  ADR-0046 joins by tracking id through `scan_commands.tsv` instead.) `<= 0` means "not supplied"
  and the writer falls back to the tracking id, warning **once** per run (it flushes under
  `ida_log_mutex_`).
  ⚠️ **Both readers of this grammar anchor on `" 0 targets"` with the leading space.** The bare
  substring also matches `"10 targets"` / `"20 targets"`, which command fan-out makes routine.
  `parseFLASHIdaLog` is one reader — since ADR-0046 it has **no production caller**, only the four
  writer round-trip sections in the logging tests, and it stays where it is by decision;
  `PrecursorSelection`'s `targeting.target_log_files` loader is the other, and it is in-scope —
  feeding an engine-written `ida.log` back in is a supported round trip.
- **Every `Mass=` line ends in `\tID=<tracking id>`, and a re-plan writes `Superseded=`** (ADR-0060).
  `ida.log` stays a decision-time record: a survey whose `replacePlan` dropped commands lists their ids on
  one `Superseded=<id> <id> …` line right after its MS1 header, and the `target_logs` loader collects every
  such id in a first pass and skips the `Mass=` lines carrying one, so a superseded decision is not an
  extra observation in `in_depth`'s `1 − ∏(1 − qscore)`. A log without `ID=` loads as before.
  `scan_results`' MS1 `child_ids` likewise stays decision-time and still lists a superseded id.
- **Every TSV cell that carries a tracking id is written through `quoteLogCell`** (`FLASHIda/LogCell.h`,
  FDmerge `96f8ab3732`). Base-94 ids use 0x21–0x7E, so ids 8,836–17,671 begin with `"`, which a
  CSV-style reader (pandas, Python `csv`) takes as an opening quote. A cell beginning with `"` is
  RFC 4180-quoted, and no other cell changes, so every log of a run under 8,836 commands is
  byte-identical. This covers the ids, the space-joined id lists and `scan_description`. ADR-0060's
  `decision_survey_id` is quoted on write too, and `ScanCommandJoin::parse` reads it back through
  `unquoteLogId`. `ida.log` is not a TSV and is not quoted.
  `NormalizeIdaLog` masks both id kinds, so no golden sees their values.
- **An `ida.log` target is a MASS, active within `rt_window` of a reference acquisition** (FDmerge
  `20db11cef3`, "GTM"). Under `targeting: inclusion` with `files.target_logs` and no inclusion TSV
  rows, `PrecursorSelection` makes each logged mass a target only while the current RT is within
  `precursor_selection.rt_window` of an RT at which a reference acquired it, and matches it at any
  charge. FDmerge `1fe5ee0b3e` keeps the two previous behaviours as keys, both exposed in `method.json`:
  - `log_target_rt_window: false` makes every log target a target at every RT. That was the behaviour
    from 2025-07-07 to 2026-10-21.
  - `log_target_charges: "reference"` matches only at a charge the reference acquired. Each acquisition
    uses that charge up, through `target_mass_charge_map_`, and the erase is journaled as a
    `LegacyTargetChargeErase` stamp, so a re-planned-away command spends nothing.

  The loader files `AllMass` lines under their own survey. The gate is C++: `FLASHIda_ProcessScan_test`'s
  GTM sections (FASTA tag-expanded targets, which used to add to these, are retired with ADR-0061). The `log_target_charges: "reference"` path is also exercised by `FLASHIda_ChargeModes_test`'s
  two legacy-inclusion cases.
- **`Mass=` and `AllMass=` are the same rendering, and that is asserted** (ADR-0035 decision 5).
  Both are `std::fixed << setprecision(4)`, so a target's mass appears **verbatim** among its own
  entry's `AllMass=` tokens. `Mass=` used to say `std::defaultfloat` with no precision of its own —
  and every number on these lines goes through ONE ostream whose flags and precision are **sticky**,
  so it inherited the header's `setprecision(4)` on an entry's first target and the previous line's
  `Features` `setprecision(6)` on the rest. Under `defaultfloat` that counts *significant digits*:
  `Mass=1.235e+04` and `Mass=12383.3` against an `AllMass=` of `12351.3933` / `12383.3180`, up to
  ~5 Da wrong, on 99.4 % of targets in the field. A new field added to these lines **must set its
  own precision**; `ida_log_target_line_precision_is_pinned` pins every token's rendered width so
  the next one cannot repeat it silently.
- **`scan_commands.tsv`'s `mono_mass` had the same defect, and has its own formatter now**
  (ADR-0046 decision 10). It shared the `sc()` lambda with ten other columns, i.e. the stream
  default of six *significant* digits — `12351.4` for a commanded 12351.3933. `mass()` prints four
  decimals and keeps the stage-less rows' literal `"0"` (`num_stages == 0`), because `"0.0000"` would
  revalue every MS1/AGC row and break `commands_ms1_agc_stageless`'s raw `== "0"`.
  ⚠️ **No golden can see this.** `GoldenNumericComparer` accepts `12351.4` against `12351.3933`
  with ~1800× headroom (`RelTol 1e-3`), so the 28 `scan_commands` goldens were deliberately **not**
  recaptured and still show six significant digits. The only gate is
  `FLASHIda_Logging_test::scan_commands_mono_mass_is_written_at_four_decimals`, which compares the
  written cell **byte for byte** with the dequeued command's own `double`.
- **`ChargeRange` is the species' measured charge envelope**, `PeakGroup::getAbsChargeRange()` of
  the PeakGroup the command was built from (ADR-0035 decision 6) — not the trigger charge printed
  twice, which is what it was. `writeIDALogEntry` takes it through `ms2_sources`, a vector
  **index-parallel** to `ms2_commands` filled at push time in `FLASHIda::processScan`, rather than
  looking the mass up in the deconvolved spectrum: several PeakGroups routinely share one mass in a
  survey (ADR-0036 split envelopes — 48 of the 1324 committed golden target lines sit on a species
  with 2–4 of them), each carrying a different charge subset, so a lookup would pick one of several
  answers with nothing to notice. A degenerate `[z-z]` told every reader of the log that the
  species was seen at one charge only.
  ⚠️ This file used to say the range "reaches `FLASHDeconvFeatureFile` as columns 10–11 of
  `*_ms2.feature`". That was wrong twice: the file was **`_ms1.feature`** (the synthetic one-scan
  feature minted for a precursor rebuilt from a log row), which a single-run TopPIC search never
  opens — and since ADR-0046 no `ida.log` value reaches FLASHDeconv at all. It is also why
  `scan_commands.tsv` deliberately carries **no** charge-range column.
- **Rows are written at dequeue**, from 2 sites inside `getNextScanCommand`. So row order == dequeue
  order, and an enqueued-but-never-dequeued command never appears. Only the priority-dequeue site
  passes `precursor_id` and `ms3_proteoform` (a one-shot `takeMS3Proteoform`); the scheduled-prescan
  site (Step 1) logs `precursor_id=0` and an empty proteoform. There used to be three sites, the
  third being Step 5's fabricated prescan — the idle survey now goes out through the priority-dequeue
  site like any other command (ADR-0031).
  ⚠️ Neither site takes `analysis_mutex_` — that is ADR-0025, and it is why `IdaLogger` serialises
  each stream itself.
- **Two-stage MS3 scalars**: 11 columns are emitted as `"stage0;stage1"` **only when
  `msn_level == 3`** (mono_mass, qscore, charge_cos, charge_snr, iso_cos, snr, charge_score,
  hcd_energy, ppm_error, precursor_intensity, peakgroup_intensity). Everything else emits stage 0
  alone, byte-identical to the pre-two-stage format.
- **Stage-less rows emit literal placeholders** — for `num_stages == 0` (MS1/AGC) the per-stage
  columns are forced to `"0"` (`"none"` for activation) so a tab-splitting parser doesn't drop
  trailing empty fields and tokenize the row one column short. The test TSV parser carries the
  mirror-image guard (hand-splits to N+1 fields for N tabs).
- **Co-isolation notches add a SECOND delimiter inside the per-stage columns** (ADR-0016): `charge`,
  `precursor_mz` and `isolation_width` carry `','`-joined notches *within* each `';'` stage group, e.g.
  `"17,16,15;4,5"`. It mirrors the wire grammar, so the row records what was commanded rather than
  only the anchor window. Collision energy and activation do **not** gain the axis — all notches of a
  stage fire into one fragmentation event. With `precursor_charges`/`fragment_charges` at their
  `single` default no notch exists, so these columns stay byte-identical.
- **`scan_results`'s deconv block reuses both delimiters with different meanings** — there `';'`
  separates masses and `','` separates one mass's observed charges (`deconv_charges`,
  `deconv_intensities`); nothing to do with notches. `deconv_qscores`, between `deconv_masses` and
  `deconv_charges`, is `';'`-only: `PeakGroup::getQscore()` is one value per mass — the representative
  charge's, not an envelope aggregate — so it stays index-aligned 1:1 with `deconv_masses`, on every
  MS level.


### The chance-match p-value (ADR-0075)

Reported only, never decided on: the close-out still reads the Identification score. `FLASHIda/ChanceMatch.h` holds the
pure functions, every value a log10 (`log10PoissonTail`, `chanceLambda`, `chainMass`, `formatLog10Scientific`). Were a
scan's N masses unrelated to a chain of mass M, each would land within tolerance of one of its theoretical masses t with
probability Σ 2·tol·t / M, so λ = N · Σ 2·tol·t / M and p = P(X ≥ k). Three fill sites:

- **Own-search MS2** (`FragmentAnalysis::matchSearched_`): k and λ from `matchFromPrSM` on the search's PrSM, under every
  method (under Refine one extra, event-free `matchFromPrSM` into a scratch match). E = p × `SearchedPrSM::reference_size`
  (the proteins searched) × `shifted_windows`, which the search worker sets, with no lock held, to
  `FLASHTnTAlgorithm::shiftedCandidateCount` for a PrSM with a blind shift (`blindModCount_` > 0): the count FLASHTnT's own
  raw-score penalty uses, with its arguments (the reference protein, the hit's `Mass`, the worker's
  `maxTotalModificationMass()`, the prepared database's fixed modifications). That function and `blindModCount_`'s
  declaration are public by an owner-approved locked edit inside the FLASHTnT boundary; nothing else there changed.
- **Winner-re-matched MS2** (`ProteoformTracker::mapNonWinnerMs2_`): `confirmation_log10_p`, k the re-matched fragments,
  λ over every distinct candidate theoretical of the winner ladder (`candidatesFor`), at level 2's tolerance.
- **MS3** (`MS3FragmentMatcher::calibrateAndScore`): `confirmation_log10_p` per variant, k its tight-pass matches, λ over
  the sub-fragment ladder at the tight tolerance, M the fragment's own mass.

The board's `E` column and the summary's best-score line show the E-value of the identification a row's score is from.
Measured on the 2026-11-12 direct-infusion runs, the uniform λ is right to conservative at MS2 (chance matches on shifted
ladders 0.6–0.9 × λ) and up to about 4 × λ at MS3 even after the per-ion calibration: an MS3 confirmation p is a
confirmation, not a significance to threshold on.
## Offline coupling — FLASHDeconv reads `scan_commands.tsv` (ADR-0046)

`FLASHIda/ScanCommandJoin.h` is the one piece of FLASHIda code that runs **offline**, inside
FLASHDeconv. It is header-only, stateless, and knows nothing about FLASHDeconv — which is what lets
it be tested here (`FLASHIda_Logging_test`, the twenty `scan_commands_*` sections) when no CI job
runs FLASHDeconv at all.

```
mzML spectrum ── "scan description"[0..2] = tracking id ──► scan_commands.tsv row
row.decision_survey_id (ADR-0060) ── an ms_level 1 row ──► the survey whose ranking last chose the window
  (absent column / empty cell: row.parent_tracking_id ── walked UP the rows to the ms_level 1 row)
FLASHDeconv then takes ITS OWN PeakGroup nearest row.mono_mass from THAT survey
```

- **The tracking id survives conversion.** msconvert writes the scan-description trailer as the
  spectrum-level userParam `scan description`, which `MzMLHandler` loads as a meta value. ADR-0035
  and the glossary used to call the instrument scan number "the only one that survives"; it is not.
- **The row is a locator, not a measurement.** Commands queue, so 88 % of commanded MS2 have 1–5
  newer surveys between them and the survey their command came from. FLASHDeconv used to search the
  *latest* one and agreed with the commanded mass 40–46 % of the time, against 95 % when it happened
  to read the right one. Its deconvolution of the right survey holds the commanded mass for 99.8 %
  of MS2, so the reported mass, charge range and feature linkage stay FLASHDeconv's own.
- **A rebuilt precursor is the last resort, not the second** (decision 3, amended 2026-09-20). The
  commanded mass wins its window at the **same isotope only**; a window holding something else takes
  that species on **charge-SNR**, exactly as an uncommanded scan does; only a window holding nothing
  at all is rebuilt from the row by `peakGroupFromCommand_`. The rebuild has one peak and no isotope
  envelope, and `SpectralDeconvolution::performSpectrumDeconvolution` bounds the **fragment**
  deconvolution's max charge and max mass by the precursor PeakGroup — so preferring it over a
  measured species costs the MS2, not just the mass. `MAX_ISOTOPE_OFFSET` and `massRank`'s ladder
  stay in the header, now with no reader outside its own tests: the call site requires rank `== 0`.
- **A follow-up MS2's parent is the MS2 that triggered it**, not the survey
  (`ScanCommandQueue::buildFollowUp`), which is why `join` *walks* the parent chain rather than
  reading one hop. Since ADR-0060 the walk is the **fallback**: `surveyScanNumber_` first reads
  `decision_survey_id`, because a re-planned command's survey is the newer one, and a bought `R`'s parent
  is its `Q` while its window may have been re-chosen on a later survey.
- **The join fails closed, before any deconvolution.** Tracking ids restart in every run, so a
  foreign `scan_commands.tsv` joins 100 % of scans by id; the row-vs-spectrum check (MS level, anchor
  m/z within 0.01) is what tells them apart. A duplicate id, or a file that joins nothing, aborts too.
  A spectrum with no id or no row is an ordinary **uncommanded scan** and is never an error.
- **MS2 only.** MS3 rows are parsed and level-checked; nothing locates an MS3 yet — no acquired run
  exists to show how pwiz writes a two-stage precursor list.
- ⚠️ **`parse` resolves columns by header name and requires no column added after 2026-09.** Every
  run folder already on disk must keep parsing; `scan_commands_parse_reads_a_pre_adr0046_run_folder_file`
  pins that with three rows copied from a real one. A column *rename* in `IdaLogger` breaks
  `scan_commands_parse_roundtrips_an_engine_written_file` — which is the point.

The FLASHDeconv side (`FLASHDeconvAlgorithm::run` + `findPrecursorPeakGroupsForMSnSpectra_`, and
`FLASHDeconvTabWidget`) is **inside the no-go boundary** and was edited by an owner-approved, locked
line list. The boundary has not moved.

## Co-isolation notches (ADR-0016, ADR-0019)

**`num_stages` counts cascade stages and NEVER notches.** This is load-bearing: `processScan`'s
context gate, `syncEnergyMirrors_`, `Exploration`'s `si = num_stages - 1` and C# `ScanFactory`'s
clamp all key on it, and stay correct only while it means MSⁿ depth.

Notch descriptors live in their own **`Notch notches[18]`** (24 bytes each, geometry only), carved from
`reserved_` at @1464 alongside `pad4` (596 → 152 across two changes). **Stage k owns the fixed block
`[k * MAX_NOTCHES_PER_STAGE, + MAX_NOTCHES_PER_STAGE)`** — use `notchesForStage()` rather than
open-coding it, but write order no longer matters and the two stages cannot contend.

**Two different tens; do not collapse them.** `MAX_ISOLATION_STAGES` (10) is the `';'`-axis cascade
limit — every stage-carried iAPI key documents "a maximum of 10 values", including `ActivationType` and
`CollisionEnergy`, which have no notch axis at all. `MAX_NOTCHES_PER_STAGE` (9) + the anchor is the
`','`-axis limit, from `MSXTargets`' 10 MSX windows **per fragmentation stage**. So an MS3 can command
20 isolation windows. ADR-0017 read the first ten as a joint budget, which gave 8 slots shared across
an MS3's two stages; stage 0's inherited precursor set consumed all of them and the fragment stage got
none. `ScanCommandLayout_test` pins the two constants apart deliberately.

A `Notch` has no collision energy or activation field, which is the correct encoding of "all notches of
a stage fire into one fragmentation event" — under the old layout the writer copied `stages[k]`, so a
per-notch CE existed and merely happened to agree.

`NotchSelection.h` holds the whole policy — SNR **gate**, the authored **`charge_selection`** order
(`chargeRanksBefore`: the rule's key, then intensity, then the lower charge; ADR-0059), clamp-with-report — plus `peakGroupNotchCandidates`, shared by `PrecursorSelection` and
`buildMS2` (which writes the geometry), so the set recorded as acquired is by construction the set
isolated. Until ADR-0059 SNR admitted and **intensity** ranked, on the argument that SNR-ordering lets a
clean-but-faint charge displace an abundant one under a clamp; the ranking is now authored
(`snr | qscore | intensity | purity`, default `snr`, `purity` refused under `multiplexed`) and the same
rule picks the anchor.

⚠️ **Under `snr` and `qscore` the anchor must be a RESOLVED charge** (`fix(FLASHIda)` `463d2ba1a4`, made
visible by the FDmerge merge). `pickAnchor` (`CandidateAdmission.h`) keeps
`getRepAbsCharge()`/`getBestQScoreCharge()` only if that charge is in `resolvedCandidates_`, the same
list every other rule, the notches and the `separate` set read. A single-peak charge has no measured
window (`getMzRange` lo == hi) and is absent from that list. Since FDmerge's noise changes, such a charge
can top the per-charge SNR (measured: 6306.15 Da rep z6 at SNR 1.55 against ≤ 0.86 for the full
envelope). When the accessor's charge is unresolved, the first of `rankCharges(cands, rule)` is the anchor.
With no candidate at all, the result is `NoAnchorResolved`. When the accessor's charge is resolved, the
result is byte-identical.

**Three of the policy's inputs are now AUTHORED, and all three are inert by default** (ADR-0040,
`precursor_selection`). `snr_threshold` is the same gate as ever but reachable — it used to be
assigned a hardcoded `1.0` at the very END of `parseJSONConfig_`, after every section had been read,
so adding the key without deleting that assignment would have left it bound and inert.
`max_charge_states` bounds the SET and reaches `selectNotches` as `max_charge_states - 1`, since the
anchor is not a notch; `MAX_NOTCHES_PER_STAGE` remains the hard ceiling, so an authored value can
only lower it and the 2048-byte ABI is never at risk. `min_charge_states` is a floor applied in
`admitCandidate` — it refuses the species outright rather than trimming its set, so no
`max_precursors` slot is spent, and it **exempts a sibling**, whose set is deliberately reduced by
what the species already isolated (ADR-0036); without that exemption the floor refuses precisely the
PeakGroup supplying the last charge of a split envelope.

**The set cap has one definition**: `maxAcquisitionSetSize(max_charge_states)` in `NotchSelection.h`
(`max_charge_states`, at most `MAX_NOTCHES_PER_STAGE + 1`; `0` means that limit). Admission and the
`buildMS2` fallback pass it minus one as `max_notches`; the `separate` remainder of a re-plan
(`FLASHIda.cpp`) reads it as its cap. Before 2026-11-12 the re-plan read the raw key, so `0` was
unbounded there and `12` meant 12.

**What `min_charge_states` counts depends on the mode** (ADR-0040, amended 2026-11-12). Under `separate`
and `multiplexed` it is the size of the acquired set. Under `single` the set is always the anchor, so the
gate counts the charges the species **resolves**: `resolvedChargeCount` (`CandidateAdmission.h`) is the
anchor plus `admissibleNotches` of the species' candidates, restricted to the authored set when the row
names one. A charge acquired on an earlier survey still counts. `max_charge_states` does not cap this
count. It is computed only when `min_charge_states > 1`.

**The charge floor bounds every charge of a set, not only its anchor.** `selectNotches` takes it as a
required argument (`min_charge`, after `snr_threshold`), so a caller cannot leave it out:
`acquisitionChargeSet` and the `buildMS2` fallback pass level 1's (`min_precursor_charge`), the two
`ProteoformTracker` calls pass level 2's (`min_fragment_charge`), and `rechooseCharges` applies level 1's
to the `separate` remainder and the bought `R`. `admissibleNotches` is the gate alone (not the anchor,
isolatable, `snr >= snr_threshold`, at or above the floor): it does not clamp and does not print
`[NOTCH-CLAMP]`, which is why the `single` count reads it. Before 2026-11-12 the floor reached only the
anchor of an unrestricted species, and the notches of an ungoverned MS3 plan had none.

⚠️ **The clamp mattered most where it meant least.** `buildMS2` writes notches only under
`Multiplexed` (`ScanCommandQueue.cpp:316`), yet `separate` — which never writes one — was clamped by
the notch array's capacity all the same, because it calls `selectNotches` upstream purely to compute
its charge LIST. Measured on the `separate_charges` golden: 10 of 18 species-surveys sat at exactly
10 charges, and 100 of that mode's 135 MS2 scans came from species pinned at that coincidence.

**Both MS2 non-single modes read that one call, and this is load-bearing** (ADR-0021). `multiplexed`
turns the set into notches on one command; `separate` emits one command per member. They differ only
in scan count. `separate` used to source its charges from `charges_to_process` instead — a list that
was multi-valued only under the `charge_based_exclusion` developer flag — so acquisition GEOMETRY
depended on an exclusion-KEYING flag and the mode was inert wherever that flag sat at its default,
i.e. everywhere. The flag is deleted; exclusion is mass-keyed.

**A multiplexed scan shares its `max_it` across its stage-0 windows** (ADR-0072). The instrument applies
`MaxIT` to each fill, so a 10-window command sent at the configured value could inject for ten times it.
`sharedMaxIt_(scan_budget, cmd)` in `ScanCommandQueue.cpp` returns `scan_budget / windowsAtStage(cmd, 0)`,
and `buildMS2`, `buildMS3` and `buildFollowUp` each assign it from their own scan config once the command's
stage-0 notches are written. An MS3 divides by the precursor windows it replays and never by its own
stage-1 notches, which are one simultaneous fill. `reisolate` has no scan config: it recovers the budget
as `r.max_it * windowsAtStage(r, 0)` and shares it over the new windows. `makeMS1` and `makeAGC` are
untouched, an unset `max_it` (0) stays 0, and there is no floor. Nothing in the engine reads `cmd.max_it`
back; `scan_commands.tsv`'s `max_it` column is the value sent.

**An inclusion row may name WHICH charges** (ADR-0028), which `precursor_charges` cannot — it is
all-or-one. `PrecursorSelection::authoredChargesFor_` unions the charge sets of the RT-active rows
whose mass matches, and that **authored charge set** filters `peakGroupNotchCandidates` before
`selectNotches` runs, so the SNR gate and the intensity ranking above are unchanged — the set only
shrinks what they see. It can never extend: a charge the deconvolution did not resolve has no
measured window and is skipped, and `min_charge` still applies.

Two things move with it. The **anchor** becomes the highest-SNR *named* charge rather than
`getRepAbsCharge()`/`getBestQScoreCharge()`, resolved before the selection loop because it needs the
set; and its **own** per-charge qscore (`PeakGroup::getAllQscores()`) is what gets logged and
excluded on, rather than a sibling's. Exclusion is re-keyed to `(nominal mass, charge)` in
`authored_acquired_rt_map_` **for these species only** — everything else stays mass-keyed as
ADR-0021 left it — so `single` walks the set across surveys, one charge per survey, and the species
is skipped once every named charge is spent. A call-local guard stops non-strict inclusion's phase-1
pass taking a second named charge within one survey.

This also replaced a matcher that walked every RT-active row and took the first whose charge fell
inside the envelope **without checking that row's mass**. It was invisible because all committed
inclusion lists write `-1`, which returns on the first row before the mass ever mattered.

⚠️ **The fan-out lives at the emit loop, below the mass-level bookkeeping — never in the candidate
loop.** That bookkeeping runs once per species and both its guards
(`tqscore_exceeding_mass_rt_map_`, and the "previously acquired with higher qscore" skip) key on
`nominal_mass`, so iterating charges through it `continue`s on every sibling. Anyone "fixing"
`separate` by making the candidate loop multi-charge will reproduce the original bug.

**The MS3 side needs the fragment's whole envelope, and getting it there is the part that broke.**
`ProteoformTracker::feedScan` flattens each fragment PeakGroup to `getMaxIntensityAbsCharge()` for the
representative fields, and *also* stores every present charge as a `ChargeRecord` in
`PeakRecord::by_charge`; `upsertMappedObservation_` fills `ms2_by_charge` from that. Without the
envelope, `ms2_by_charge` held one entry per fragment per scan and `identification.characterization.fragment_charges`
was inert in **both** on-values while the spectra resolved ~38% of fragments at two or more charges.
`planNextScans` derives `separate` and `multiplexed` from **one** `selectNotches` call, so the two modes
acquire the same set and differ only in scan count. That call takes `min_fragment_charge` as its floor,
so no charge below it is a notch or a `separate` scan of an ungoverned plan (a governed or yield-gated
plan floors its candidates earlier, in `gateAdmitsCharge`).

Two things notches deliberately do **not** carry: per-notch **scores** (the slots hold geometry only,
and the `*_s1` block is two-stage, so `charge_snr`/`precursor_intensity` remain the anchor's), and a
redefined **`window_snr`** (it is carried onto the follow-up command and into the MS2/MS3 proteoform
contexts, so it stays the anchor window's purity rather than a union; since ADR-0048 no metric scores
from it — the precursor share is read from the raw scan instead).
- **Delimiters: `';'` everywhere EXCEPT id lists, which use a space.** `child_ids` and
  `contributing_scan_ids` join with `' '` because base-94 tracking ids draw from 0x21–0x7E, which
  **includes `';'`** (an id like `!!;` would collide); space is the only printable char the alphabet
  excludes. Readers on both sides must split those two columns on `' '`.
- **Column order was permuted for legibility and goldens were NOT recaptured.** The C# comparison
  resolves columns **by header name**, so a further pure reorder is free — but a rename/add/drop
  fails closed. The C++ `FLASHIda_LoggingFields_test` hard-codes the **new** indices, so any reorder
  *does* require editing that test.

## Tests

Registering a test in `src/tests/class_tests/openms/executables.cmake` is **not enough to run it**.
A C++ test executes **at all** only if it appears in **both** places in `../.github/workflows/flashida-ci.yml`: the build
`--target` list **and** the `ctest -R` alternation. Miss the first and it never builds; miss the
second and it builds but never runs. Neither list is reproduced in prose here or anywhere else —
any copy you find is stale by construction; `../CLAUDE.md` owns how the two are read out of the yml.

Tests read fixtures from `../../FlashIDA/test-data` relative to the **build directory**, so the
FlashIDA submodule must be checked out and the build dir must sit directly under `OpenMS/` with
`FlashIDA/` beside it — that is what makes the `"../../FlashIDA/test-data/…"` literals resolve.
Test exes land in `<build>/src/tests/class_tests/bin/`, not `<build>/bin/`. **Staging the 5-DLL set
beside them (including `zlib.dll`) is Windows-only** — on Linux the equivalent is
rpath/`LD_LIBRARY_PATH` to `libOpenMS.so`, and there is no `zlib.dll` at all.

**Drive acquisitions only through `FLASHIda_TestHelpers.h::runInterleaved`** — the C++ mirror of
C# `PushScanAndDrainFull`. One contract: pull a command → classify idle vs workload → feed one
response scan stamped with that command's own engine-emitted description → repeat; terminate on
`idle >= 3`. You cannot fabricate a scan id (gate 3 above), which is precisely why hand-rolled
drive loops don't work. `FLASHIda_ProcessScan_test` pins this with
`processScan_ms1_gate_rejects_unrequested_id`.

**Every drive waits for the search pool after every feed** (ADR-0062 d.7): `runInterleaved` and the
hand-driven helpers (`pushAllScans`, `bootstrapExplorationGroup`, `driveStaleSurvey`) call
`ida->waitForSearchIdle()` after each `processScan`, so the next pull sees exactly what a synchronous
`processScan` left queued; `runInterleaved` adds its return to `ms2_feed_returns`. The sections that pin the
pool itself (`FLASHIda_ProcessScan_test`, `FLASHIda_exploration_test`) deliberately skip the wait and use
`waitUntil` and `SearchGate` (`FLASHIda_TestHelpers.h`) with the `FLASHIda_TestAccess.h` seams
`onBeforeSearch`, `hasSearchPool`, `searchThreads`, `searchesFinished`, `appliesFinished` and
`setSearchBacklog`.

Division of labour with the C# suite: **C++ ctests assert plausibility ranges** (stable across
engine bumps); **C# NUnit asserts exact goldens**. Put a new numeric expectation on the C# side
unless it is genuinely range-based.

**One toolchain builds this code: MSVC, Release** — `OPENMS_PRECONDITION` and debug asserts are
compiled out, an accepted tradeoff to match the production toolchain. No gcc build exists (ADR-0050).
Do not count on assertions here anyway: `grep OPENMS_PRECONDITION` over the whole TOPDOWN tree returns
**zero**, and the single `assert(` there is commented out.
