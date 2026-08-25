# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

**Scope: the FLASH real-time IDA engine only.** This is a C++20 fork of OpenMS, checked out as a
git submodule of the `flashida-development` workspace. It exists to serve FLASHIda; this file
documents the real-time acquisition engine and nothing else. Upstream OpenMS subsystems
(FORMAT, CHEMISTRY, the TOPP tools, pyOpenMS) are out of scope — treat them as a vendored library.

The parent `../CLAUDE.md` owns the bridge/ABI contract, CI, build and test commands, and the
config flow, and wins on any conflict. `../FlashIDA/CLAUDE.md` owns the C# side.

**Do not build this project unless explicitly asked** — it is resource-intensive and CI handles it.

## Where the code lives

Headers under `src/openms/include/OpenMS/ANALYSIS/TOPDOWN/`, sources under
`src/openms/source/ANALYSIS/TOPDOWN/` — note the source tree has **no** `OpenMS/` path segment
(`source/OpenMS/…` does not exist). `sources.cmake` is **not** a reliable header inventory:
`Ms2Params.h` and `ProteoformTracker.h` are live but unregistered. Glob the directory instead.

### The untouchable boundary

The directory split is a near-perfect proxy for what you may modify.

| | Files | Status |
|---|---|---|
| **FLASHDeconv** (flat `TOPDOWN/`) | `DeconvolvedSpectrum`, `SpectralDeconvolution`, `FLASHDeconvAlgorithm`, `FLASHHelperClasses`, `MassFeatureTrace`, `PeakGroup`, `PeakGroupScoring`, `Qvalue`, `TopDownIsobaricQuantification` | **OFF-LIMITS** |
| **FLASHTnT** (flat `TOPDOWN/`) | `FLASHTnTAlgorithm`, `FLASHTaggerAlgorithm`, `FLASHGappedTaggerAlgorithm`, `FLASHExtenderAlgorithm` | **OFF-LIMITS** |
| **FLASHIda engine** | `FLASHIda.{h,cpp}`, `FLASHIdaBridgeFunctions.{h,cpp}` (flat) + everything in `TOPDOWN/FLASHIda/` | Fair game |

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
  Tuning knobs arrive via the `flashtnt` config section, not through `FLASHTnTAlgorithm`'s own
  defaults.
  **Every Param handed to either algorithm is built in one of two places** —
  `FragmentAnalysis::buildTaggerParam` / `buildExtenderParam` (:197 / :213, which construct a
  throwaway instance purely for `getDefaults()`). All three call sites route through them, so a new
  knob is added once, not per site. `fixed_mod` is set **unconditionally** there: an empty config
  list is passed verbatim as an empty Param list, *not* left at the algorithms' declared `{""}`
  placeholder — pinned by `FragmentAnalysis_test`'s `build*Param_sets_fixed_mod_when_empty`.
  ⚠️ Two keys that look like FLASHTnT tuning are **not**: `max_ptm_count` and
  `max_flanking_mass_diff` are neither Params nor read by either algorithm on the IDA path. They
  come from `precursor_selection.tag_expansion` and drive FLASHIda's own FASTA target expansion.

> **Name trap.** `TopDownIsobaricQuantification` (flat, off-limits, used only by
> `FLASHDeconvAlgorithm`) and `FLASHIda/Quantification` (IDA, fair game) are unrelated. The IDA one
> hardcodes `TMTSixPlexQuantitationMethod` + `IsobaricChannelExtractor` with a
> `// TODO: Variable channel extractors` — the real-time quant test is **6-plex only**, whatever
> label was actually run.

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

### Enqueue discipline: return, don't push

`FLASHIda.cpp` is the **sole enqueuer** — exactly nine `queue_.push(` sites; no `.push(` on a
`ScanCommandQueue` exists anywhere else under `TOPDOWN/`. The `ScanCommandQueue&` that
`Exploration::initiate/feedResult/initiateNextLevel` receive is a **command builder, not an
enqueue target**: `buildMS2`/`buildMS3`/`buildFollowUp` allocate a tracking id and return a
`ScanCommand` by value. Any new command-producing component must follow the same contract.

They also do **not** register in the pending-scan map — that happens only in `dequeue()` and in
explicit `registerPending()` calls, which `FLASHIda.cpp` makes itself on the drain path.
⚠️ The doc comment at `ScanCommandQueue.h:59-60` claiming builders register is **stale**.

## `processScan` — five admission gates

`FLASHIda.cpp:79-169`. Before any MS-level branching, a scan must clear all of these; each returns
`0`. This is the always-on **engine-id-echo** contract: a spectrum is analysed only if the engine
itself minted and dispatched its id.

| # | Gate | Trace |
|---|---|---|
| 1 | `scan_description` (or null) shorter than 3 chars | silent |
| 2 | `desc[3] == 'A'` — AGC calibration scan; resolves pending, returns | silent |
| 3 | `peekPending(decode(desc[0..2]))` empty — id never emitted | `[TRACK-RESOLVE] status=not_found` |
| 4 | `resolvePending` empty after a successful peek (upstream race) | `status=context_lost_race` |
| 5 | context support: `required_stages` = 0/1/2 for ms_level 1/2/3; rejects if ms_level ∉ {1,2,3}, or `parent_ctx.msn_level != ms_level`, or `parent_ctx.num_stages < required_stages` (a **less-than**, not an equality) | `status=context_unsupported` |

A sixth early return follows the gates: MS1 with `config_.level(1).selection == None` returns 0
before selection. None of these paths reaches the TSV writers, so an un-commanded scan produces no
log row — but gates 3–5 do leave a stdout trace.

### Scan-description wire format

`{3-char base-94 tracking id}{1-char type marker}{payload}`. Markers, per the authoritative switch
in `IdaLogger::scanTypeFromDescription_`:

`S` survey MS1 · `A` AGC calibration · `R` "recording" — any data-acquiring MS2 *or* MS3 ·
`F` quantification follow-up MS2 · `C` tagging conditional follow-up MS2 · `E` exploration variant.

`processScan` branches on the marker twice: `desc[3]=='A'` (gate 2) and
`is_follow_up_scan = desc[3]=='F' || desc[3]=='C'`.

> The ABI field is `char scan_description[256]`, but **every writer caps at `snprintf(dst, 16, …)`
> — 15 chars + NUL.** That is why MS2/MS3 mass tokens go through
> `ScanCommandQueue::formatMassToken`, which computes a token budget and degrades decimal precision
> from 6 down to 0 so the trailing `@{charge}{ion_type}{index}` is never truncated.
> Known dead store: `FLASHIda.cpp:661` sets `[3]='C'` but :665 overwrites the buffer, so
> cycle-time-forced MS1s emit as `'S'`.

## `getNextScanCommand` — never returns 0

Five steps: (1) scheduled AGC prescan → (2) cycle-time MS1 → (3) cleanup → (4) priority dequeue →
(5) idle fallback. **Every path returns 1.** When all four queues are empty, Step 5 mints a fresh
priority-3 survey MS1, pushes it, and re-enters `dequeue()` to return it — the instrument is never
observably starved. The only `return 0` in the whole file belongs to `processScan`.
`FLASHIda.h:99`'s "0 if error" comment is imprecise.

Callers must therefore terminate on `msn_level == 1 && priority == 3`, or bound their iterations.
See `../CLAUDE.md`.

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
  critical sections. `analysis_mutex_` is held only around the three `writeScanCommandRow` calls.
  The only lock-free `processScan` → `getNextScanCommand` channel is two atomics
  (`exploration_active_`, `current_faims_cv_`), written release / read acquire.

This matters because the two bridge calls genuinely run on **different threads**: `processScan` on
the C# TPL Dataflow `ActionBlock` thread (serialized with itself, `MaxDegreeOfParallelism = 1`),
`getNextScanCommand` on the instrument's event thread.

### Priority ladder (0 highest → 3 lowest)

| P | Commands |
|---|---|
| 0 | FAIMS CV-transition MS1; cycle-time MS1; **both follow-up kinds** (quant `F`, conditional `C`) |
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

## Where deconvolution actually runs — 3 call sites, 2 engines

- **MS1** is *not* deconvolved from `processScan`: it happens inside
  `PrecursorSelection::filterAndRank` → `deconv_.deconvolveMS1(…)`, just before mass-exclusion.
- **Regular MS2/MS3** — `deconv_.deconvolveMSn(…)`, called directly from `processScan`.
- **Exploration variants use a different engine.** `Exploration` owns a *second* `Deconvolution`
  instance built with `config.explorationToleranceList()` rather than `config.toleranceList()`.
  So `deconv_.storedMS2()` is empty/stale during an exploration group, and exploration results can
  use a different ppm tolerance than production scans.

`target_mode` semantics (and the fact that 2/3 are swapped in every doc comment) are documented in
`../CLAUDE.md` — mode **2 is deep/in-depth**, mode **3 is exclusion**.

## Characterization / MS3

**`ProteoformTracker::planNextScans` is the only producer of MS3 targets.** `initiateNextLevel`
builds commands exclusively inside `if (model != nullptr && !model->proteoform_sequence.empty())`;
every other path returns zero commands. No tracker, no finalized model, or an unidentified
precursor ⇒ **zero MS3 by design**, with no intensity or legacy fallback.

- **`characterization.mode`** (`off | ambiguity | coverage | exhaustive`) is the single MS3 switch
  and also picks the targets: each of its **three** on-values *is* an objective, one for one
  (ADR-0013, extended by ADR-0023). `Off` means no MS3 at all.
  **Unknown values throw** — the old `objective` parse was `if (== "coverage") … else Ambiguity`, so
  `"Coverage"` or any typo silently meant ambiguity, and with `mode` now carrying the on/off bit a
  typo'd `"Off"` would have silently *enabled* MS3. The struct still exposes
  `characterization().objective`, derived from `mode`, so every existing read site is unchanged.
  ⚠ **`mode` itself has zero read sites outside `Config.cpp`** — it is parsed, projected and quoted
  in error text there and nowhere else; the engine branches on `.objective`
  (`ProteoformTracker.cpp:421`, `:524`). A new mode that does not also assign an objective therefore
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
  Its pool floor is the new `characterization.min_target_mass` (Da, `0` = off), and an unassigned
  target's CE sweep is forced onto `RemainingPrecursor` whatever the config asked for — gated
  `msn_level >= 3`, because `initiate`'s `ion_type` defaults to `'\0'` and an ungated test would
  drag every MS2 exploration group with it (ADR-0023 decisions 5 and 11, D-f, D-g).
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
  `characterization.max_targets` / `.min_fragment_charge` and projected onto level 2 by
  `Config::applyCharacterizationMode_`, so the read sites are unchanged while the keys finally sit
  where the feature does. `selection_strategy.ms3.max_targets` / `.min_charge` were parsed and
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

### One un-fragmented reference per swept activation (ADR-0029)

`exploration.activations` is a **list**, and `buildVariants_` emits each activation's whole sweep as a
contiguous block before starting the next. A **baseline** heads each block: that activation's own
variant with the **swept axis alone** turned off — CE for HCD/CID, reaction time for ETD, both for
EThcD — and every other parameter left as its siblings carry it, so an ETD baseline keeps
`base_config.collision_energy` rather than dropping to 0.

⚠️ **The two axes turn off at different values, and the asymmetry is the instrument's.** A collision
energy of 0 is commandable and simply does not fragment; a reaction time of 0 is **rejected**, so "no
reaction" is `Config.h`'s `MIN_REACTION_TIME_MS` (0.03 ms, ~300× shorter than a working ETD time).
Zeroing both would emit an ETD reference the instrument refuses to acquire. The authored sweep *grid*
is deliberately **not** floored — `Config::validate` rejects a `reaction_time_min` below the floor
when an ETD-family activation is swept, which is what lets the grid's first point still coincide with
the baseline and be suppressed. Gate it on `needsReactionTime`, never unconditionally: every
CE-only config leaves `reaction_time_min` at its 0 default.

Three consequences that are easy to get backwards:

- **It is suppressed, not duplicated, when the sweep already contains its turn-off point.** With
  `ce_min: 0` (or `reaction_time_min: 0.03`) the block's first sweep point *is* the reference, so it
  is flagged in place and nothing is synthesized. `buildVariants_` decides this by comparing the two **emitted commands**, not by
  testing `ce_min == 0` — which is what makes the flagged and synthesized forms provably the same scan
  and covers EThcD's cross-product for free. `initiate` prepends **nothing**; it used to insert one
  baseline stamped `base_config.activation`, which need not be a swept activation at all.
- **An activation that sweeps neither axis gets no baseline** and competes normally. A baseline is
  "this scan with fragmentation off"; with no swept axis there is nothing to turn off. Drop that guard
  and such an activation's single variant is identical to its own baseline, gets flagged, and can never
  win.
- **`ExplorationGroup::baseline_intensity` is a `std::map<std::string, double>` carrying three
  states**: absent = that activation's reference has not returned; `> 0` = usable denominator;
  `<= 0` = returned empty. `has_baseline` and `baseline_failed` are gone — both were group-scalar facts
  that are now per-activation.

⚠️ **An empty reference bars one activation; it cancels nothing.** Its variants are still acquired and
score **`-1.0`**, and that value is load-bearing rather than stylistic: winner selection seeds
`best_score = -1.0` and takes `score > best_score`, so `-1.0` is excluded by the existing comparison
with no new flag — while `0.0` **wins the group at zero**, and at MS3 a measuring metric then fires a
production scan at that coin-flip CE. Siblings are unaffected and supply the winner; only when every
activation is de-referenced does `best_idx` stay `-1`. `ScanCommandQueue::cancelByScanIds` therefore
has no production caller left and is kept as a tested queue primitive.

### Reading vs measuring exploration metrics — the MS3 evidence rule (ADR-0020)

**Only `FragmentCount` identifies its own MS3 pre-scans.** It scores a variant *by* matching it, so
`ProteoformTracker::scoreCalibratedVariants` — the only matcher that works on an MS3 spectrum, since
it is seeded with the proteoform context and the parent ion index — runs inside a
`metric == FragmentCount` guard in `Exploration.cpp`. `MassCount` and `RemainingPrecursor` are
**measuring**: they score from peak counts and window intensity and never invoke a matcher.
`isMeasuringMetric()` (`Config.h`, beside the enum) is the classification; do not open-code
`== FragmentCount`.

Three things follow, and the third is the trap:

- `ExplorationVariant::identification_result` has **one write site**, inside that guard. Under a
  measuring metric it is empty for every variant.
- The inline fold at MS3 group completion is gated on that field, so under a measuring metric it
  cannot fire — which is why a completed MS3 sweep used to contribute *nothing*: no
  `identification.tsv` row, no `foldMs3`, no pooled evidence, after paying for every pre-scan.
- Therefore a **measuring metric at MS3 always emits a production scan**, regardless of `overrides`.
  The post-winner dispatch fires on `!overrides.empty() || (msn_level >= 3 && isMeasuringMetric(...))`.
  The follow-up returns on the *regular* MS3 path and is matched there.

`overrides` means *"the pre-scans ran at degraded settings"* — variants are built from
`scans[0] + overrides`, the production scan from `scans[0]` alone. That is the whole reason the
original gate keyed on it. **MS2 is exempt** from the metric condition: `computeFragmentMatch_`'s
whole-protein matcher is correct at MS2, so every MS2 variant is identified under every metric and an
MS2 group cascades instead of re-acquiring.

`planNextScans` emits `[MS3-PLAN] … reason=…` naming which of its zero-target causes fired
(`all_mods_localized` is the common, *correct* one: ambiguity mode with everything already localized).
stdout only — like every engine marker, it is discarded during instrument acquisition.

## Config (`FLASHIda/Config.cpp`)

Unknown keys are hard-rejected by one free function `rejectUnknownKeys(obj, allowed, path)` called
~18 times, each with its own hand-written allowed-set, throwing `std::invalid_argument` naming the
offenders and pointing at `FlashIDA/test-data/config_schema_reference.json`. Two exemptions:
`exploration.overrides` (a dynamic string→string map) and a top-level `ms3`, which gets a dedicated
**migration** error before the root check runs.

> **Resolved:** `selection_strategy` is deleted, and with it the one `std::runtime_error` in the
> file. A *present* `selection_strategy` now throws `std::invalid_argument` with a migration message
> naming all seven destinations, so every config error in this file is now the same exception type.

Traps worth internalizing before touching config code:

- **`.value(key, default)` fallbacks are dead in production.** C# `ToCppJson` emits every key
  unconditionally, so the C++ literals never fire and several disagree with the effective C#
  defaults. Don't read a C++ default and believe it.
- **`kScanKeys` is one lenient 17-key union validating every scan object, and one parser reads all
  17 at all five sites.** Validation and parsing share a single source of truth (`Config.cpp:70-73`
  records why: when the two lists drifted, nine keys passed validation and were then discarded —
  notably `reaction_time` on a `follow_up_scan`, which made an ETD follow-up unconfigurable). The
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
  `ScanCommand`). Only the two `follow_up_scan` blocks keep `"Orbitrap"`.
- **Selection is never authored per level any more** — it is projected from `characterization.mode`
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
| `ida.log` | — | free-text MS1 summary (not a TSV) |
| `scan_commands.tsv` | 34 | one row per **dequeued** command; the wide MS3-fragment stream (32→34 ADR-0026 `first_mass`/`last_mass`, between `faims_enabled` and the trailing `enqueue_ts`) |
| `scan_results.tsv` | 32 | pure acquisition-**event** log per `processScan` (34→29 identification payload moved out, →28 per-charge deconv restructure, →29 `deconv_qscores`; →32 the identification-YIELD block `tag_count`/`fragment_count`/`tic_coverage` after `remaining_ratio`) |
| `identification.tsv` | 34 | per-scan MS2/MS3 identification leaf (32→34: `tag_count` beside `flash_extender_score`, `fragment_qscores` inside the aligned fragment-mass table) |
| `pooled_identification.tsv` | 19 | per-precursor cumulative proteoform trajectory |

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

`NotchSelection.h` holds the whole policy — SNR **gate**, descending-**intensity** order (charge as the
tiebreak), clamp-with-report — plus `peakGroupNotchCandidates`, shared by `PrecursorSelection` and
`buildMS2` (which writes the geometry), so the set recorded as acquired is by construction the set
isolated. SNR admits and intensity ranks: SNR is a purity measure, so ordering by it would let a
clean-but-faint charge displace an abundant one under a clamp and trade away the ion current the fill
exists to harvest.

**Both MS2 non-single modes read that one call, and this is load-bearing** (ADR-0021). `multiplexed`
turns the set into notches on one command; `separate` emits one command per member. They differ only
in scan count. `separate` used to source its charges from `charges_to_process` instead — a list that
was multi-valued only under the `charge_based_exclusion` developer flag — so acquisition GEOMETRY
depended on an exclusion-KEYING flag and the mode was inert wherever that flag sat at its default,
i.e. everywhere. The flag is deleted; exclusion is mass-keyed.

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
envelope, `ms2_by_charge` held one entry per fragment per scan and `characterization.fragment_charges`
was inert in **both** on-values while the spectra resolved ~38% of fragments at two or more charges.
`planNextScans` derives `separate` and `multiplexed` from **one** `selectNotches` call, so the two modes
acquire the same set and differ only in scan count.

Two things notches deliberately do **not** carry: per-notch **scores** (the slots hold geometry only,
and the `*_s1` block is two-stage, so `charge_snr`/`precursor_intensity` remain the anchor's), and a
redefined **`window_snr`** (it feeds exploration's `RemainingPrecursor` scoring, so it stays the anchor
window's purity rather than a union).
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
A C++ test executes in CI only if it appears in **both** places in
`../.github/workflows/flashida-ci.yml`: the build `--target` list **and** the `ctest -R`
alternation. Miss the first and it never builds; miss the second and it builds but never runs.

Tests read fixtures from `../../FlashIDA/test-data` relative to `OpenMS/build`, so the FlashIDA
submodule must be checked out. Test exes land in `build/src/tests/class_tests/bin/`, not
`build/bin/`, and need the 5-DLL set staged beside them (including `zlib.dll`).

**Drive acquisitions only through `FLASHIda_TestHelpers.h::runInterleaved`** — the C++ mirror of
C# `PushScanAndDrainFull`. One contract: pull a command → classify idle vs workload → feed one
response scan stamped with that command's own engine-emitted description → repeat; terminate on
`idle >= 3`. You cannot fabricate a scan id (gate 3 above), which is precisely why hand-rolled
drive loops don't work. `FLASHIda_ProcessScan_test` pins this with
`processScan_ms1_gate_rejects_unrequested_id`.

Division of labour with the C# suite: **C++ ctests assert plausibility ranges** (stable across
engine bumps); **C# NUnit asserts exact goldens**. Put a new numeric expectation on the C# side
unless it is genuinely range-based. Note CI builds **Release**, so `OPENMS_PRECONDITION` and debug
asserts are compiled out — an accepted tradeoff to match the production toolchain.
