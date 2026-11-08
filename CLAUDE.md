# CLAUDE.md

Guidance for Claude Code in the OpenMS submodule. **Scope: the FLASH real-time IDA engine only.**
This is a C++20 fork of OpenMS serving FLASHIda; upstream subsystems (FORMAT, CHEMISTRY, TOPP tools,
pyOpenMS) are a vendored library. The parent `../CLAUDE.md` owns the working rules, the bridge/ABI
contract, CI and landing, build and test commands and the config flow, and wins on any conflict;
`../FlashIDA/CLAUDE.md` owns the C# side. Each rule below is one line; the mechanism, measurements
and history are in `doc/FLASHIda-engine-reference.md`, under the same heading.

**Nothing builds this locally.** Push the branch: the `build` job compiles it on MSVC, Release, on
the self-hosted runners (ADR-0050). This repo is a **public** fork.

## Where the code lives

Headers in `src/openms/include/OpenMS/ANALYSIS/TOPDOWN/`, sources in
`src/openms/source/ANALYSIS/TOPDOWN/` (**no** `OpenMS/` segment). `FLASHIda.cpp` and
`FLASHIdaBridgeFunctions.cpp` sit in `TOPDOWN/`; the runtime helpers in `TOPDOWN/FLASHIda/`:
`Config`, `Deconvolution`, `Exploration`, `FAIMS`, `FragmentAnalysis`, `IdaLogger`, `KineticModel`,
`MS3FragmentMatcher`, `Ms2Params`, `PrecursorSelection`, `ProteoformTracker`, `ProteomeSearch`, `Quantification`,
`Refine`, `ScanCommand`, `ScanCommandQueue`, `ScanRole`, `SearchPool`, plus the offline `ScanCommandJoin`.
`sources.cmake` is not a reliable header inventory; glob the directory.

### The untouchable boundary

| | Files | Status |
|---|---|---|
| FLASHDeconv (flat `TOPDOWN/`) | `DeconvolvedSpectrum`, `SpectralDeconvolution`, `FLASHDeconvAlgorithm`, `FLASHHelperClasses`, `MassFeatureTrace`, `PeakGroup`, `PeakGroupScoring`, `Qvalue`, `TopDownIsobaricQuantification` | **OFF-LIMITS** |
| FLASHTnT (flat `TOPDOWN/`) | `FLASHTnTAlgorithm`, `FLASHTaggerAlgorithm`, `FLASHGappedTaggerAlgorithm`, `FLASHExtenderAlgorithm` | **OFF-LIMITS** |
| FLASHIda engine | `FLASHIda.{h,cpp}`, `FLASHIdaBridgeFunctions.{h,cpp}`, all of `TOPDOWN/FLASHIda/` | Fair game to propose |

Every change under `src/` (engine, bridge, tests, `.cmake`) needs the owner's approval before it is
written; an approved plan that names the files and describes the change **is** that approval. A
departure from the plan goes back to the owner. Ambiguous cases count as source changes.

- **Three owner-approved edits sit inside the boundary**: `PeakGroup::updateSNR_` (`9721595ca4`, a
  charge with no signal has SNR 0), `FLASHTnTAlgorithm`'s **prepared database** (ADR-0061 d.3):
  `prepareDatabase()` holds `run()`'s database side moved as it was, and `setPreparedDatabase()` makes
  `run()` reuse it read-only, and `FLASHTnTAlgorithm.cpp:1131` `if (! prepared_) setLogType(CMD);` (d.3 as
  amended for push 4c: a prepared `run()` keeps its caller's log type). Everything else in FLASHDeconv/FLASHTnT is exactly as Kyowon's FDmerge has it.
- **A prepared `run()` must never write the shared state**: it binds the mass ladders `const` and skips
  `preparePrecursorWindows`/`clearPrecursorWindows` (an unlisted precursor mass is swept: the same
  windows). `run()` keeps its results until `releaseResults()`, and FLASHTnTAlgorithm's copy
  constructor copies only its `Param`: construct each instance and call `setParameters`.
- `OptimizationMetadata` is IDA code in FLASHDeconv territory (a sanctioned hook with its own test).
- `FragmentAnalysis.cpp` and `PrecursorSelection.cpp` reach into FLASHTnT, driving Tagger and Extender
  directly (every Param built in `buildTaggerParam`/`buildExtenderParam`); `ProteomeSearch.cpp` drives
  `FLASHTnTAlgorithm::run` itself, one MS2 at a time on the prepared database. Nothing calls it yet
  (ADR-0061 push 4b wires it into the search pool); `FLASHIda_ProteomeSearch_test` pins it to FLASHTnT.
- `TopDownIsobaricQuantification` (off-limits) and `FLASHIda/Quantification` are unrelated.
  `channel_normalization` is applied in `Quantification::evaluate()` only, after the empty check,
  and never written into the logged channels (ADR-0057).

### Ownership and enqueueing

- **`FLASHIda` owns ten subsystems by value, and declaration order is load-bearing**: several bind
  references to earlier siblings, so reordering members is undefined behaviour.
- **`FLASHIda.cpp` is the sole enqueuer.** Builders (`buildMS2`/`buildMS3`/`buildFollowUp`) return
  a `ScanCommand`; they never push or register pending. Every push site is in the table
  `FLASHIda_Replan_test::push_site_inventory` pins: a new site fails it until named.
- **A survey's MS2s are a plan until the drain commits them** (ADR-0060). Selection journals stamps
  and rolls memory back; `replacePlan` swaps the plan per species; `dequeue_` commits under
  `queue_mutex_` and hands stamps to the next survey's fold. A committed species is never dropped.
  `wouldBar` must mirror `admitCandidate` (a drift-guard test pins them). Depth-1 harnesses never
  re-plan; `FLASHIda_Replan_test` forces the race through a seam.
- **An MS2's MS3s are pushed by a search worker's apply**, not by `processScan` (ADR-0062): the
  "MS3 targeting" and "exploration MS2 return" push sites live in `finishRegularMS2_` and
  `finishVariantMS2_`. The pool exists only with an `identification` section, and `~FLASHIda` stops
  it before any member an apply touches is destroyed.

## `processScan`

- **Five admission gates** before any branching (description < 3 chars; AGC echo; id never
  emitted; resolve race; unsupported context). Only an id the engine minted and sent is analysed.
  Gate 3 also refuses a result older than `MAX_IN_FLIGHT` (10 000) sent commands (`[TRACK-LOST]`).
- **Retention time is SECONDS here** (ADR-0051): the export converts the host's minutes; C++ tests
  pass seconds. `scan_results` `rt` and `ida.log` `RT` write minutes back.
- **Scan description**: `{3-char base-94 id}{marker}{payload}`, at most 15 chars. The marker is
  decoded only in `FLASHIda/ScanRole.h` (`S A M R Q E`; `C` retired with tag targeting, ADR-0061); its table is the source of the golden
  `scan_type` strings, and a new enumerator without a row fails to compile. **An undecodable
  description resolves as a survey**: do not "fail closed" there.
- **The decision gates read the engine's own command (`parent_ctx`)**; only the observe gate reads
  the echo (ADR-0042).
- **A matched MS2 is searched beside the analysis thread** (ADR-0062). `processScan` copies what
  the match needs (regular MS2, or an MS2 exploration variant) and returns; a `SearchPool` worker
  (`identification.search_threads`) runs `Exploration::searchNextLevel`/`searchVariant` with no
  lock, then `applyRegularSearch_`/`applyVariantSearch_` under `analysis_mutex_`. MS3 variants,
  trap variants and `Q` scans stay inline. A sweep completes when its searches return; searched
  rows are written in completion order. A full backlog (4 per worker, provisional) drops the oldest
  waiting search (`[SEARCH-DROPPED]`): its scan identifies nothing.
- **Quantification** (ADR-0038/0039/0056): `Q` is the scan measured; a verdict buys an ordinary `R`.
  `quantification.identify` decides which verdicts buy; `enriched_in` is evaluated on
  `condition_means`, **never on `fold_change`** (−1 is a sentinel for an absent condition).

## `getNextScanCommand` — never returns 0

Steps: (1) scheduled AGC prescan → (2) cycle-time MS1 → (2b) monitor scan → (3) cleanup →
(4) priority dequeue → (5) idle survey. Every path returns 1.

- **Step 5 emits no prescan** and re-enters `dequeue()` (ADR-0031). Prescans come only from Step 1.
- **Step 2b is not an `else` of Step 2**: `takeMonitorDue()` must run on every drain to re-arm.
- **`recordMS1Time()` carries a `roleDecides` conjunct**, or a monitor scan starves `cycle_time`.
- **No `analysis_mutex_` in the drain** (ADR-0025). It has two acquisition sites, `processScan` and
  a search worker's apply (ADR-0062). Lock order: `analysis_mutex_` → { `queue_mutex_`,
  `precursor_map_mutex_`, the `IdaLogger` stream mutexes, the `SearchPool` mutex }, all leaves.
  The steps are not atomic; the drain and `processScan` run on different threads.
- **Priorities**: 0 = CV-transition/cycle-time MS1, bought `R`, monitor `M`;
  1 = MS3; 2 = MS2; **3 = the idle survey and nothing else** (C# drain loops stop on it). Strict
  FIFO within a priority.

## Deconvolution: four call sites, three engines

MS1 inside `PrecursorSelection::filterAndRank`; regular MSn in `processScan`; exploration variants
through `Exploration`'s own `Deconvolution` (exploration tolerances); monitor scans through the
private `monitor_deconv_` (ADR-0042). **Do not merge the monitor engine into `deconv_`.**

## Characterization / MS3

- **`ProteoformTracker::planNextScans` is the only producer of MS3 targets.** No identified model
  (or `exhaustive`, ADR-0055) ⇒ zero MS3, by design.
- **`identification` gates matching; its `characterization` block gates MS3** (ADR-0061).
  `Config::identification()` (`IdentificationConfig`: `sequence`, `method`, `refine`,
  `has_refine_block`, `flashtnt`, `enabled()` = non-empty sequence) replaces
  `characterization().protein_sequence`/`.method`/`.refine` and `targeting()`'s seven FLASHTnT
  fields; the three MS2-matching gates in `FLASHIda.cpp` read `identification().enabled()`. Root
  `characterization`/`flashtnt` and the moved keys throw with a migration message.
- **`identification.characterization.mode`** is the objective (`ambiguity | coverage | fragments |
  exhaustive | discriminating`), required when the block is present; `"off"` and unknown values
  throw. `CharacterizationMode::Off` survives only as the internal no-block state. No mode is a set cover (ADR-0052);
  full-length ions are dropped where the pools are built. The engine branches on `.objective`, so a
  new mode without an objective ships as `ambiguity`.
- **`Config::applyCharacterizationMode_` must assign level 1 too**, or every MS1 is skipped. Budget
  and fragment-charge floor are read off level 2.
- **`fragmentContains` targets; `fragmentBrackets` narrows.** Never swap them (ADR-0005).
- **An empty MS3 clears the stored deconvolution** (`else deconv_.clearStoredMS2()`); no unit test
  pins it. Do not remove it as dead code.
- **Any activation feeds an MS3; only a b/y precursor's MS3 is annotated** (ADR-0054):
  `isKnownIonClass` (has a frame) and `isAnnotatableMs3Precursor` (b/y) answer different questions;
  do not merge them.
- **`identification.method`**: `getTopFragmentMatches` is the one switch between FLASHTnT and
  Refine (ADR-0049). `Refine` is pure static functions; `RefineSettings.h` exists to break an include
  cycle. The tracker owns the discriminating round; the MS3 return closes it; `Config::validate`
  refuses `refine`/`discriminating` under `flashtnt` and `discriminating` with level-3 exploration.
  Markers are stdout only.
- **A sweep is exactly its grid** (index arithmetic, top edge appended) and every pre-scan carries
  its **precursor share** (ADR-0048). The reaction-time floor (0.03 ms) still applies to swept ETD.
- **`kinetic_model`**: no variant wins; the follow-up is the fitted modeled maximum at every level
  (ADR-0047). The fit is the study's **lattice search**, not Levenberg-Marquardt. `k2` is never
  defaulted; one CE-swept activation per kinetic sweep; a failed fit aborts with no fallback.
- **Only `fragment_count` identifies MS3 pre-scans** (ADR-0020); a measuring metric at MS3 always
  emits a production scan. Use `isMeasuringMetric()`, never `== FragmentCount`.
- **A trap pre-scan is measured, never identified** (ADR-0045): decided from the command's own
  `analyzer`; only `remaining_precursor` may score a trap sweep.

## Config (`FLASHIda/Config.cpp`)

- Unknown keys throw via `rejectUnknownKeys` (~18 hand-written allowed-sets).
- **`scheduling.target_depth` is allowed and never read, correctly**: only the host sizes the
  instrument queue. Do not wire it up.
- **`.value(key, default)` fallbacks are dead in production**; don't trust a C++ default.
- **`kScanKeys` (17 keys) both validates and parses every scan object**; the binding constraint is
  the C# emit set (ADR-0011).
- `ScanConfig.analyzer` is a closed set (`Orbitrap | IonTrap | ""`), checked on the raw JSON.
- `applyOverrides` ignores unknown keys and needs string values. `deconvolution.tol` needs ≥ 3.

## Logging (`FLASHIda/IdaLogger.cpp`)

Five streams in `runtime.log_dir` (empty = nothing opens; the host creates the directory):
`ida.log` (free text, decision record), `scan_commands.tsv` (36 columns, one row per dequeued
command; what FLASHDeconv reads), `scan_results.tsv` (56), `identification.tsv` (46),
`pooled_identification.tsv` (20). `IdaLogger` serialises each stream itself.

- **Every value written must set its own precision**: the ostream is sticky. `Mass=` and `AllMass=`
  are asserted identical (ADR-0035); `scan_commands` `mono_mass` has four decimals, and no golden can
  see either (the comparer's `RelTol 1e-3`).
- `ida.log`'s `MS1 Scan#` is the instrument scan number; `Access ID` is the tracking id. Every
  `Mass=` line ends in `ID=`; a re-plan writes `Superseded=` (ADR-0060). Readers anchor on
  `" 0 targets"` with the leading space.
- Every TSV id cell goes through `quoteLogCell` (ids starting with `"` are RFC 4180-quoted).
- An `ida.log` target is a mass active within `rt_window` of a reference acquisition, with the old
  behaviour behind `log_target_rt_window` / `log_target_charges`.
- Rows are written at dequeue. Two-stage MS3 scalars are `"stage0;stage1"`; stage-less rows write
  `"0"`. Notches add `','` inside the per-stage columns. **Id lists join with a space**, because the
  id alphabet contains `';'`.
- A column reorder needs no recapture (goldens compare by name) but does need
  `FLASHIda_LoggingFields_test`'s hard-coded indices updated.

## Offline coupling (ADR-0046)

`FLASHIda/ScanCommandJoin.h` runs offline inside FLASHDeconv: a spectrum's tracking id finds its
`scan_commands.tsv` row, `decision_survey_id` (else the parent chain) finds the survey, and
FLASHDeconv takes its **own** PeakGroup for that mass from that survey. The join fails closed
(foreign file, duplicate id, nothing joined). MS2 only. `parse` resolves columns by name and must
keep reading pre-2026-09 run folders. The FLASHDeconv side is inside the boundary.

## Co-isolation notches (ADR-0016, ADR-0019)

- **`num_stages` counts cascade stages, never notches.**
- Notches live in `Notch notches[18]`, stage k owning a fixed block of 9; use `notchesForStage()`.
  `MAX_ISOLATION_STAGES` and `MAX_NOTCHES_PER_STAGE + 1` are **different tens**.
- `NotchSelection.h` holds the policy (SNR gate, authored `charge_selection`, clamp). Under `snr` and
  `qscore` the anchor must be a **resolved** charge (`463d2ba1a4`).
- `snr_threshold`, `max_charge_states` (→ `max_charge_states - 1` notches) and `min_charge_states`
  (exempts a sibling) are authored and inert by default (ADR-0040).
- `separate` and `multiplexed` read **one** `selectNotches` call and differ only in scan count
  (ADR-0021). **Fan out at the emit loop, never in the candidate loop.**
- An inclusion row's authored charge set restricts, never extends, and re-keys exclusion to
  `(nominal mass, charge)` for those species only (ADR-0028).
- MS3 needs the fragment's whole envelope (`PeakRecord::by_charge`).

## Tests

- A C++ test runs only if it is in **both** the yml's `--target` list and its `ctest -R`
  alternation; `executables.cmake` alone gets you nothing.
- Fixtures load from `../../FlashIDA/test-data` relative to the build dir; test exes land in
  `<build>/src/tests/class_tests/bin/` and need the 5-DLL set (incl. `zlib.dll`) on Windows.
- **Drive acquisitions only through `FLASHIda_TestHelpers.h::runInterleaved`**, the mirror of C#
  `PushScanAndDrainFull`; a hand-rolled drive loop is a defect even if correct.
- **Every drive calls `waitForSearchIdle()` after every `processScan`** (ADR-0062), so it sees the
  order a synchronous `processScan` produced; the bridge's sixth export `WaitForSearchIdle` exists
  for that and production never calls it. Never call it holding `analysis_mutex_`. Only the
  sections that pin the pool skip the wait (`waitUntil`, `SearchGate`, the `FLASHIda_TestAccess.h`
  search seams).
- C++ ctests assert **plausibility ranges**; C# NUnit asserts exact goldens.
- Release only: `OPENMS_PRECONDITION` and asserts are compiled out, and the TOPDOWN tree has none
  anyway.
