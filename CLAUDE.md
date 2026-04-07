# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Important Warning

**DO NOT build OpenMS unless explicitly asked.** The build is extremely resource-intensive and time-consuming.

## Project Overview

OpenMS is an open-source C++20 library for LC-MS data management and analysis (proteomics, metabolomics). This is a development fork (`FIdevelop` branch) with FLASH/FLASHIda extensions for top-down proteomics — deconvolution, tagging, and tandem mass tag algorithms.

## Build System

CMake-based (minimum 3.21). Dependencies managed via vcpkg.

```bash
# Configure (from a build directory)
cmake -DCMAKE_PREFIX_PATH=<vcpkg-installed-path> <source-dir>

# Build
cmake --build <build-dir> --config Release

# Key CMake options:
#   -DPYOPENMS=ON          Build Python bindings
#   -DWITH_GUI=OFF          Skip GUI (TOPPView) if Qt not needed
#   -DMT_ENABLE_OPENMP=ON   OpenMP parallelization (default ON)
#   -DWITH_HDF5=ON          HDF5 support
#   -DWITH_PARQUET=ON       Parquet support
#   -DENABLE_STYLE_TESTING=ON  cpplint/cppcheck checks
```

Default build type is `Release`. Uses C++20 standard.

## Testing

Tests use CTest with a custom OpenMS ClassTest framework (`<OpenMS/CONCEPT/ClassTest.h>`).

```bash
# Run all tests
ctest

# Run a single test by name
ctest -R ClassName_test

# Run tests matching a pattern
ctest -R FLASH
ctest -R pyopenms
```

### Test types and locations
- **Class tests** (unit): `src/tests/class_tests/openms/source/ClassName_test.cpp` — registered in `src/tests/class_tests/openms/executables.cmake`
- **TOPP tests** (integration): `src/tests/topp/` — tests CLI tools with input/output comparison via FuzzyDiff
- **GUI tests**: `src/tests/class_tests/openms_gui/`
- **Python tests**: `src/pyOpenMS/tests/` (unittests, integration, memoryleak)

To add a new class test: create `ClassName_test.cpp` in `src/tests/class_tests/openms/source/` and add `ClassName_test` to the appropriate list in `executables.cmake`.

Note: FLASH test entries (`FLASHDeconvAlgorithm_test`, `FLASHDeconvHelperStructs_test`, `FLASHTaggerAlgorithm_test`) are currently commented out in `executables.cmake`.

## Code Architecture

### Layered structure
1. **Data Structures** (KERNEL): `Peak1D`, `MSSpectrum`, `MSExperiment`, `Feature`, `FeatureMap`, `ConsensusMap`, `PeptideIdentification`
2. **Algorithms** (PROCESSING, ANALYSIS): signal processing, feature detection, quantification, identification
3. **TOPP Tools** (`src/topp/`): ~150 CLI tools subclassing `TOPPBase` with standardized parameter handling (CTD scheme)
4. **pyOpenMS** (`src/pyOpenMS/`): Python bindings via Cython/autowrap

### Source layout
- Headers: `src/openms/include/OpenMS/<MODULE>/ClassName.h`
- Sources: `src/openms/source/<MODULE>/ClassName.cpp`
- Modules map to subdirectories: KERNEL, FORMAT, PROCESSING, ANALYSIS, CHEMISTRY, MATH, FEATUREFINDER, COMPARISON, ML, SYSTEM, CONCEPT, etc.

### FLASH-specific code (this fork)
Located in `ANALYSIS/TOPDOWN/`:
- `FLASHDeconvAlgorithm` — deconvolution
- `FLASHTaggerAlgorithm` — amino acid tagging
- `FLASHExtenderAlgorithm` — sequence extension
- `FLASHGappedTaggerAlgorithm` — gapped tagging (new, empty files)
- `FLASHTnTAlgorithm` — tandem mass tags
- `FLASHIda` / `FLASHIdaBridgeFunctions` — IDA integration
- Helper classes: `DeconvolvedSpectrum`, `PeakGroup`, `MassFeatureTrace`, `PeakGroupScoring`

### pyOpenMS binding workflow
1. Declare C++ classes in `.pxd` files (`src/pyOpenMS/pxds/`)
2. Autowrap generates Cython code; manual extensions go in `src/pyOpenMS/addons/`
3. Cython compiles `.pyx` to C++, then compiled into a Python module
4. Annotations: `wrap-ignore`, `wrap-as`, `wrap-iter-begin` for customization

### TOPP tool development
Subclass `TOPPBase`, implement `main_()`. Tools get: parameter handling, logging, progress monitoring, I/O handling. Registration via `install_tool()` macro in CMake.

## Code Style

- clang-format configured (`.clang-format`): LLVM-based, 150 column limit, 2-space indent, Allman brace style, `SpaceAfterLogicalNot: true`
- Naming: PascalCase for classes, camelCase for variables
- Test file naming: `ClassName_test.cpp`
- All files need SPDX BSD-3-Clause license header with `$Maintainer` and `$Authors` tags

## Key Dependencies

Boost (multiple components), Qt6, Xerces-C (XML), Eigen3 (linear algebra), SQLite, COIN-OR (LP solvers), LibSVM. Bundled: nlohmann_json, IsoSpec, SQLiteCpp.

## CI

GitHub Actions in `.github/workflows/`. Main pipeline: `openms_ci_matrix_full.yml` (Windows MSVC, macOS Intel/ARM, Ubuntu GCC/Clang). CDash dashboard at cdash.seqan.de.

## FLASHIda Development (Phase 7 — Exploration Engine)

This fork's `flashida-v9-bridge` branch is the active C++ development branch for FLASHIda. All C++ changes go here. The `build-dlls` workflow auto-triggers on push and produces `OpenMS.dll` for the C# side.

### Key files for Phase 7

All in `src/openms/{include,source}/OpenMS/ANALYSIS/TOPDOWN/`:
- **FLASHIda.h** — Main header. Contains `ScanCommand` struct (1248 bytes, `static_assert` at line ~109), enums, and class definition. Phase 7 adds: `SelectionMetric`, `ExplorationMetric` enums, `MSLevelConfig` struct, `ExplorationGroup`/`ExplorationVariant` structs, `level_configs_` map, `variant_tracking_to_group_` map.
- **FLASHIda.cpp** — Implementation (~3800 lines). `processScan()` is the main entry point (MS1/MS2 routing). Phase 7 adds: `parseLevelConfig_()`, `initiateExploration_()`, `feedExplorationResult_()`, `computeExplorationScore_()`, `initiateNextLevel_()`, `buildCEVariants_()`, `applyOverrides_()`.
- **FLASHIdaBridgeFunctions.cpp/.h** — C bridge exports consumed by C# P/Invoke. Phase 7 does not add new exports.
- **OptimizationMetadata.h** — 18-field struct for recording exploration outcomes. Already exists from Phase 2. `DeconvolvedSpectrum` carries `optional<OptimizationMetadata>`.

### Test file for Phase 7

- New: `src/tests/class_tests/openms/source/FLASHIda_exploration_test.cpp`
- Register in: `src/tests/class_tests/openms/executables.cmake`
- Existing test binaries (must all remain in CI): `DeconvolvedSpectrum_OptimizationMetadata_test`, `FLASHIdaQueueTracking_test`, `FLASHIda_ProcessScan_test`, `ScanCommandLayout_test`, `FLASHIdaFAIMS_test`

### Critical constraints

- `IsolationStage.collision_energy` is `double` — all CE variables must be `double`, not `int`
- `IsolationStage.activation_type` is `char[32]` — use 32-byte buffers
- `ScanCommand` is 1248 bytes — do not change without the 6-file lockstep rule
- MSVC `/WX` treats warnings as errors — use `(void)var;` for unused-in-release variables
- `DeconvolvedSpectrum(int scan_number)` — constructor takes scan_number, NOT ms_level
- `toSpectrum()` returns `MSSpectrum` by value and requires at least one PeakGroup
- Test helpers follow `nameForTest()` pattern: `encodeBase36ForTest()`, `pushCommandForTest()`, `decodeBase36ForTest()`, `updateCVSkipForTest()`, `getCVSkipAmountForTest()`

### Build batching

Phase 7 + Phase 8 = Build #4 (final). Batch all changes before pushing to `flashida-v9-bridge`. Each push triggers a ~40 min DLL build. Verify GCC compilation locally before pushing (MSVC has stricter warnings).
