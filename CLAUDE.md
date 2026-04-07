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

## FLASHIda Development (Phase 8 — Cleanup + Documentation)

This fork's `flashida-v9-bridge` branch is the active C++ development branch for FLASHIda. All C++ changes go here. The `build-dlls` workflow auto-triggers on push and produces `OpenMS.dll` for the C# side.

### Key files for Phase 8

All in `src/openms/{include,source}/OpenMS/ANALYSIS/TOPDOWN/`:
- **FLASHIdaBridgeFunctions.h** — Remove 18 `extern "C" OPENMS_DLLAPI` declarations. Leave exactly 5: `CreateFLASHIda`, `DisposeFLASHIda`, `ProcessScan` (with `double faims_cv`), `GetNextScanCommand`, `GetNextTrackingId`.
- **FLASHIdaBridgeFunctions.cpp** — Remove 18 function bodies. Leave 5 keepers only.
- **FLASHIda.h** — Remove orphaned private method declarations (`getPeakGroupSize_`, `fillIsolationWindows_`, `deconvolveMS2_`, `resetScanState_`, `parseLegacy_`). **Do NOT remove** Phase 7 exploration methods (16 production + 6 ForTest helpers — see Phase 8 plan Step 4).
- **FLASHIda.cpp** — Remove `parseLegacy_()` definition + `else parseLegacy_(arg)` branch in constructor. Remove orphaned private method bodies. Throw `std::invalid_argument` for non-JSON input.
- **OptimizationMetadata.h** — 19-field struct (field 19 `exploration_metric` added in Phase 7). No changes in Phase 8.

### Test files

- **Existing (6 binaries, must all remain in CI):** `DeconvolvedSpectrum_OptimizationMetadata_test`, `FLASHIdaQueueTracking_test`, `FLASHIda_ProcessScan_test`, `ScanCommandLayout_test`, `FLASHIdaFAIMS_test`, `FLASHIda_exploration_test`
- **New (Phase 8):** `src/tests/class_tests/openms/source/FLASHIda_LegacyConfig_test.cpp` — P8-U04 (legacy config rejection)
- **Register in:** `src/tests/class_tests/openms/executables.cmake` — add `FLASHIda_LegacyConfig_test`

### Critical constraints

- `ScanCommand` is 1248 bytes — Phase 8 does NOT modify the struct
- MSVC `/WX` treats warnings as errors — removing bridge functions may expose unused parameters (`C4100`) or unused variables (`C4189`); fix with `(void)var;` cast, never by commenting out singleton initializers
- `parseJSONConfig_()` has `std::stringstream ss` at line ~3198 — use descriptive variable names, not `ss`/`it`/`cfg`
- Test helpers follow `nameForTest()` pattern
- `DeconvolvedSpectrum(int scan_number)` — constructor takes scan_number, NOT ms_level

### Build #4 (final)

Phase 7 C++ already shipped. Phase 8 C++ changes (dead code removal, legacy config rejection) are the remaining batch. Push once to `flashida-v9-bridge`, wait ~40 min for DLL build. Verify GCC compilation locally before pushing (MSVC has stricter warnings).
