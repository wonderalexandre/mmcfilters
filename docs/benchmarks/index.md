# Benchmark artifact index

This index separates reproducible measurements from historical evidence whose
input or code provenance is incomplete. Timing values are machine-specific and
must not be interpreted as portable performance guarantees.

## Artifact catalog

| Artifact | Code snapshot | Protocol | Provenance status |
| --- | --- | --- | --- |
| [Two-tree residual, pre-cleanup](min-max-residual-icdar-480p-first10-2026-08-04.md) | base `ff913c1` plus an uncommitted residual-tree integration | first 10 480p images, 3 alternating repetitions, median | historical only; exact source tree cannot be reconstructed from Git |
| [Two-tree residual, post-cleanup](min-max-residual-icdar-480p-first10-2026-08-04-post-cleanup.md) | results committed with `e70c51f` | first 10 480p images, warm-up plus 3 alternating repetitions, median | code/results recorded; dataset provenance incomplete |
| [Six-family comparison, pre-cleanup HTML](tree-construction-comparison-icdar-first10-2026-08-04.html) | results committed with `e70c51f`; run predates cleanup | 10 images at 480p/720p/1080p, 5 rotated repetitions | code/results recorded; dataset provenance incomplete |
| [Six-family comparison, post-cleanup HTML](tree-construction-comparison-icdar-first10-2026-08-04-post-cleanup.html) | `e70c51f` | 10 images at 480p/720p/1080p, warm-up plus 5 rotated repetitions | code/results recorded; dataset provenance incomplete |

The JSON files beside the HTML reports are portable report manifests. The raw,
per-image, summary, scaling, paired-ratio, and validation data are under
`benchmarks/results/` with matching stems.

## Recorded environment

- Date of the recorded runs: 4 August 2026.
- Hardware: Apple M4.
- Compiler: Apple Clang 21 where recorded.
- Build: Release/NDEBUG, sequential benchmark execution.
- Connectivity: 8-neighbourhood for the six-family comparison; shared
  4-neighbourhood for the two-tree residual report.
- Timing scope after cleanup: external `std::chrono::steady_clock` around the
  complete public construction call. Image loading and tree destruction are not
  included.

The operating-system version, memory capacity, CPU power mode, exact CMake
version, and full compiler build string were not captured. New benchmark runs
must record them.

## Build and execution commands

Configure and build the benchmark executables from the repository root:

```bash
cmake -S . -B build-bench \
  -DMMCFILTERS_BUILD_PYTHON=OFF \
  -DMMCFILTERS_BUILD_BENCHMARKS=ON \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build-bench --parallel --target \
  mmcfilters_min_max_residual_tree_benchmark \
  mmcfilters_tree_construction_comparison_benchmark
```

Run the complete six-family comparison:

```bash
./build-bench/benchmarks/mmcfilters_tree_construction_comparison_benchmark \
  /path/to/icdar-root \
  benchmarks/results/tree-construction-comparison-new-raw.csv \
  5
```

The data root must contain `icdar_480p/`, `icdar_720p/`, and `icdar_1080p/`,
each with `test_000.png` through `test_009.png`. Run the focused residual
benchmark once per image:

```bash
./build-bench/benchmarks/mmcfilters_min_max_residual_tree_benchmark \
  /path/to/icdar-root/icdar_480p/test_000.png \
  3
```

Generate the comparison summaries and portable report with the checked-in
analysis programs; use `--help` to inspect their current argument contract:

```bash
python benchmarks/analyze_tree_construction_comparison.py --help
python benchmarks/build_tree_construction_report.py --help
```

## Dataset provenance gap

The recorded files identify their inputs only as an "ICDAR" set with generated
480p, 720p, and 1080p directories. The original dataset title/version, download
URL, redistribution license, resize/color-conversion procedure, and per-image
checksums were not committed. Several ICDAR datasets use overlapping naming, so
this repository cannot infer those facts safely.

Until that information is supplied, these reports are valid historical local
measurements but are not independently reproducible dataset artifacts. A new
run must add a manifest containing, at minimum:

```text
dataset_name
dataset_version
source_url
license_name_and_url
preprocessing_command
file_path,sha256,rows,cols,dtype
```

Generate the checksum column only from the actual benchmark inputs, for example
with `shasum -a 256` or an equivalent SHA-256 tool. Do not copy checksums from a
different resolution or dataset release.

## Acceptance checklist for new results

- Commit hash and clean/dirty worktree state are recorded.
- Complete compiler, CMake, operating system, and hardware identifiers are
  recorded.
- Dataset URL, license, preprocessing, and SHA-256 manifest are present.
- Warm-up, repetition order, sample count, aggregation, and timing boundaries
  are explicit.
- Exact reconstruction, deterministic node counts, positive timings, and raw
  row completeness pass before summaries are generated.
- Raw data and the commands that produced every derived table remain available.
