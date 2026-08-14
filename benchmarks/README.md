# Scientific API benchmark suite

These benchmarks are public for scientific reproducibility and project
maintenance, but they are not part of the public API contract, general user
documentation, or user-facing performance guarantees.

The focused CHECKED/UNCHECKED protocol is documented alongside the benchmark
sources in [Validation-sensitive scientific benchmarks](validation-sensitive.md).

`mmcfilters_api_benchmark` measures representative public API paths while
preserving the same mathematical workload in `CHECKED` and `UNCHECKED` builds.
It complements the smaller validation-sensitive benchmark with construction,
attributes, filters, interoperability, editing, CASF, and complete scientific
pipelines.

## Measurement scopes

Every scenario states one of two scopes in its name and output record:

- `end_to_end` includes construction of the data structure and all work needed
  to produce the scientific result;
- `established_input` prepares a valid tree or mutable state outside the timer
  and measures one public operation over that established input.

An untimed warm-up precedes the samples. Checksumming and destruction of the
result are outside the timed interval. The output contains median, median
absolute deviation, minimum, maximum, workload dimensions, and a deterministic
checksum. `--emit-samples` additionally emits every timed repetition; it does
not alter the measured operation.

## Profiles and inputs

| Profile | Default image | Repetitions | Intended use |
| --- | ---: | ---: | --- |
| `smoke` | 48 x 48 | 2 | Fast correctness and integration check |
| `core` | 192 x 192 | 5 | Routine API performance evaluation |
| `publication` | 512 x 512 | 15 | Broad experiment for a reported result |

The `structured`, `noise`, `ramp`, `geometric`, and `flat` deterministic inputs
can be selected independently. `--rows`, `--columns`, and `--repetitions` override
the profile defaults. `--input-file` accepts a real image, converts it to one
grayscale channel, and verifies the declared dimensions without implicit
resampling.

## Reproducible workload manifests

The versioned manifest at
`benchmarks/workloads/scientific-api.ini` defines named workloads. Each entry
can pin the profile, synthetic generator or image file, dimensions, repetitions,
suites, CASF quantiles, input checksum, and scientific attribute bundles.

```text
./build-benchmark-checked/benchmarks/mmcfilters_api_benchmark \
  --manifest benchmarks/workloads/scientific-api.ini \
  --workload core_structured
```

Real-image paths are resolved relative to the manifest. File workloads must
declare `rows` and `columns`; a mismatch aborts instead of resizing the data. The
reported `input_checksum` is platform-independent FNV-1a over the ASCII
`<rows>x<columns>:` header and grayscale bytes. An `input_checksum` field in the
manifest pins those bytes, so changing dimensions or image contents fails
before any timed scenario.

Bundle fields use scalar names and explicit group references:

```text
attribute_bundle.radiometric_size=AREA,MEAN_GRAY_LEVEL,VOLUME,GRAY_LEVEL_HEIGHT
attribute_bundle.shape_boundary=group:BOUNDARY,AREA,BOX_WIDTH,BOUNDING_BOX_HEIGHT,MAX_DIST
```

Every bundle is compared with its deduplicated sequence of scalar API calls.
Overlapping groups are therefore allowed, while semantic equivalence remains a
mandatory checksum condition.

The manifest also pins three real images already used by the repository's
numerical-validation tests: `publication_real_lena`,
`publication_real_brain`, and `publication_real_wrist`. Their bundles represent
actual API workloads for connected filtering, moment/shape descriptors,
`MAX_DIST`, area-based CASF, and bounding-box CASF. The names deliberately
describe algorithms rather than asserting undocumented article provenance.

## Suites

| Suite | Representative paths |
| --- | --- |
| `construction` | max-tree, min-tree, Tree of Shapes, residual trees, and additional altitude types |
| `attributes` | scalar attributes, heterogeneous requests, altitude-based node-attribute samples, public groups, and manifest-defined bundles |
| `filters` | direct, subtractive, pruning, Viterbi, preservation-mask adjustment, UAO, extinction, MSER, and depth stability |
| `interoperability` | reconstruction, Higra import/export and attribute projection, saliency maps |
| `editing` | safe single edits, staged editor commits, and editing sequences |
| `casf` | paired-tree construction, incremental sequences, a step after an established prefix, bounding-box attributes, and typed C++ paths |
| `pipelines` | tree construction, attribute computation, preservation-mask construction, and filtering in one timed scientific workflow |

### Attribute groups

The benchmark treats groups as first-class workloads: `GRAY_LEVEL`, `SHAPE`,
`MOMENTS`, `BOUNDARY`, `TREE_TOPOLOGY`, and `ALL`. In the `core` profile, each
group except `ALL` is compared with the equivalent sequence of public scalar
calls; `publication` also adds the sequential `ALL` reference. A canonical
checksum sorts columns by scalar attribute identity, so it verifies semantic
equivalence without depending on layout order. The benchmark aborts if the
grouped and sequential results differ.

This comparison separates two scientifically relevant questions: the cost of
requesting a coherent family of attributes in one dependency-aware operation,
and the cost paid by a caller that requests the same attributes independently.

### CASF

CASF has separate measurements for construction, construction plus a threshold
sequence, a sequence over an established CASF state, and light, medium, and
heavy steps with their preceding state prepared outside the timer. Thresholds
are empirical quantiles of the selected attribute over the non-root nodes of
both initial trees; the defaults are `0.10,0.50,0.90` and can be changed through
the manifest or `--casf-quantiles`. Quantile derivation is experimental setup
and remains outside the timed region; the reported pipeline time begins with
CASF construction and uses the already fixed thresholds.

Every CASF record reports the effective thresholds, nodes removed from each
tree (`primary` is the min-tree and `secondary` is the max-tree), and
incremental and complete validation-commit counts. A quantile step
with no structural edit is retained and reported as zero work instead of being
silently replaced. The complete sequence must exercise incremental editing
when candidates exist, and any complete-validation commit in the hot path
aborts the benchmark. The checksum includes the output image and full ownership
state of both component trees.

## Build and run

```text
cmake -S . -B build-benchmark-checked -DCMAKE_BUILD_TYPE=Release \
  -DMMCFILTERS_BUILD_PYTHON=OFF -DMMCFILTERS_BUILD_BENCHMARKS=ON \
  -DMMCFILTERS_CONTRACT_MODE=CHECKED
cmake --build build-benchmark-checked --target mmcfilters_api_benchmark

cmake -S . -B build-benchmark-unchecked -DCMAKE_BUILD_TYPE=Release \
  -DMMCFILTERS_BUILD_PYTHON=OFF -DMMCFILTERS_BUILD_BENCHMARKS=ON \
  -DMMCFILTERS_CONTRACT_MODE=UNCHECKED
cmake --build build-benchmark-unchecked --target mmcfilters_api_benchmark

./build-benchmark-checked/benchmarks/mmcfilters_api_benchmark \
  --manifest benchmarks/workloads/scientific-api.ini \
  --workload core_structured
```

Suites may be repeated or comma-separated. `--format jsonl` emits one metadata
record followed by one record per scenario.

Compare modes with alternating process order and mandatory checksum equality:

```text
python3 benchmarks/compare_validation_modes.py --process-runs 5 \
  --capture-samples \
  --output-dir benchmarks/results/publication-structured \
  --experiment-name publication-structured \
  build-benchmark-checked/benchmarks/mmcfilters_api_benchmark \
  build-benchmark-unchecked/benchmarks/mmcfilters_api_benchmark -- \
  --manifest benchmarks/workloads/scientific-api.ini \
  --workload publication_structured

python3 benchmarks/analyze_scientific_experiment.py \
  benchmarks/results/publication-structured
```

Everything after the standalone `--` is forwarded to both benchmark
executables. Without `--output-dir`, the comparison remains stdout-only for
compatibility. `--capture-samples` forwards `--emit-samples`, checks that the
raw repetitions reconstruct every reported median, MAD, minimum, and maximum,
and persists `samples.jsonl` and `samples.csv`. With `--output-dir`, the runner
also retains every raw process output and writes `experiment.json`, normalized
`raw.jsonl`, `summary.jsonl`, and `summary.csv`.
The experiment metadata includes Git state, CPU, operating system, executable
hashes, CMake caches, compiler flags, command line, and timestamps. The analysis
step adds `scope-summary.csv`, a checksum-paired `scope-pairs.csv`,
`attribute-pairs.csv`, `casf.csv`, `sample-statistics.csv`, and a concise
`analysis.md` report. The paired scope table compares the complete CASF area
pipeline with the same fixed threshold sequence over an already established
CASF state.

The inferential unit is one process invocation, not one timed repetition. For
each scenario and contract mode, the analysis first computes one median per
process and then uses a deterministic 20,000-replicate cluster bootstrap over
those process medians to estimate a 95% interval for the relative difference.
It also reports Cliff's delta over process medians; negative values favor
UNCHECKED and positive values favor CHECKED. A result is classified as faster
or slower only when the interval excludes zero. The pooled repetition median is
retained as a descriptive value, never used to inflate the independent sample
count.

Multiple completed experiment directories can be combined without mixing Git
commits, CPUs, profiles, repetition counts, or compilers:

```text
python3 benchmarks/summarize_scientific_campaign.py \
  benchmarks/results/publication-campaign \
  benchmarks/results/publication-structured \
  benchmarks/results/publication-real-lena \
  benchmarks/results/publication-real-brain \
  benchmarks/results/publication-real-wrist
```

The campaign report retains workload-level results, all observations, scenarios
common to every workload, timing-scope summaries and pairs, attribute grouping,
and CASF edit outcomes as separate CSV tables. When every source experiment has
captured samples, it also writes `campaign-sample-statistics.csv` and summarizes
supported, opposed, and inconclusive intervals without pooling repetitions
across workloads.

The legacy positional runner for
`mmcfilters_validation_sensitive_algorithms_benchmark` remains supported.
When tests and benchmarks are enabled together, CTest registers
`benchmark_api_manifest_smoke`; it checks the pinned input, bundle equivalence,
CASF determinism, and output generation without imposing a performance limit.
The Linux CI smoke job builds both compile-time contract modes, runs that CTest
in each build, and alternates two process executions to require equality of all
scenario checksums and structural outcomes.

## Python facade

`benchmarks/python_api_benchmark.py` measures the public Python binding as a
separate layer, including tree construction, scalar and grouped attributes,
filtering, CASF, and an end-to-end pipeline. It also verifies grouped results
against sequential scalar calls in `core` and `publication`.

```text
PYTHONPATH=build/python python3 benchmarks/python_api_benchmark.py \
  --profile core --input structured --casf-quantiles 0.10,0.50,0.90 \
  --format jsonl
```

Use the same Python interpreter ABI used to compile the extension. Python
results should not be mixed with the C++ timings: they answer the distinct
question of public binding and conversion overhead.

## Reporting rule

For an article, retain the compiler and flags, CPU, operating system, contract
mode, profile, input generator, dimensions, repetitions, process-run count,
all checksums, `raw.jsonl`, `samples.jsonl`, and `samples.csv`. Interpret timing
differences only after the `CHECKED` and `UNCHECKED` checksums match for every
scenario, and distinguish descriptive point estimates from process-level
confidence intervals.
