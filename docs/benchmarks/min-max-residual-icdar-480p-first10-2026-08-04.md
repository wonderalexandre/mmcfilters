# Two-tree residual construction on the first ten ICDAR images

Date: 4 August 2026.

Provenance status: historical only. The exact code was `ff913c1` plus
uncommitted changes and therefore cannot be reconstructed from Git. Dataset
source, license, preprocessing, and checksums were not recorded. See the
[benchmark artifact index](index.md).

> [!NOTE]
> This is the historical **pre-cleanup** run. Its two phase-level timing
> columns were collected by instrumentation inside the production builder and
> are retained only for reproducibility. The current builder no longer records
> internal durations. See the
> [post-cleanup external benchmark](min-max-residual-icdar-480p-first10-2026-08-04-post-cleanup.md)
> for measurements of the current implementation.

## Purpose

This benchmark validates the unrestricted and saturated synchronized
max-tree/min-tree implementations after their integration into
`MorphologicalAttributeFilters`. Both modes use the same agenda, incremental
dual-tree updater and residual assembler. They differ only in the eligibility
policy applied to current regional extrema.

## Protocol

- images `test_000.png` to `test_009.png` from the ICDAR 480p set;
- 480 x 853 pixels, 8-bit grey scale;
- shared 4-adjacency (`radius=1`);
- `ContrastInvariantSpatial` tie policy;
- exterior pixel 0;
- three alternating repetitions per mode and image;
- median of the three repetitions;
- total construction time, including the initial max-tree and min-tree;
- Release/NDEBUG build, Apple Clang 21, Apple M4;
- local branch `saliency-map`, base commit `ff913c1` plus the uncommitted
  residual-tree integration.

The executable was built with:

```bash
cmake -S . -B build-bench \
  -DMMCFILTERS_BUILD_PYTHON=OFF \
  -DMMCFILTERS_BUILD_BENCHMARKS=ON \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build-bench \
  --target mmcfilters_min_max_residual_tree_benchmark -j2
```

## Results

| Image | Unrestricted | Saturated | Ratio |
|---|---:|---:|---:|
| `test_000` | 626.025 ms | 942.959 ms | 1.506x |
| `test_001` | 528.781 ms | 677.754 ms | 1.282x |
| `test_002` | 383.921 ms | 560.045 ms | 1.459x |
| `test_003` | 506.499 ms | 1081.153 ms | 2.135x |
| `test_004` | 484.018 ms | 648.897 ms | 1.341x |
| `test_005` | 378.572 ms | 534.652 ms | 1.412x |
| `test_006` | 279.352 ms | 403.124 ms | 1.443x |
| `test_007` | 462.333 ms | 775.116 ms | 1.677x |
| `test_008` | 259.475 ms | 361.694 ms | 1.394x |
| `test_009` | 407.732 ms | 617.198 ms | 1.514x |

Aggregates:

| Measure | Unrestricted | Saturated |
|---|---:|---:|
| Sum of per-image medians | 4.317 s | 6.603 s |
| Mean per image | 431.671 ms | 660.259 ms |
| Median across images | 435.033 ms | 633.048 ms |
| Mean candidate preparation/certification time | 59.037 ms | 146.892 ms |
| Share of total time recorded in that stage | 13.68% | 22.25% |

The final two rows above describe the removed internal instrumentation. They
must not be interpreted as metrics exposed by the current builder or compared
with externally sampled construction phases.

The saturated construction adds 228.588 ms per image on average, or 52.95%
over the unrestricted construction. Its aggregate ratio is 1.530x. The
per-image ratios range from 1.282x (`test_001`) to 2.135x (`test_003`).

Both modes reconstructed every input exactly during the benchmark warm-up.
The unrestricted mode rejected no extrema and performed no complement
traversal. The complete measurements are available in
`benchmarks/results/min-max-residual-icdar-480p-first10-2026-08-04.csv`.
