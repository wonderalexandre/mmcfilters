# Two-tree residual construction after profiling cleanup

Date: 4 August 2026.

Code/results snapshot: `e70c51f`. Dataset source, license, preprocessing, and
checksums were not recorded; see the [benchmark artifact index](index.md).

## Purpose

This benchmark measures the current unrestricted and saturated synchronized
max-tree/min-tree implementations after removal of production timing
instrumentation and non-essential work counters. It verifies that performance
measurement is external to the builder while hierarchy construction remains
unchanged.

## Protocol

- images `test_000.png` to `test_009.png` from the ICDAR 480p set;
- 480 x 853 pixels, 8-bit grey scale;
- shared 4-adjacency (`radius=1`);
- `ContrastInvariantSpatial` tie policy;
- exterior pixel 0;
- one untimed warm-up per mode and image;
- three alternating repetitions per mode and image;
- median of the three repetitions;
- external `std::chrono::steady_clock` around the complete construction call;
- total construction time, including the initial max-tree and min-tree;
- Release/NDEBUG build, Apple M4.

The production builder does not sample clocks. Its retained statistics are
semantic diagnostics used by consistency checks and regression tests, not a
phase profiler.

## Results

| Image | Unrestricted | Saturated | Ratio |
|---|---:|---:|---:|
| `test_000` | 577.742 ms | 868.461 ms | 1.503x |
| `test_001` | 474.768 ms | 628.675 ms | 1.324x |
| `test_002` | 345.412 ms | 531.393 ms | 1.538x |
| `test_003` | 475.291 ms | 1011.167 ms | 2.127x |
| `test_004` | 439.468 ms | 602.846 ms | 1.372x |
| `test_005` | 341.986 ms | 479.981 ms | 1.404x |
| `test_006` | 252.543 ms | 365.027 ms | 1.445x |
| `test_007` | 414.691 ms | 710.786 ms | 1.714x |
| `test_008` | 238.148 ms | 327.966 ms | 1.377x |
| `test_009` | 357.749 ms | 554.186 ms | 1.549x |

Aggregates:

| Measure | Unrestricted | Saturated |
|---|---:|---:|
| Sum of per-image medians | 3.918 s | 6.080 s |
| Mean per image | 391.780 ms | 608.049 ms |
| Median across images | 386.220 ms | 578.516 ms |

The saturated construction adds 216.269 ms per image on average. Its aggregate
ratio is 1.552x. Relative to the historical pre-cleanup run, the observed mean
time decreased by 9.24% for the unrestricted mode and 7.91% for the saturated
mode. These differences are empirical observations from one rerun and include
normal machine variability; they are not presented as isolated causal estimates
of instrumentation overhead.

Both modes reconstructed every input exactly during warm-up. Node counts and
saturation rejections matched the historical run for all ten images. Complete
measurements are available in
`benchmarks/results/min-max-residual-icdar-480p-first10-2026-08-04-post-cleanup.csv`.
