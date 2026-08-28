# Scientific benchmark builds

mmcfilters selects caller-contract validation at compile
time without changing the public C++ API. Normal builds use `CHECKED`.
Controlled scientific benchmarks may use `UNCHECKED` when the experiment has
already established every input domain.

This guide describes how to configure both modes and the validation boundaries
that keep their scientific workloads comparable. The benchmark material is
public for reproducibility, but it is not part of the public API contract or a
user-facing performance guarantee.

## Choosing a contract mode

| Mode | Intended use | Caller-contract diagnostics |
| --- | --- | --- |
| `CHECKED` | Normal applications, development, and experiments that include public API validation | Enabled |
| `UNCHECKED` | Controlled experiments whose inputs and domains are established before measurement | Disabled |

`UNCHECKED` removes only caller-contract diagnostics. Invalid ids, shapes,
pointers, or domains then violate preconditions and may cause undefined
behavior. Scientific capability checks and structural invariants remain
enforced because disabling them could change the defined scientific result or
permit an incoherent data structure.

All algorithms in one comparison must use the same contract mode, compiler
options, and optimization level.

## Configuring independent builds

Use separate build directories so that one experiment cannot mix objects
compiled with different contract modes:

```bash
cmake -S . -B build-checked \
  -DMMCFILTERS_BUILD_PYTHON=OFF \
  -DMMCFILTERS_BUILD_BENCHMARKS=ON \
  -DMMCFILTERS_CONTRACT_MODE=CHECKED

cmake -S . -B build-unchecked \
  -DMMCFILTERS_BUILD_PYTHON=OFF \
  -DMMCFILTERS_BUILD_BENCHMARKS=ON \
  -DMMCFILTERS_CONTRACT_MODE=UNCHECKED

cmake --build build-checked \
  --target mmcfilters_scientific_pipeline_benchmark
cmake --build build-unchecked \
  --target mmcfilters_scientific_pipeline_benchmark
```

The complete API benchmark suite and its reproducible execution protocol are
documented in the
[Scientific API benchmark suite](https://github.com/wonderalexandre/mmcfilters/blob/main/benchmarks/README.md).
Named, checksum-pinned workloads are provided in the
[`scientific-api.ini` manifest](https://github.com/wonderalexandre/mmcfilters/blob/main/benchmarks/workloads/scientific-api.ini).

## Validation boundaries

A public entry point validates its complete caller-owned input once and then
enters a validation-free execution core. It must not call another public
operation from a per-node or per-pixel loop.

Internal kernels therefore use committed tree, grid, and image access after one
public boundary. They do not repeat per-node or per-pixel validation even in a
`CHECKED` build.

### Classification

Every condition belongs to one of the following categories:

| Category | Examples | Policy |
| --- | --- | --- |
| Caller contract | valid node or pixel id, non-null pointer, buffer shape, finite threshold, positive window radius, supported padding, valid call sequence | `MMCFILTERS_CONTRACT_REQUIRE` or `MMCFILTERS_CONTRACT_CHECKED_ONLY` |
| Scientific capability | globally monotone altitude order required by MSER/extinction, adjacency information required by a projection | Always enforced when the algorithm cannot produce the defined scientific quantity without it |
| Structural invariant | coherent parent/child topology, valid smallest-node support assignment, monotone owned altitude after mutation, projection conservation, representable output | Always enforced |
| Established kernel precondition | live-node domain, valid parent/child links, dense buffer indexing, compatible image domain already established by the caller | No runtime validation inside the kernel |

### Validation-free execution cores

Validation-free cores live in `detail::kernel` and keep short operation names
such as `computeArea`, `filterDirect`, or `traversePostOrder`. Public functions
retain their existing names and validate before entering a kernel. Code named
`Impl` does not implicitly receive this guarantee.

Kernels use established spans and the narrow
`CommittedTreeAccess`/`CommittedImageAccess` facades for primitive tree and
image access. Helpers outside `detail::kernel` may also use those facades after
their own dominating boundary has established the same domain.

For example, an MSER computation validates the altitude buffer, window radius, and
hierarchy capability once. Its ancestor search then reads altitude slots and
parent links directly for every live node. Likewise, a reconstruction filter
validates its node-preservation mask and output image once before traversing committed child
and proper-part ranges.

### Redundant checks

A defensive check is removed rather than wrapped when its truth follows from a
dominating boundary or from an always-maintained structural invariant. Typical
examples are rechecking that every node returned by the live-node iterator is
alive, revalidating an altitude span inside every ancestor search, and checking
the same topology mutation version in both an adapter and the implementation it
immediately calls.

The validation-boundary audit rejects contract checks, defensive public
accessors, and explicit exceptions inside every `detail::kernel` namespace.
Tests run in both compile-time modes: invalid-input expectations are enabled
only in `CHECKED`, while scientific result checks and invariant tests run in
both modes.

## Scientific comparison rules

A benchmark must report the selected contract mode and verify an identical
scientific checksum between `CHECKED` and `UNCHECKED` builds before timing
results are interpreted. The comparison must preserve the same mathematical
workload, compiler options, optimization level, and measurement scope for every
algorithm.

The full API benchmark runner alternates process order and checks deterministic
cross-mode results. Its manifest separates end-to-end scenarios from operations
over established inputs, allowing validation costs to be reported without
changing the scientific operation being measured.

## Distance-transform scalar studies

The sensitivity study aligns all 29 attributes in `DIST_TRANSF` with their
`DIST_TRANSF_EXACT` counterparts on the same live nodes:

```bash
PYTHONPATH=build/python python3 benchmarks/distance_transform_sensitivity.py \
  /path/to/icdar benchmarks/results/distance-transform-sensitivity \
  --count 10 \
  --trees max,min,tos,residual_unrestricted,residual_saturated \
  --radius 1.5 --infinity-pixel 0
```

The timing study measures established-tree attribute calls at 480p, 720p, and
1080p. It compares both complete groups and the isolated maximum attributes,
with three repetitions per image/tree/mode:

```bash
PYTHONPATH=build/python python3 benchmarks/distance_transform_timing.py \
  /path/to/icdar benchmarks/results/distance-transform-timing \
  --count 10 \
  --trees max,min,tos,residual_unrestricted,residual_saturated \
  --resolutions 480p,720p,1080p --repetitions 3 \
  --expected-contract-mode CHECKED
```

The ToS case uses the default canonical minimum-4/maximum-8 complementary
immersion with `TopographicDomainExtension.NONE`, `infinity_pixel=0`, and
`uint8` altitudes. Resizing and tree construction are outside the attribute
timer; construction time and node counts are recorded separately. Input hashes,
result hashes, native-module metadata, observation order, and raw per-call
times are retained with the generated analysis.

Before interpreting timings, run the checksum-pinned
`benchmarks/workloads/distance-transform-contract.ini` workload in optimized
`CHECKED` and `UNCHECKED` builds. The grouped and sequential scalar checksums
must agree within each build, and the scientific checksums must match across
contract modes.

These scripts measure elapsed execution time; they do not perform hotspot
profiling. Runtime attribution must use an optimized build and external native
profiling tools, keeping profiler output outside public source and benchmark
artifacts. Production kernels must not gain timers, counters, or diagnostic
branches solely for profiling.
