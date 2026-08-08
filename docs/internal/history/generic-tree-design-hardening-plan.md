# Generic Morphological-Tree Design Hardening Plan

## Goal

Complete the separation between generic partial-partition topology and concrete
image-tree producers without regressing the runtime or memory behavior of
existing construction, adjustment, attribute, and reconstruction paths.

The central model remains unchanged:

- internal nodes and proper parts use independent dense id domains;
- each proper part has exactly one direct owner;
- a live node may own no direct proper part, but its full subtree support must
  be non-empty;
- `MorphologicalTree` owns topology;
- `WeightedMorphologicalTree<T>` owns topology plus typed node altitudes;
- geometry and adjacency remain optional capabilities.

This plan does not add a 3D domain and does not add Tree-of-Shapes sidecar
metadata.

## Performance Contract

Every phase must preserve:

- the asymptotic complexity of public methods;
- the number of dominant full-domain passes in hot paths;
- the max-tree/min-tree union-find ordering and node-id assignment;
- the ToS builder algorithm and its output buffers;
- allocation behavior on per-pixel/per-node inner loops;
- the canonical public C++ and Python call shapes.

A phase is blocked when a repeatable benchmark regression remains after warmup.
Timing comparisons use medians from at least five process runs. A difference
inside normal machine noise is not treated as a semantic performance change,
but no new linear pass, per-neighbour allocation, or global validation may be
hidden behind that tolerance.

Correctness tests do not contain timing thresholds. Benchmarks report timings
and structural counters separately.

## Phase 0 — Baseline And Permanent Performance Coverage

Deliverables:

- retain the generic native construction/commit benchmark;
- add equivalent deep and branching Higra-import cases;
- retain the broader construction/attribute/CASF/oracle benchmark;
- record the exact compiler mode and representative baseline medians.

Gate:

- full configured tests pass before structural changes;
- Release benchmarks produce stable checksums and output schemas.

Initial Release baseline on 2026-07-23 (`c++ -O3 -DNDEBUG`, Apple Clang,
100,000 internal nodes, 10 empty checked commits):

| case | warm construction median | warm commit batch median |
| --- | ---: | ---: |
| native deep | 1.23 ms | 8.33 ms |
| native branching | 1.07 ms | 8.02 ms |
| Higra deep | 1.55 ms | 11.38 ms |
| Higra branching | 1.57 ms | 11.31 ms |

The broader baseline also covers max/min construction, ToS, `GRAY_HEIGHT`,
`MAX_DIST`, checked commit, and CASF.

## Phase 1 — Extract Concrete Import Policy Without Behavioral Change

Deliverables:

- introduce a small internal `HigraImportLayout`;
- move leaf-first id validation and id conversion out of
  `MorphologicalTree`;
- keep public factories as delegating producer/import facades;
- preserve the exact number of loops, allocations, link operations, ids, and
  diagnostics.

Gate:

- Higra malformed-input, round-trip, projection, and edit-invalidation tests
  pass;
- Higra import benchmark does not regress.

## Phase 2 — Proof-Carrying Native Materialization

Deliverables:

- replace the `Checked`/`TrustedProducer` policy enum with a move-only validated
  native-hierarchy representation;
- centralize raw-buffer validation so public imports and internal producers use
  the same invariant definition;
- consume owning buffers by move where producers already own them;
- retain debug post-materialization validation as an oracle;
- ensure each production path performs at most one complete structural
  validation.

For producers whose invariants are established during construction, validation
facts must be accumulated in their existing passes. This phase must not add a
second full-domain pass to ToS.

Gate:

- invalid cycles, detached nodes, invalid owners, support-empty nodes,
  incompatible grid/adjacency, and invalid altitude order are rejected;
- ToS buffers remain identical;
- construction benchmarks do not regress.

## Phase 3 — Generic Incremental Edit Proof

Deliverables:

- replace concrete adjuster friendship and `commitTrusted()` with a generic,
  move-only edit proof;
- maintain a delta ledger in existing edit primitives for the safe generic
  path;
- validate touched nodes, changed parent paths, changed proper-part ownership,
  and changed altitude edges;
- use full checked commit as fallback for operations not yet supported by the
  incremental proof;
- let internal algorithms that already establish every invariant in their
  existing passes issue the same proof through one generic
  established-by-construction boundary;
- keep debug full validation as an oracle after an accepted delta proof.

The safe Release path must remain proportional to the mutation delta. An
established-by-construction hot path must not duplicate the algorithm's existing
local checks or scan every node or proper part after each adjustment step.

Gate:

- both dual min/max adjusters produce exactly the same results as full checked
  commit;
- instrumentation reports zero full scans in covered hot commits;
- CASF and paired-oracle benchmarks do not regress.

## Phase 4 — Public Mutation And RAII Safety

Deliverables:

- move unchecked altitude setters behind the internal mutation boundary;
- make the default public edit recoverable when abandoned;
- preserve a snapshot-free internal editor for hot algorithms;
- prefer an undo log proportional to the mutation delta over a complete tree
  snapshot;
- keep the public edit surface limited to checked publication paths.

Gate:

- exceptions and abandoned public editors restore a usable committed owner;
- public methods cannot publish invalid altitude order or topology;
- ordinary single-step edits do not acquire a full-tree copy.

## Phase 5 — Complete Producer Separation

Deliverables:

- add an owning native hierarchy representation complementary to
  `NativeHierarchyView<T>`;
- move max-tree/min-tree image construction into a component-tree producer;
- route the producer through the native materialization boundary;
- move Higra import/export projection policy toward an interoperability module;
- keep current producer/import factories as delegating facades.

Gate:

- union-find ordering, node ids, owners, altitudes, and exported Higra buffers
  are byte-identical;
- max/min construction and Higra import benchmarks do not regress.

## Phase 6 — Reentrant 2D Adjacency

Deliverables:

- make stored grid adjacency immutable;
- move traversal cursor state into iterator/range objects;
- remove mutable adjacency access from committed trees;
- name the capability explicitly as regular-grid 2D adjacency;
- avoid type erasure, virtual dispatch, and allocation in neighbour loops.

Gate:

- adjacency, contour, BitQuad, `MAX_DIST`, saliency, and CASF tests pass;
- neighbourhood iteration microbenchmarks and end-to-end consumers do not
  regress.

## Phase 7 — Exact ToS Altitudes Without Sidecar Metadata

Deliverables:

- consolidate interpolation/connectivity, optional exterior padding, and
  \f$p_\infty\f$ in `TreeOfShapesProducerOptions`;
- expose `TreeOfShapesProducer` as the producer name;
- add an opt-in typed ToS factory that preserves doubled construction levels;
- retain the quantized `uint8_t` convenience factory;
- keep proper parts and `GridDomain2D` on the original image domain;
- use exact altitude order for directional-adjacency consumers;
- do not add ToS-specific fields to `MorphologicalTree`.

Gate:

- original-domain reconstruction remains exact;
- padded and unpadded construction both publish only the original proper-part
  domain;
- quantized output is unchanged;
- exact output removes quantization-induced equal-arc ambiguity;
- ToS construction does not gain another image-domain pass.

## Phase 8 — API And Installation Cleanup

Deliverables:

- migrate internal use from image/pixel names to grid-domain/proper-part names;
- remove deprecated compatibility spellings after migrating internal callers;
- keep Python owner-oriented by removing the dead topology-only registration
  and overloads;
- keep current Python factory/import return types;
- define an explicit public-header manifest and stop installing internal
  implementation headers;
- make attribute capability requirements part of the descriptor metadata
  table.

Gate:

- installed C++ consumer and complete Python suite pass;
- public examples compile against the canonical API;
- attribute kernels and scheduling benchmarks do not regress.

## Validation Matrix

After every phase:

```bash
cmake --build /tmp/mmcfilters-tos-python-build -j4
ctest --test-dir /tmp/mmcfilters-tos-python-build --output-on-failure
python -m compileall -q python/mmcfilters unit-tests/python
git diff --check
```

Performance gates:

```bash
c++ -std=c++20 -O3 -DNDEBUG -I. \
  benchmarks/tree_validation_benchmark.cpp \
  -o /tmp/mmcfilters_tree_validation_benchmark

/tmp/mmcfilters_tree_validation_benchmark 100000 10
```

The broader benchmark additionally covers max/min, ToS, attributes, checked
commit, and CASF.

## Implementation Status

- Phase 0: completed; Release baselines recorded and permanent native/Higra
  benchmark cases added.
- Phase 1: completed; `detail::HigraImportLayout` now owns leaf-first id-domain
  validation and constant-offset conversion. The existing import loops,
  allocations, public API, ids, and diagnostics are preserved.
- Phase 1 validation:
  - complete configured C++/Python suite passes (`53/53`);
  - Python bytecode compilation and `git diff --check` pass;
  - seven interleaved before/after Release runs measured Higra deep
    construction at `1.579 ms` versus `1.595 ms`, and Higra branching at
    `1.639 ms` versus `1.665 ms`; the approximately 1–2% differences are inside
    the observed process noise and no pass, allocation, or inner-loop
    operation was added.
- Phase 2: completed.
  - `detail::ValidatedNativeHierarchy<T>` is the move-only owning
    materialization representation;
  - `detail::NativeTopologyProof` replaces the
    `Checked`/`TrustedProducer` policy;
  - public `NativeHierarchyView<T>` imports receive one complete centralized
    validation before their buffers are copied once and then transferred;
  - ToS records connected, supported nodes and proper-part owners in its
    existing projection loops;
  - producer-owned parent, owner, and altitude vectors are transferred into
    the weighted owner instead of being copied at the factory boundary;
  - assertion-enabled materialization retains the full
    `validateConnectedRootedTree()` oracle.
- Phase 2 validation:
  - complete configured C++/Python suite passes (`53/53`);
  - all malformed native-input checks remain covered, including cycles,
    invalid owners, empty subtree support, grid/adjacency mismatch, altitude
    shape, and altitude order;
  - seven interleaved Release runs measured ToS construction at a median
    `3.944 ms` before and `3.987 ms` after (about `1.1%`, inside process
    noise);
  - the native deep-import median improved from `1.272 ms` to `0.532 ms`, and
    the branching-import median from `1.104 ms` to `0.406 ms`, because
    validated owning buffers are no longer recopied during materialization;
  - checksums and public return types are unchanged.
- Phase 3: completed.
  - `TreeEditor::IncrementalProof` is move-only and bound to the editor, tree,
    and exact mutation version;
  - the safe generic editor records unique touched nodes, detached-node and
    unsupported-leaf balances, ownership endpoints, and changed parent paths
    in storage proportional to the edit delta;
  - supported edits validate only that ledger, while unsupported primitives
    fall back to the complete validator;
  - `WeightedTreeEditor<T>` additionally validates monotone altitude on touched
    arcs and can use strict altitude order as generic acyclicity evidence;
  - both dual min/max adjusters now publish through the generic proof protocol;
    concrete adjuster friendship and `commitTrusted()` were removed;
  - the adjusters use a generic internal established-by-construction entrypoint
    because forcing duplicate ledger work into that already-validated loop
    produced a repeatable Release regression during implementation;
  - assertion-enabled builds retain complete topology and altitude oracles;
  - `TreeEditValidationStatistics` makes complete versus incremental
    publication strategies observable without timing assertions.
- Phase 3 validation:
  - complete configured C++/Python suite passes (`53/53`);
  - tests reject copied proofs at compile time, stale proofs, changed-parent
    cycles, and invalid touched altitude arcs;
  - CASF tests report incremental publications and zero complete publications
    on the covered hot path;
  - seven interleaved Release runs measured the CASF filter median at
    `13.171 ms` before and `13.354 ms` after (about `1.4%`, inside process
    noise), while total CASF time changed from `16.379 ms` to `16.322 ms`;
  - reverse-order runs did not reproduce an oracle regression: the adjusted
    paired-tree median was `9.931 ms` after versus `10.260 ms` before;
    checksums remained identical.
- Phase 4: completed.
  - every public `TreeEditor` now provides explicit and destructor-driven
    rollback;
  - the journal is allocated lazily before the first mutation and captures each
    affected node/proper-part/altitude slot only on its first write;
  - dense node-vector growth and free-list pushes/pops are restored without a
    complete-tree copy;
  - abandoned editors and exception unwinding restore topology, ownership,
    sibling/proper-part order, root, live/free slots, mutation versions, and
    weighted altitudes;
  - `edit()` uses the delta-journal behavior instead of cloning the complete
    owner;
  - the generic internal established-by-construction editor remains
    journal-free;
  - public `setAltitudeUnchecked()` and
    `setAltitudeBufferUnchecked()` were removed; staged altitude changes must
    cross the weighted editor validation/proof boundary.
- Phase 4 validation:
  - complete configured C++/Python suite passes (`53/53`);
  - rollback tests cover append/reuse/release, reparenting, child and
    proper-part splices, root promotion, subtree pruning, parent merge,
    explicit rollback, destruction, exception unwinding, and invalid weighted
    altitude;
  - ten interleaved/reverse-order Release pairs measured the CASF filter median
    at `13.401 ms` before and `13.428 ms` after (about `0.2%`);
  - empty checked-commit median changed from `0.665 ms` to `0.656 ms`;
  - total CASF time remained within observed process variation
    (`16.420 ms` before, `16.704 ms` after), and all checksums were identical.
- Phase 5: completed.
  - `detail::ComponentTreeProducer<T>` now owns max-tree/min-tree image
    construction and emits proven native parent, owner, and altitude buffers;
  - the historical union-find pixel ordering, dense node ids, child order,
    direct ownership, altitude representatives, and construction mutation
    version are preserved;
  - max/min factories delegate to the common owning native materialization
    boundary, and the image-specific topology and altitude-inference
    constructors were removed from the generic owners;
  - `detail::adaptHigraHierarchy<T>` now owns leaf-first id validation,
    node/owner projection, typed altitude projection, and structural proof;
  - Higra imports also cross the native materialization boundary, while
    `MorphologicalTree` retains only a generic optional external node-id offset
    used by preserved Higra-domain projection;
  - export and attribute projection continue to share
    `detail::ExportedHigraLayout`.
- Phase 5 validation:
  - complete configured C++/Python suite passes (`53/53`);
  - exact fixture tests pin max/min parent, proper-part-owner, and altitude
    buffers, existing tests pin traversal order and Higra round trips, and a
    signed-zero plateau test verifies the historical first-row-major altitude
    representative;
  - five forward and five reverse interleaved Release pairs found no max/min
    construction regression; in both run orders the phase-5 medians were
    faster than the phase-4 binaries, while CASF and unrelated ToS cases
    remained within observed process variation;
  - seven warm Release runs measured Higra deep and branching construction at
    medians of approximately `0.685 ms` and `0.812 ms` for 100,000 internal
    nodes;
  - broad benchmark checksums are byte-identical to Phase 4.
- Phase 6: completed.
  - `RegularGridAdjacency2D` is the explicit adjacency capability;
  - traversal state lives in small forward-range/iterator values, so nested and
    interleaved traversals of the same const relation are reentrant;
  - forward-only offset indices are precomputed once and hot neighbor loops
    retain static dispatch, contiguous stencil storage, and zero allocations;
  - committed trees expose only const uniform and directional adjacency
    pointers through the canonical `*GridAdjacency2D` accessors;
  - component-tree, saliency, shape-space, and CASF
    consumers now borrow adjacency as const, and saliency no longer creates
    defensive adjacency copies;
  - the same immutable concrete relation supports Euclidean-radius disks,
    centered rectangles, digital lines, and centrally symmetric structuring
    elements; max/min factories accept an explicit relation;
  - adjacency-inducing structuring elements require one origin and central
    symmetry, preserving the undirected forward-edge contract;
  - BitQuad capability validation rejects noncanonical stencils instead of
    silently treating every non-4-connected relation as 8-connected.
- Phase 6 validation:
  - complete configured C++/Python suite passes (`54/54`), including the
    installed-consumer check;
  - new tests cover const access, C++20 forward-range conformance, exact legacy
    order, bounds, nested/interleaved traversal, rectangles, horizontal,
    vertical, and oriented lines, arbitrary symmetric elements, invalid
    elements, explicit max-tree construction, and exact equivalence between a
    custom cross and radius-based 4-connectivity;
  - nine warm Release runs preserved the neighborhood checksum and improved
    median full traversal from `1.818 ms` to `1.782 ms`; forward-only traversal
    improved from `1.184 ms` to `1.011 ms`;
  - the broad benchmark checksum is unchanged; max/min construction improved
    in the sampled medians and ToS, attributes, checked commit, CASF, and
    the paired-tree oracle remained within observed process variation.
- Phase 7: completed.
  - `TreeOfShapesProducerOptions` owns the three built-in interpolation choices,
    exterior/no-padding policy, and \f$p_\infty\f$ coordinates in the selected
    transient immersion domain;
  - `TreeOfShapesProducer` is the producer name;
  - both padding choices publish exactly one proper part per source pixel and
    the original `GridDomain2D`; no producer option or ToS sidecar was added to
    `MorphologicalTree`;
  - `createTreeOfShapesExact()` returns
    `WeightedMorphologicalTree<ToSGrayLevel>` in doubled gray units, while the
    existing factory retains its `uint8_t` return type and quantization;
  - altitude encoding is a compile-time producer policy, so no `std::function`,
    type erasure, or virtual dispatch was introduced in construction loops;
  - C++ and Python cover unpadded construction for all three interpolation
    choices, original-domain reconstruction, seed bounds, and retained
    half-level structural nodes;
  - the complete configured suite passes (`54/54`);
  - seven clean interleaved Release runs preserved the broad benchmark checksum
    (`3731633`) and measured median default ToS construction at `4.200 ms`
    before versus `4.248 ms` after (about `1.14%`, within process noise).
- Phase 8: completed.
  - generic topology internals and public accessors now use proper-part and
    regular-grid-domain terminology; the adjacency type is
    `RegularGridAdjacency2D`;
  - canonical adjacency traversal is expressed in grid indices
    (`getAdjacentIndices()`, `getNeighborIndices()`, and
    `getForwardNeighborIndices()`);
  - tree, edit, adjacency, and Tree-of-Shapes APIs use one canonical spelling;
  - Python factory/import return types remain `WeightedMorphologicalTree`.
    No separate topology view or unconstructible topology-only registration is
    exposed;
  - `cmake/mmcfiltersPublicHeaders.cmake` explicitly declares the stable include
    entry points and their exact transitive closure. Template-support `detail`
    headers remain installable without becoming independent API contracts,
    while unused helpers are no longer installed;
  - each `AttributeMetadata` row now owns its complete capability requirements;
    the parallel classification switch was removed, ordinal alignment is
    checked at compile time, and Python exposes the complete contract including
    `canonical4Or8Adjacency`.
- Phase 8 validation:
  - the complete configured C++/Python suite passes (`54/54`), including the
    installed-consumer test with canonical APIs, exact/unpadded ToS
    construction, and assertions that removed compatibility and internal
    oracle headers are absent;
  - every configured example target compiles; Python package/test sources pass
    `compileall`, and the installed manifest equals the complete transitive
    closure of the declared public entry points (no missing or extra support
    header);
  - seven clean interleaved Release runs preserved the broad checksum
    (`3731633`). Max/min construction changed by `+0.23%`/`+0.50%`, default ToS
    by `-0.60%`, `MAX_DIST` by `+0.23%`, and CASF
    filter/total by `-0.14%`/`+0.09%`;
  - the two sub-millisecond microcases moved by `0.013 ms` (gray height) and
    `0.030 ms` (checked commit) and reversed or narrowed across repeated
    seven-run sets, consistent with measurement noise rather than a changed
    hot-path implementation.
  - after removing the compatibility surface and renaming adjacency traversal
    to `Indices`, the Release benchmark's complete `__text` section remained
    byte-identical to the Phase 8 baseline and all 14 interleaved executions
    preserved checksum `3731633`; observed timing differences therefore reflect
    measurement noise rather than generated-code changes.
