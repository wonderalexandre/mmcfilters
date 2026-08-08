# Contour Trace API Implementation Plan

Status: initial C++/Python implementation completed on `contour-trace-api`;
this document remains the design and follow-up plan for the definitive
geometric contour API.
Created: 2026-05-31.

## Motivation

The current `ContoursComputedIncrementally` API materializes one flat set of
support pixel ids per tree node. A pixel enters the contour when at least one
4-neighbour side is exposed. This is compact and useful for scalar descriptors,
but it loses side-level geometry:

- one contour pixel may expose several sides;
- exposed sides may touch different complement components;
- external and internal boundaries cannot be separated reliably from the pixel
  set alone;
- loop geometry and orientation are not represented.

The contour trace API should become the definitive geometric contour API for
new contour work in this repository. It is not a transitional wrapper around
the pixel-contour API. Its primitive should be a grid edge, not a contour
pixel.

## Goals

- Add a definitive contour-trace API without breaking the current
  pixel-contour API.
- Represent boundaries as oriented grid edges attached to support pixels.
- Trace materialized edges into ordered loops.
- Classify loops as external or internal from orientation under a documented
  convention.
- Preserve the incremental tree workflow: compact local deltas first, lazy
  materialization per node, optional full prefetch.
- Keep total contour-edge count compatible with existing perimeter attributes.
- Expose a usable Python API after the C++ contract is stable.
- Make the full all-node extraction/materialization/tracing target
  `O(N + P + sum_v B_v)`, where `B_v` is the boundary-edge count of node `v`.

## Non-Goals

- Do not remove `ContoursComputedIncrementally` in the first implementation;
  keep it as a compatibility and scalar-descriptor path while trace becomes the
  canonical geometric contour path.
- Do not change existing `CONTOUR_*` attribute semantics.
- Do not add a per-node complement flood fill as the primary method.
- Do not infer loop separation from `getContour(node)` pixel lists.
- Do not optimize for GPU execution in this iteration.

## Mathematical Model

Let the image domain be a rectangular grid of pixels. For a tree node `v`, let
`S(v)` be the set of proper parts in the connected component represented by the
subtree rooted at `v`.

A boundary edge is a pair `(p, s)` where:

- `p` is a pixel in `S(v)`;
- `s` is one of `North`, `West`, `East`, or `South`;
- the neighbour of `p` across side `s` is either outside the image domain or
  not in `S(v)`.

The boundary edge is oriented so that the support pixel remains on a fixed side
of the directed edge. The exact sign convention must be fixed in tests, because
image coordinates use rows growing downward. Once fixed, external and internal
loops have opposite signed area.

## Public C++ API Sketch

The API should live under `mmcfilters/contours/` and mirror the current
incremental-contour style.

```cpp
enum class ContourTraceSide {
    North,
    West,
    East,
    South
};

enum class ContourLoopKind {
    External,
    Internal
};

struct ContourTraceEdge {
    int pixel = -1;
    ContourTraceSide side = ContourTraceSide::North;
};

struct ContourLoopInfo {
    ContourLoopKind kind = ContourLoopKind::External;
    uint32_t edgeOffset = 0;
    uint32_t edgeCount = 0;
    int signedArea2 = 0;
};

class ContourTraceComputation {
public:
    struct IncrementalContourTraces {
        LoopRange getLoops(NodeId node) const;
        EdgeRange getLoopEdges(const ContourLoopInfo& loop) const;
        EdgeRange getEdges(NodeId node) const;
        void materializeAll() const;
        bool isMaterialized() const;
        bool isNodeTraced(NodeId node) const;
    };

    static IncrementalContourTraces extract(const MorphologicalTree& tree);

    template<AltitudeValue T>
    static IncrementalContourTraces extract(const WeightedTreeView<T>& tree);
};
```

The exact range classes can follow the current `ContourRange` and
`ContoursByNodeRange` pattern. If nested range lifetimes become awkward for
pybind11, the C++ API can remain range-based while Python returns lists of loop
objects for one node.

## Internal Data Model

Use packed side ids in compact stores:

```text
packed_edge = 4 * pixel + side_index
side_index: North=0, West=1, East=2, South=3
```

Suggested files:

```text
mmcfilters/contours/ContourTraceComputation.hpp
mmcfilters/contours/detail/ContourTraceDeltaStore.hpp
pybinds/ContourTraceBindings.cpp
pybinds/ContourTraceBindings.cpp
unit-tests/contours/test_contour_traces_on_morphological_tree.cpp
```

The first implementation keeps geometry helpers inside
`ContourTraceComputation.hpp`. A separate `ContourTraceGeometry.hpp` can still
be introduced later if the tracing rules grow beyond the current header.

`ContourTraceDeltaStore` should be analogous to `ContourDeltaStore`, but its
values are packed boundary edges rather than pixels:

```text
addValues        = compact local edge additions
addSpans[node]   = offset/size into addValues
removeValues     = compact local edge removals
removeSpans[node]= offset/size into removeValues
```

The materialized cache should keep:

```text
cachedEdgeValues_[]
cachedEdgeOffset_[node]
cachedEdgeSize_[node]
cachedEdgeReady_[node]

cachedLoopInfos_[]
cachedLoopInfoOffset_[node]
cachedLoopInfoSize_[node]
cachedLoopReady_[node]
```

Keep edge materialization and loop tracing as separate readiness flags. Some
callers may only need unordered boundary edges or perimeter-compatible counts.
After every live edge cache is ready, discard compact local edge deltas and the
edge-mark scratch buffer because later `getEdges`/`getLoops` calls can read the
materialized slices directly.
The vertex-to-outgoing-edge tracing scratch should be allocated only when loop
tracing is requested, so extraction and edge-only workloads do not pay the
image-vertex table cost.
For isolated node-local `getLoops(node)` calls, the implementation should use a
sparse outgoing-vertex table when that table is smaller than the dense
image-vertex table. Full `materializeAll()` should keep the dense table because
all-node tracing repeatedly looks up grid vertices and the dense array is
faster. A long sequence of node-local traces should also switch back to the
dense table after a small number of sparse queries, preserving interactive
memory savings without making all-node random access pay hash-table overhead.
When tracing a node, reuse its `cachedEdgeValues_` segment for ordered loop
edges whenever the parent edge cache is already materialized. If the node is
traced before its parent can safely consume the original edge order, append the
ordered loop edges to `cachedEdgeValues_` as a fallback slice. This avoids a
persistent second global loop-edge cache while preserving lazy access-order
independence.
During full `materializeAll()`, reserve and fill global loop metadata directly.
During node-local `getLoops(node)`, keep a local loop-info buffer before
committing to the global cache so random query order does not inflate retained
capacity.
The loop follower should take a direct successor fast path for vertices with a
single outgoing boundary edge; only ambiguous vertices need the full
right/straight/left/back tie-breaking scan.
Successor selection must return only unvisited edges. With that invariant, the
inner loop can mark and consume the selected edge directly without repeating
the visited test on every traced edge.

## Edge Delta Extraction

The current contour algorithm tracks whether a pixel has at least one exposed
side. The trace API should instead track each side independently.

For every direct proper part `p` owned by the currently processed node `ownerP`,
scan the four side neighbours:

1. If the neighbour is outside the image domain:
   - add edge `(p, side)` at `ownerP`;
   - there is no removal event.
2. If the neighbour `q` is inside the image domain:
   - compute `entry = properPartEntryNode(tree, p, q)`;
   - if `entry == InvalidNode`, ignore the event;
   - if `entry == ownerP`, `q` is already in the support at `ownerP`, so no
     boundary edge is active there;
   - otherwise add edge `(p, side)` at `ownerP` and remove it at `entry`.

This models the interval of tree nodes where `p` is visible but the neighbour
across that side is not yet visible. It removes the need for the pixel-level
`ncount` bookkeeping used by the old API.

The extraction traversal should still be post-order and should still validate:

- non-empty image domain;
- live root;
- stable `WeightedTreeView` for weighted overloads;
- live dense internal node ids for public reads.

## Lazy Edge Materialization

Materialization follows the existing contour cache pattern:

1. Walk the requested subtree in post-order.
2. Reuse already materialized child edge caches.
3. Accumulate child edges into a temporary vector with a generation-marked edge
   set of size `4 * P`.
4. Add local edge additions.
5. Apply local edge removals.
6. Commit the compact edge slice for the node.

The edge mark buffer is scratch storage. It should use the same generation
strategy as the pixel-contour cache, but indexed by packed edge id.

## Loop Tracing

Tracing should operate on the materialized edge set of one node.

For each packed edge:

1. Convert `(pixel, side)` to two integer grid vertices.
2. Orient the directed edge according to the support-side convention.
3. Build a temporary vertex-to-outgoing-edge adjacency.
4. Follow unvisited directed edges until the start edge is reached.
5. Emit one `ContourLoopInfo` plus its ordered edge slice.

Ambiguous vertices are possible when boundaries touch at a grid corner. The
implementation must choose a deterministic half-edge successor rule. The rule
should be documented and tested. A good default is a wall-following rule that
keeps the support on the same side of the directed edge and uses a fixed turn
priority at vertices. The exact priority should be validated on 2x2 diagonal
fixtures before being treated as stable.

After each loop is traced, compute its doubled signed area from the grid
vertices. Under the chosen orientation convention, one sign corresponds to
external loops and the opposite sign corresponds to internal loops. The sign
mapping should be calibrated by a one-pixel fixture:

- one isolated pixel must produce one external loop;
- a 3x3 ring must produce one external loop and one internal loop.

## Python API Sketch

Initial Python exposure can prioritize clarity over zero-copy nested ranges:

```python
traces = mmcfilters.ContourTraceComputation.extraction(tree)

for loop in traces.getLoops(node_id):
    print(loop.kind, loop.signedArea2)
    for edge in loop.edges:
        print(edge.pixel, edge.side)
```

The Python objects may materialize small lists for one node while the C++ object
keeps compact caches internally. This avoids fragile iterator lifetime issues
for nested loop/edge ranges.

## Compatibility Checks

For every tested node:

- `sum(loop.edgeCount for loop in loops) == CONTOUR_PERIMETER`;
- the set of `edge.pixel` values equals the current pixel-contour set, after
  projecting edges to pixels and removing duplicates;
- the current `CONTOUR_SIDE_NORTH/WEST/EAST/SOUTH` attributes equal directional
  counts from traced edges;
- `BITQUADS_NUMBER_HOLES` should match the number of internal loops for simple
  connected supports under the matching connectivity convention.

The last check should be treated carefully for Tree of Shapes or ambiguous
digital topology cases, because bitquad scalar projection has its own
connectivity convention.

## Complexity

Let:

- `P` be the number of image pixels;
- `N` be the number of internal node slots;
- `N_live` be the number of live nodes;
- `B(S)` be the total number of boundary edges committed while materializing a
  subtree `S`;
- `M(S)` be the number of missing live nodes visited during materialization;
- `B_v` be the number of boundary edges of node `v`;
- `L_v` be the number of loops of node `v`.

Extraction with valid tree-query caches:

```text
O(N + P)
```

Cold-cache extraction:

```text
O(T_tree_cache + N + P)
```

where `T_tree_cache` is the same ancestry/LCA preprocessing cost described in
`docs/contours.md`.

First edge materialization for a subtree:

```text
O(M(S) + B(S) + D(S))
```

where `D(S)` is the number of compact local edge additions/removals read for
the newly materialized part of `S`. Because the side stencil is fixed, the
global delta count is bounded by a constant factor of `P`.

Tracing one materialized node:

```text
O(B_v + L_v)
```

Materializing and tracing all nodes should meet the definitive API target:

```text
O(N + P + sum_v B_v)
```

This absorbs loop bookkeeping because every non-empty loop contains at least
one boundary edge, so `sum_v L_v <= sum_v B_v`. Any design that requires
reconstructing or flood-filling the full image domain independently for each
node would fall back to `O(NP)` and should be rejected for this API.

Memory:

- local edge deltas: `O(P)` with a larger constant than pixel contours;
- edge scratch mark buffer: `O(4P)`;
- local edge deltas and the edge scratch mark buffer can be released after all
  edge caches are materialized;
- materialized edge cache: `O(sum materialized B_v)`;
- loop cache: `O(sum traced L_v)` when loop edges reuse materialized edge
  segments, with possible fallback duplicate slices for nodes traced before
  their parents;
- one-node tracing scratch: `O(B_v)`.

## Implementation Milestones

1. Add edge packing and geometry helpers.
   - Implement side enum, packing/unpacking, neighbour lookup, vertex mapping,
     and orientation conversion.
   - Add small unit tests independent of trees.

2. Add compact edge-delta extraction.
   - Implement `ContourTraceDeltaStore`.
   - Implement `ContourTraceComputation::extractTraceDeltas`.
   - Validate additions/removals against hand-computable masks.

3. Add lazy edge materialization.
   - Mirror the current `IncrementalContours` cache pattern.
   - Add stale-tree checks and live-node validation.
   - Compare unordered edge counts against `CONTOUR_PERIMETER`.

4. Add loop tracing.
   - Build directed edge adjacency from one materialized node.
   - Implement deterministic successor rule for ambiguous vertices.
   - Classify loops by signed area after calibrating sign convention.

5. Add public C++ API.
   - Expose `getEdges(node)`, `getLoops(node)`, loop-edge access, diagnostics,
     and `materializeAll`.
   - Keep existing contour APIs unchanged.

6. Add Python bindings.
   - Bind enums, edge objects, loop objects, and extraction factory.
   - Prefer list materialization per queried node if nested ranges are fragile.

7. Add docs and benchmark.
   - Document API and tracing convention.
   - Add a benchmark beside `examples/contour_benchmark.cpp` or extend it with
     trace-specific timings.

## Test Plan

C++ unit tests:

- one isolated pixel: four edges, one external loop;
- solid rectangle: one external loop, expected perimeter and signed area;
- 3x3 ring: one external loop, one internal loop;
- two-hole fixture: one external loop and two internal loops;
- shape touching image border: border edges handled as external;
- diagonal-touch ambiguity fixture: deterministic loop decomposition;
- max-tree and min-tree fixtures from existing contour tests;
- Tree of Shapes fixture from existing contour tests;
- mutation after extraction rejects stale contour-trace reads;
- `WeightedTreeView` stale topology is rejected.

Python tests:

- import smoke;
- extraction from topology and weighted tree;
- loop list contains stable edge objects with pixel and side;
- simple ring fixture exposes external/internal loop counts.

Cross-check tests:

- total traced edge count equals `CONTOUR_PERIMETER`;
- directional edge counts equal `CONTOUR_SIDE_*`;
- projected edge pixels equal existing `getContour(node)` pixels;
- simple internal loop count agrees with `BITQUADS_NUMBER_HOLES`.

## Risks and Open Decisions

- Orientation sign must be fixed under image row/column coordinates and locked
  with tests.
- Ambiguous vertex successor rules can change loop decomposition without
  changing perimeter. This needs explicit design and regression fixtures.
- Tree of Shapes projected supports may expose cases not seen in ordinary
  component trees; tests should include ToS early.
- Full materialization can use much more memory than pixel contours because it
  stores side-level geometry and ordered loop caches.
- Python nested ranges may be difficult to make ergonomic and safe; materialized
  per-node Python lists may be a better first API.

## Implementation Status

The definitive implementation has been completed with C++ and Python coverage:

1. edge packing and geometry helpers;
2. edge delta extraction;
3. lazy edge materialization;
4. perimeter and directional-count cross-checks.

Loop tracing, Python bindings, internal profiling, benchmark coverage, sparse
node-local adjacency, and scratch-release optimizations were added after the
edge-set cross-checks were in place. Further performance work should be driven
by profiler evidence from optimized builds, because broad data-layout changes
have already shown benchmark regressions.
