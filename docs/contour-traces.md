# Contour traces

Use `ContourTraceComputation` when an operation needs the oriented contour edges
of a node or its ordered external boundaries and holes. Edge extraction and
boundary tracing are lazy: requesting unordered edges does not pay the cost of
ordering them. Use [Pixel contours](contours.md) when only the set of contour
pixels is required.

## C++ API

```cpp
#include <mmcfilters/contours/ContourTraceComputation.hpp>

using namespace mmcfilters;

ContourTraceComputation contourTraces(tree);

for (ContourEdge edge : contourTraces.edges(nodeId)) {
    // edge.pixel is a row-major support-pixel identifier.
    // edge.side is North, West, East, or South.
}

for (const ContourBoundary& boundary : contourTraces.boundaries(nodeId)) {
    if (boundary.kind == ContourBoundaryKind::External) {
        // External boundary.
    }

    for (ContourEdge edge : contourTraces.boundaryEdges(boundary)) {
        // Edges follow boundary order.
    }
}

ContourBoundary external = contourTraces.externalBoundary(nodeId);
for (ContourEdge edge : contourTraces.boundaryEdges(external)) {
    // Edges follow the unique external boundary in order.
}

for (PixelId pixel : contourTraces.boundaryPixels(external)) {
    // Equivalent to edge.pixel for each ordered edge.
}
```

`edges(node)` returns a borrowed range over the unordered contour edges of one
node. `boundaries(node)` traces that node and returns an independent vector of
`ContourBoundary` descriptors. Each descriptor records its kind, its edge range,
and its doubled signed area. `externalBoundary(node)` returns the unique external
boundary and rejects a node support with zero or multiple external boundaries.
`boundaryEdges(boundary)` returns the ordered edge range identified by a
descriptor from the same computation. `boundaryPixels(boundary)` projects that
range onto `edge.pixel` without allocating storage; the same pixel can occur
more than once when different edges occupy different sides of it. `traceAll()`
traces every live node.

The computation references its source tree. The tree must outlive the
computation and must remain unchanged while its results are used.

For a tree of shapes built with complementary adjacencies, construct the
computation from the valued view so it can determine whether each node is a
lower or upper shape:

```cpp
ContourTraceComputation contourTraces(valuedTree.asView());
```

The constructor copies the connectivity choice for every node and does not
retain the altitude buffer. Construct another `ContourTraceComputation` after an
altitude change that alters shape polarity. Construction from topology alone can
always provide `edges(node)`. If boundary tracing reaches a diagonal ambiguity
without enough information to select the foreground connectivity,
`boundaries(node)` throws `std::invalid_argument`.

## Python API

```python
contour_traces = mmcfilters.ContourTraceComputation(tree)

edges = contour_traces.edges(node_id)
boundaries = contour_traces.boundaries(node_id)
external = contour_traces.external_boundary(node_id)

for boundary in boundaries:
    boundary_edges = contour_traces.boundary_edges(boundary)
    boundary_pixels = contour_traces.boundary_pixels(boundary)
```

Python returns independent lists for edge and pixel queries. The
`ContourTraceComputation` object keeps its source tree alive.

The public geometry types are:

- `ContourSide`, which identifies the side of a support pixel;
- `ContourEdge`, which contains a pixel identifier and a side;
- `ContourBoundaryKind`, which distinguishes external and internal boundaries;
- `ContourBoundary`, which contains an ordered edge range and
  `doubled_signed_area` in Python or `doubledSignedArea` in C++.

## Geometry convention

Each contour edge is one side of a support pixel. Image rows increase downward
and columns increase to the right. A directed edge keeps the support pixel on
its right.

Under this convention, external boundaries have positive doubled signed area
and internal boundaries have negative doubled signed area. The sign determines
`ContourBoundaryKind`.

At a diagonal contact, tracing uses the foreground connectivity retained from
tree construction. A 4-connected foreground is paired with an 8-connected
background, so the successor turns right around the foreground pixel. An
8-connected foreground is paired with a 4-connected background, so the
successor turns left around the background pixel.

A shared construction adjacency applies to every node. Under complementary
topographic adjacencies, lower shapes use `minAdjacency` and upper shapes use
`maxAdjacency`; node and parent altitudes determine the shape polarity. A tree
of shapes built with `SelfDualSpanImmersion` uses a 4-connected foreground for
both lower and upper shapes after projection onto the original image grid. This
matches the scalar bitquad projection, and each foreground uses the
complementary 8-connected background.

Every directed edge has one successor determined by this connectivity. A
boundary closes when the successor is its first directed edge. A boundary may
visit the same vertex more than once, including its initial vertex, but each
edge belongs to exactly one boundary. A missing successor or an edge repeated
before closure is an error.

The API preserves its right-side orientation and the channel order `North`,
`West`, `East`, and `South`. The SIBGRAPI slides use left-side orientation and
the order `East`, `West`, `South`, and `North`. Reversing the traversal and the
sign used by the shoelace formula converts between the conventions. Both
conventions report a positive area for external boundaries.

## Relation to scalar attributes

For a node `v`, let `B_v` be the number of its contour edges.

- `B_v == CONTOUR_PERIMETER(v)`.
- Counts by side equal `CONTOUR_SIDE_NORTH`, `CONTOUR_SIDE_WEST`,
  `CONTOUR_SIDE_EAST`, and `CONTOUR_SIDE_SOUTH`.
- Projecting the edges to unique pixel identifiers produces the same pixel set
  as `ContourComputation::contour(v)`.

`BITQUAD_NUMBER_HOLES` provides a consistency check when its foreground
connectivity matches the connectivity used for tracing. Area sums and edge
coverage alone cannot validate that connectivity: an incorrect split at a
self-contact can preserve both values while creating a spurious hole.

## Storage and lifetime

Construction stores compact edge additions and removals for every node. The
first request for a node combines the missing descendants of its subtree and
caches their final edge sets. Tracing is a separate step, so an operation that
only needs unordered edges does not allocate ordered boundary descriptors.

The computation captures the topology mutation version. Create a new
computation after any topology edit.

`hasCachedEdges(node)` and `hasTracedBoundaries(node)` expose the lazy state in
C++. Python provides `has_cached_edges(node)` and
`has_traced_boundaries(node)`. `hasTracedAllBoundaries()` in C++, or
`has_traced_all_boundaries` in Python, becomes true after every live node has
been traced.

## Complexity

Let:

- `P` be the number of image pixels;
- `N` be the number of internal node slots;
- `M(S)` be the number of nodes without cached edges in a requested subtree `S`;
- `B(S)` be the number of contour edges cached for those nodes;
- `D(S)` be the number of compact edge changes read for those nodes;
- `B_v` be the number of contour edges of node `v`.

Construction takes `O(N + P)` time and performs no LCA queries. The first edge
request for subtree `S` takes `O(M(S) + B(S) + D(S))`; the complete delta store
is `O(P)` because each pixel has four sides.

Tracing a node takes `O(B_v)` time after initializing the direct vertex index,
or expected `O(B_v)` time with the sparse hash index. The direct index requires
`O(P)` initialization once, then resets only the vertices touched by each node.
Tracing every live node is output-sensitive:

```text
O(N + P + sum_v B_v)
```

Cached edge and boundary storage is proportional to the output requested so far.
Even `edges(root)` caches edge sets for all descendants, although it does not
trace their ordered boundaries. The API retains these results, so tracing all
nodes can store `O(sum_v B_v)` output in addition to `O(P + N)` working buffers.
An LCA cache on the same tree belongs to another operation and is not part of
this computation.

Profiling uses external timers or sampling tools. The library does not expose
profiling counters or a separate instrumented build mode.

## Related guides

- [Pixel contours](contours.md) covers compact unordered contour pixels.
- [Attribute catalog](attribute-catalog.md) defines contour side attributes.
- [Morphological trees](trees.md) defines supports and the mapping from pixels
  to their smallest nodes.
- [Editing API](editing-api.md) defines the lifetime of derived data.
