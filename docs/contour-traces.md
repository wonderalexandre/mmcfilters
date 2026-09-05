# Contour traces

`ContourTraceComputation` extracts contour edges and orders them into closed
boundaries for one node or all live nodes. It derives edges directly from node
supports, without taking a `ContourComputation` result as input. Use
[Pixel contours](contours.md) for distinct foreground contour pixels.

## Geometry and connectivity

A `ContourEdge` identifies a support pixel and a `ContourSide` facing the
support complement or image exterior: `North`, `West`, `East`, or `South`.
A `ContourBoundary` records the sequence's offset (`edgeOffset`), edge count
(`edgeCount`), `ContourBoundaryKind`, and doubled signed area (`doubledSignedArea`).
Coordinates refer to the original image grid, with rows increasing downward
and columns to the right. Directed edges keep the support pixel on the right.
`External` boundaries have positive signed area; `Internal` boundaries around
holes have negative signed area.

Each edge belongs to exactly one boundary. Consecutive edges share a vertex;
the last connects to the first without repeating it in the output. Vertices
and pixels can repeat. Boundary order and each boundary's starting edge are
unspecified.

At a diagonal contact, a 4-connected foreground uses an 8-connected background
and a right turn; an 8-connected foreground uses a 4-connected background and
a left turn. Connectivity determines boundary grouping; it does not change
the exposed pixel-side edges of a fixed support.

## C++ access and lifetime

```cpp
#include <mmcfilters/contours/ContourTraceComputation.hpp>

using namespace mmcfilters;

ContourTraceComputation contourTraces(tree);

for (auto [nodeId, trace] : contourTraces) {
    for (const ContourBoundary& boundary : trace.boundaries()) {
        for (ContourEdge edge : trace.boundaryEdges(boundary)) {
            // Consume edges in boundary order.
        }
    }
}

contourTraces.forEachTrace([](NodeId nodeId, ContourTraceView trace) {
    // Consume the borrowed view or copy it with ContourTrace(trace).
});

ContourTrace nodeTrace = contourTraces.trace(nodeId);
```

Iteration and callbacks visit live nodes in post-order. The C++20 `input_range`
has single-pass iterators: copies share one position, and each `begin()` starts
an independent traversal. A yielded `ContourTraceView` and its ranges expire
when that iterator or a copy advances, or when the traversal is destroyed.
Callback views expire on return.

`trace(nodeId)` accepts a live internal node and returns an owned `ContourTrace`.
`ContourTrace(view)` copies a borrowed result. Owned results can outlive the
computation and tree; their ranges borrow the result's storage. Repeated
queries and traversals recompute results. Node queries preserve active iterators.

| Method on a trace or view | Result |
| --- | --- |
| `boundaries()` | Borrowed span of boundary descriptors |
| `edges()` | All edges concatenated in boundary order |
| `externalBoundary()` | Unique external boundary; throws `std::logic_error` if there are zero or multiple |
| `boundaryEdges(boundary)` | Ordered edges selected by a descriptor from this trace |
| `boundaryPixels(boundary)` | One `edge.pixel` per ordered edge, including repetitions |

Edge and pixel ranges allocate no storage. Use `boundaryEdges` for a continuous
sequence: `edges()` can cross between separate boundaries. A projected support
can have multiple external boundaries. When there is exactly one:

```cpp
ContourBoundary external = nodeTrace.externalBoundary();
for (ContourEdge edge : nodeTrace.boundaryEdges(external)) {
    // Ordered external-boundary edges.
}
for (PixelId pixel : nodeTrace.boundaryPixels(external)) {
    // Different sides can yield the same pixel.
}
```

The source tree must outlive the computation and its iterators. Its topology
and pixel mapping must remain unchanged. Iterators retain indexes if the
computation is destroyed. Rebuild the computation after editing topology or
pixel mapping. Borrowed views do not check storage validity. Use `break` to
stop iteration; callback exceptions propagate and release traversal buffers.

## Input requirements

The source needs a committed rooted topology, a non-empty regular 2D domain,
and one live smallest node per pixel. A shared canonical 4- or 8-neighbor
construction adjacency applies to every node. For a tree of shapes with
complementary adjacencies, pass the valued view:

```cpp
ContourTraceComputation contourTraces(valuedTree.asView());
```

Lower shapes use `minAdjacency`; upper shapes use `maxAdjacency`. Altitudes
below or above the parent's determine polarity. Equal altitudes leave it
unresolved when these connectivities differ. With `SelfDualSpanImmersion`,
projected supports use a 4-connected foreground and an 8-connected background
for both polarities.

Construction copies the connectivity choices without retaining altitudes.
Rebuild the computation if an altitude change alters polarity. Reaching a
diagonal contact with unresolved connectivity throws `std::invalid_argument`.
Missing successors or edges revisited before closure are also errors.

## Python access

```python
contour_traces = mmcfilters.ContourTraceComputation(tree)

for node_id, trace in contour_traces:
    for boundary in trace.boundaries():
        boundary_edges = trace.boundary_edges(boundary)
        boundary_pixels = trace.boundary_pixels(boundary)

contour_traces.for_each_trace(process)  # process(node_id, trace)
root_trace = contour_traces.trace(tree.root)
external = root_trace.external_boundary()
```

Python accepts a valued tree. Iteration and callbacks copy each emitted view
into an owned `ContourTrace`; node queries also return owned traces.
`ContourTraceView` is C++ only. Edge, pixel, and boundary collections are
independent lists. `external_boundary()` returns one descriptor and raises
`RuntimeError` unless exactly one external boundary exists.

Descriptor fields use snake_case: `edge_offset`, `edge_count`, and
`doubled_signed_area`. The computation and iterators keep the tree alive;
retained traces and lists are independent of it.

## Algorithm and cost

An edge's lifetime starts at its support pixel's smallest node and ends just
before the LCA with its neighbor's smallest node. Empty lifetimes produce no
events; image-border edges remain through the root. Construction submits
horizontal and vertical pixel pairs to `MorphologicalTree`'s batch LCA.
`ContourEdgeDeltaStore` groups events by counting, prefix sums, and direct
placement. Each edge has at most one addition and removal; no deduplication
pass is needed.

`ContourTraceTraversal` reuses the child buffer with the most edges, merges
and releases the others, and applies the node's changes. `ContourBoundaryTracer`
uses a dense vertex index to follow successor edges into closed boundaries.
Scratch buffers are reused; per-edge geometry retains only the end vertex and
area contribution. Direction is decoded at diagonal contacts.

`trace(nodeId)` enumerates support through the node's subtree, tests neighbors
with ancestry intervals, and orders exposed edges with a sparse vertex index.

Let `P` count image pixels, `N` internal node slots, `P_v` pixels in the support
of node `v`, `T_v` its subtree's live-node count, and `B_v` its contour edges.
Let `B = sum_v B_v` over live nodes.

- Construction: `O(P + N)` time plus LCA preparation and `O(P + N)` storage;
  see [Pixel contours](contours.md) for batch LCA costs and cache exclusions.
- Full traversal: `O(P + N + B)` for assembly and ordering, plus `O(P R)` for
  membership-table resets. The 16-bit generation counter resets about once per
  65,535 visited nodes; `R` counts those resets. Each traversal uses `O(P + N)`
  working storage and shares the immutable indexes.
- `trace(v)` after construction: expected `O(T_v + P_v + B_v)` time and
  `O(T_v + B_v)` auxiliary storage, including the support iterator's stack.
  If ancestry intervals are absent, their first construction adds `O(N)` time
  and tree-owned storage. The result occupies `O(B_v)`.

Costs exclude source-tree construction and consumer work. Copying a trace costs
`O(B_v)`; retaining all traces requires `O(B)` storage. Python also allocates
and fills each requested list.

## Relation to scalar attributes

Edge counts equal `CONTOUR_PERIMETER`; side counts equal the corresponding
`CONTOUR_SIDE_NORTH`, `CONTOUR_SIDE_WEST`, `CONTOUR_SIDE_EAST`, and
`CONTOUR_SIDE_SOUTH` attributes. Deduplicating pixels from all edges gives the
same set as `ContourComputation::contour(v)`.

Summed `doubledSignedArea` values equal twice the support area. Internal-boundary
counts agree with `BITQUAD_NUMBER_HOLES` when connectivity matches. Edge coverage
and signed area alone cannot detect an incorrect pairing at a diagonal contact.

## Related guides

- [Pixel contours](contours.md): distinct foreground contour pixels.
- [Attribute catalog](attribute-catalog.md): contour and bitquad attributes.
- [Morphological trees](trees.md): supports, smallest nodes, and LCA queries.
- [Editing API](editing-api.md): mutations and derived-data lifetime.
