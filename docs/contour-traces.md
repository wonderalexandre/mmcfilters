# Contour traces

`ContourTraceComputation` represents node boundaries as oriented pixel sides and
traces those sides into ordered loops. Use it when an operation needs side-level
geometry, loop order, signed area, or external/internal boundary separation. For
compact boundary-pixel sets, use [Pixel contours](contours.md).

## C++ API

```cpp
auto traces = ContourTraceComputation::extract(tree);

for (ContourTraceEdge edge : traces.getEdges(nodeId)) {
    // edge.pixel is the row-major support-pixel ID.
    // edge.side is North, West, East, or South.
}

for (const ContourTraceLoop& loop : traces.getLoops(nodeId)) {
    if (loop.kind == ContourLoopKind::External) {
        // outer boundary
    }

    for (ContourTraceEdge edge : traces.getLoopEdges(loop)) {
        // ordered loop edge
    }
}
```

`getEdges(node)` materializes the unordered side boundary of one node.
`getLoops(node)` materializes edges when necessary, traces only that node, and
returns an owning `std::vector`. Retaining the returned loop metadata is safe
while other nodes are traced. `materializeAll()` traces every live node.

The C++ result references its source tree, which must outlive subsequent lazy
access.

## Python API

```python
traces = mmcfilters.ContourTraceComputation.extraction(tree)

edges = traces.getEdges(node_id)
loops = traces.getLoops(node_id)

for loop in loops:
    loop_edges = traces.getLoopEdges(loop)
```

Python node-local queries return lists. The returned `ContourTraces` object keeps
its source tree alive.

## Geometry convention

Each boundary primitive is one side of a support pixel. Image rows grow downward
and columns grow rightward. Directed edges are oriented with the support pixel
on their right.

With this convention:

- external loops have positive doubled signed area;
- internal loops have negative doubled signed area;
- `ContourLoopKind` is determined from the signed-area sign.

At grid vertices with several outgoing boundary edges, tracing uses the
deterministic priority right turn, straight, left, then back. This separates
diagonal-touching supports into stable loop components.

## Relation to scalar attributes

For a node `v`, let `B_v` be its number of traced boundary sides.

- `B_v == CONTOUR_PERIMETER(v)`.
- Directional side counts equal `CONTOUR_SIDE_NORTH`,
  `CONTOUR_SIDE_WEST`, `CONTOUR_SIDE_EAST`, and `CONTOUR_SIDE_SOUTH`.
- Projecting edges to unique pixel IDs equals the pixel set returned by
  `ContoursComputedIncrementally::getContour(v)`.

`BITQUADS_NUMBER_HOLES` can provide a consistency check for simple connected
supports, but its digital-connectivity convention is separate from loop
classification.

## Materialization and lifetime

Extraction stores compact local edge changes. Edge materialization combines the
missing descendants of a requested subtree and caches each final node boundary.
Loop tracing is a separate lazy step, so callers that need side geometry do not
pay for ordered loops.

The result captures the source topology mutation version. After a topology edit,
create a new trace result.

## Complexity

Let:

- `P` be the number of image pixels;
- `N` be the number of internal node slots;
- `M(S)` be the number of missing nodes in a requested subtree;
- `B(S)` be the number of boundary sides committed for that subtree;
- `D(S)` be the compact local edge changes read during materialization;
- `B_v` be the boundary-side count of node `v`.

With valid tree-query caches, extraction is `O(N + P)`. First edge
materialization of `S` is `O(M(S) + B(S) + D(S))`; the total compact delta count
is `O(P)` because each pixel has four sides. Tracing one materialized node is
`O(B_v)`.

Materializing and tracing every live node is output-sensitive:

```text
O(N + P + sum_v B_v)
```

Cached edge and loop storage is proportional to the materialized output.

## Related guides

- [Pixel contours](contours.md): compact unordered boundary pixels.
- [Attribute catalog](attribute-catalog.md): contour-side attributes.
- [Morphological trees](trees.md): support and proper-part ownership.
- [Editing API](editing-api.md): derived-state lifetime.
