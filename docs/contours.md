# Pixel contours

`ContoursComputedIncrementally` represents each node boundary as a set of
support-pixel IDs. Use it for compact pixel contours and compatibility with
`CONTOUR_*` scalar attributes. Use [Contour traces](contour-traces.md) when an
operation needs oriented sides, ordered loops, or external/internal separation.

## Contract

For a node, the contour contains every support pixel that exposes at least one
side to the 4-neighbor complement. Pixel IDs are row-major indices in the
tree's regular 2D pixel domain.

The side relation is fixed to 4-neighbor geometry, independently of the
adjacency used to construct the morphological tree. This matches the public
`CONTOUR_PIXELS` and `CONTOUR_PERIMETER` definitions.

Extraction requires:

- a committed rooted topology;
- a non-empty regular 2D pixel domain;
- one valid direct smallest node for every pixel.

The result captures the tree mutation version and rejects reads after topology
mutation.

## C++ API

```cpp
auto contours =
    ContoursComputedIncrementally::extractCompactContours(tree);

for (int pixel : contours.getContour(nodeId)) {
    // use one contour pixel
}

for (auto [nodeId, contour] : contours.contoursByNode()) {
    for (int pixel : contour) {
        // use every live-node contour
    }
}
```

`getContour(nodeId)` returns a cache-aware range. The first iteration
materializes the requested subtree as needed; later iterations over an already
materialized node scan cached contiguous values.

`contoursByNode()` iterates live nodes and uses the same node-local
materialization path. `materializeAll()` is an optional prefetch for workloads
that will revisit many contours; it is not required for ordinary iteration.

## Python API

```python
contours = mmcfilters.ContourComputation.extraction(tree)

root_contour = list(contours.get_contour(tree.root))

for node_id, contour in contours.contours_by_node():
    pixels = list(contour)
```

Python uses `get_contour(node_id)`, `contours_by_node()`, and `materialize_all()` with
the same semantics as C++.

## Materialization and lifetime

Extraction stores compact node-local changes. Accessing a node combines missing
descendant results, applies the node's local additions and removals, and caches
the final unique pixel set. Materializing one subtree does not require
materializing unrelated branches.

A cached result references its source tree in C++. The source must therefore
outlive contour access. Python keeps the source tree alive through the returned
contour object.

After a topology edit, create a new contour result. Altitude-only changes do not
change support geometry; topology changes require a new result under the source
tree's mutation contract.

## Complexity

Let:

- `P` be the number of image pixels;
- `N` be the number of internal node slots;
- `M(S)` be the number of not-yet-materialized live nodes visited in subtree
  `S`;
- `C(S)` be the number of contour-pixel values committed while materializing
  `S`.

Lifetime extraction performs an iterative offline LCA pass in
`O((P + N) alpha(N))` time and `O(P + N)` memory, where `alpha` is the inverse
Ackermann function. It deliberately does not build the tree's persistent
`O(N log N)` Euler/RMQ LCA cache. Copying lifetime events into the compact
contour store is `O(P + N)` and does not change this bound.

The public upper bound for first materialization of `S` is
`O(M(S) + C(S) + P)`. The tighter output-sensitive form replaces `P` with the
compact local changes read for that subtree. Accessing an already materialized
node is `O(1)` before iteration; iterating its range is linear in that contour's
size.

Materialized storage is output-sensitive because each cached node owns its own
contiguous contour slice.

## Related guides

- [Contour traces](contour-traces.md): oriented boundary sides and loops.
- [Morphological trees](trees.md): smallest-node mapping and `NodeId` domains.
- [Attribute catalog](attribute-catalog.md): `CONTOUR_*` attributes.
- [Editing API](editing-api.md): lifetime after topology mutation.
