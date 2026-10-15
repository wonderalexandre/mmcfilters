# Pixel contours

`ContourComputation` extracts distinct foreground contour pixels for one node
or all live nodes. Use [Contour traces](contour-traces.md) for oriented edges
and ordered external and internal boundaries.

## Definition and input

A node's contour contains each support pixel with at least one side facing
the support complement or the exterior of the image. This is the foreground
A4 contour. Pixel identifiers refer to the original image grid in row-major
order; each pixel occurs once.

Membership uses four side neighbors regardless of construction adjacency.
For a tree of shapes, it uses the support projected onto the original grid.
The source must have a committed rooted topology, a non-empty regular 2D
domain, and one live smallest node for every pixel.

## C++ access and lifetime

```cpp
#include <mmcfilters/contours/ContourComputation.hpp>

using namespace mmcfilters;

ContourComputation contours(tree); // MorphologicalTree or current valued view

for (auto [nodeId, contourPixels] : contours) {
    // Borrowed pixels; valid until this iterator advances.
}

contours.forEachContour(
    [](NodeId nodeId, std::span<const PixelId> contourPixels) {
        // Consume or copy the pixels during this callback.
    });

std::vector<PixelId> nodeContour = contours.contour(nodeId);
```

Iteration and callbacks visit every live node in post-order: children precede
their parent. Pixel order and sibling order are unspecified. The C++20
`input_range` has single-pass iterators: copies share one position, and each
`begin()` starts an independent traversal. Its node sequence need not match
that of `ContourTraceComputation`.

An emitted span expires when its iterator or a copy advances, or when the
shared traversal is destroyed. Callback spans expire on return. To stop early,
use `break` or throw from the callback; exceptions propagate.

`contour(nodeId)` accepts a live internal `NodeId` and returns an owned vector
without invalidating active traversals. Repeated queries and traversals
recompute their results.

The tree must outlive the computation and its iterators; its topology and pixel
mapping must remain unchanged. Iterators retain indexes if the computation is
destroyed; owned vectors are independent of both. Rebuild the computation after
editing topology or pixel mapping. Checked access rejects stale computations;
borrowed spans perform no validity checks. Traversal buffers are released with
the last iterator sharing them.

## Python access

```python
contours = mmcfilters.ContourComputation(tree)

for node_id, contour_pixels in contours:
    process(node_id, contour_pixels)

contours.for_each_contour(process)
node_contour = contours.contour(node_id)
```

Every result is an owned NumPy array, including callback arguments. Each
emitted contour is copied once, so retained arrays survive iterator advancement
and destruction. The computation and its iterators keep the source valued
tree alive; the same topology and pixel-mapping requirements apply.

## Algorithm

`ContourLifetimeIndex` assigns each pixel an interval along the ancestor chain
of its smallest node. For an interior image pixel, the interval ends just
before the lowest common ancestor (LCA) of its smallest node and those of its
four side neighbors.
When both endpoints coincide, the pixel never belongs to a contour. Image-border
pixels remain on the contour through the root.

One batch of LCA queries to `MorphologicalTree` determines the endpoints.
Additions and removals are grouped by node. `NodeSupportIndex` stores every
pixel once and identifies each node support by a contiguous interval.

`ContourTraversal` processes the child with the largest support last, so the
parent can reuse its contour buffer and position table. It appends the other
child contours, releases their buffers, and applies additions and removals.
Removal swaps with the last pixel and updates its position in constant time.
`contour(nodeId)` instead scans the indexed support and tests pixel lifetimes.

The extinction contour map and exact distance transform consume this
traversal; the approximate distance transform shares the lifetime index.

## Cost and storage

Let `P` count image pixels, `N` internal node slots, `P_v` pixels in the support
of node `v`, and `C_v` its contour pixels. Let `C = sum_v C_v` over live nodes.

| Operation | Time | Additional storage |
| --- | --- | --- |
| Build indexes | `O(P + N)` plus LCA preparation | `O(P + N)` |
| Visit all contours | `O(N + P log(P + 1))`; reading or copying all pixels adds `O(C)` | `O(P + N)` per traversal |
| `contour(v)` after construction | `O(P_v)` | `O(C_v)` for the result |

The batch reuses an existing Euler tour and range minimum query (RMQ) cache.
Otherwise, it uses depth-first search (DFS) intervals for comparable pairs and
selects RMQ or offline Tarjan by estimated storage. Tarjan takes `O((P + N) alpha(N))` time,
where `alpha` is the inverse Ackermann function. RMQ preparation takes
`O(N log(N + 1))` time and space and is selected only when its estimate fits
the linear batch storage budget. See [Morphological trees](trees.md).

Each transfer from a child other than the one with the largest support places
a pixel in a support at least twice as large, giving at most `O(log P)` transfers
per pixel. Independent traversals share indexes but have separate working
buffers.

Bounds exclude source-tree construction, a pre-existing LCA cache, and retained
outputs. Keeping every contour requires `O(C)` storage, up to `O(PN)`.

## Related guides

- [Contour traces](contour-traces.md): ordered edges and boundaries.
- [Morphological trees](trees.md): supports, smallest nodes, and identifier domains.
- [Attribute catalog](attribute-catalog.md): `CONTOUR_*` scalar attributes.
- [Editing API](editing-api.md): mutations and derived-data lifetime.
