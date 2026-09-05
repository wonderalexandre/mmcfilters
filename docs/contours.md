# Pixel contours

Use `ContourComputation` to visit the foreground contour pixels of every tree
node or to compute the contour of one node. A complete traversal keeps
`O(P + N)` working storage and does not retain previously emitted contours. Use
[Contour traces](contour-traces.md) when an operation needs oriented edges,
ordered boundaries, signed area, or separation between external boundaries and
holes.

## Contour definition

For each live node, the contour contains every pixel in the node support that
has at least one side adjacent to the 4-neighbor complement or to the exterior
of the image. This is the foreground A4 contour. Pixel identifiers follow
row-major order in the original regular 2D domain.

The contour does not depend on the adjacency used to construct the tree. For a
tree of shapes, the node support and its contour are projected onto the original
image grid.

The input tree must have:

- a committed rooted topology;
- a non-empty regular 2D domain;
- one valid live smallest node for every pixel.

The topology must remain unchanged while a `ContourComputation`, one of its
iterators, or a borrowed contour span is in use. Checked builds reject access
after a topology change.

## C++ API

```cpp
#include <mmcfilters/contours/ContourComputation.hpp>

mmcfilters::ContourComputation contours(tree); // topology or valued view

for (auto [nodeId, contourPixels] : contours) {
    // contourPixels is valid until this iterator advances.
}

contours.forEachContour(
    [](mmcfilters::NodeId nodeId,
       std::span<const mmcfilters::PixelId> contourPixels) {
        // Copy contourPixels only when it must outlive this callback.
    });

std::vector<mmcfilters::PixelId> nodeContour = contours.contour(nodeId);
```

Iteration follows post-order, so every child is emitted before its parent.
Pixel order and sibling order are unspecified. The range satisfies the C++20
`input_range` requirements:

- copies of one iterator share a single-pass position;
- each call to `begin()` starts an independent traversal;
- advancing an iterator invalidates its previous span and spans obtained from
  copies of that iterator.

A caller can stop iteration with `break` or by throwing from the callback. The
working buffers are released when the last iterator for that traversal is
destroyed. Iterators keep the contour indexes alive if the
`ContourComputation` object is destroyed, but the source tree must outlive both.

`contour(node)` scans only the requested support and returns an independent
`std::vector`. Calling it does not invalidate an active traversal.

## Python API

```python
contours = mmcfilters.ContourComputation(tree)

for node_id, contour_pixels in contours:
    process(node_id, contour_pixels)

contours.for_each_contour(process)
node_contour = contours.contour(node_id)
```

Python returns an independently owned NumPy array for every contour, including
contours passed to callbacks. Retaining an array is therefore safe after the
iterator advances or is destroyed. The copy occurs once for each emitted
contour. The Python object and its iterators keep the source tree alive.

## Implementation

Each pixel enters the contour at its smallest node. An interior pixel leaves the
contour at the lowest common ancestor of its smallest node and the smallest
nodes of its four side neighbors. A pixel on the image boundary remains on the
contour through the root.

`ContourLifetimeIndex` stores these start and end events once. It submits all
required lowest common ancestor queries to `MorphologicalTree` as one batch.
The tree reuses an existing Euler tour and range minimum query cache when one is
available. Otherwise, it chooses between building that cache and running the
iterative offline Tarjan algorithm from their estimated storage requirements.
Equality and ancestor relations are resolved directly from depth-first search
intervals. See [Morphological trees](trees.md) for the batch LCA contract.

`ContourTraversal` selects the child with the largest support as the heavy
child. It processes that child immediately before its parent, allowing the
parent to reuse the child's contour vector. Contours from the other children
are transferred and released before the traversal applies the node's start and
end events. A position table supports constant-time removal by swapping with the
last element.

The exact distance transform uses the same traversal and adds only the contour
changes and heavy-path information required to update its distance field. The
approximate distance transform keeps its level schedule and shares the contour
lifetime index.

A call to `contour(node)` scans the node support and tests the lifetime of each
pixel. Repeated calls recompute the contour. The extinction contour map uses the
incremental traversal and preserves the priority of selected nodes where their
contours overlap.

## Complexity

Let `P` be the number of pixels, `N` the number of internal node slots, and
`C = sum_v |contour(v)|` the total number of emitted contour pixels.

- Index construction takes `O(P + N)` time when the tree already has an
  Euler/RMQ cache. Without that cache, it uses either RMQ construction or
  `O((P + N) alpha(N))` offline Tarjan. The strategy admits RMQ only when its
  estimated storage does not exceed the linear Tarjan buffers for this batch.
- A complete traversal takes `O(N + P log(P + 1))` time to maintain the contour
  sets. Reading or copying every emitted pixel adds `O(C)` time. Each transfer
  from a light child moves a pixel into a support at least twice as large, which
  bounds the number of transfers per pixel by `O(log P)`.
- Each active traversal uses `O(P + N)` working storage. Independent simultaneous
  traversals have independent working buffers and share the immutable indexes.
- `contour(v)` takes `O(|support(v)|)` time and returns
  `O(|contour(v)|)` owned storage after index construction.

The implementation does not cache every emitted contour. If the caller retains
all outputs, their storage is `O(C)` and can reach `O(PN)`. The working-storage
bound excludes retained outputs and any LCA cache that already belongs to the
tree.

## Related guides

- [Morphological trees](trees.md) defines the mapping from pixels to their
  smallest nodes and the `NodeId` domains.
- [Attribute catalog](attribute-catalog.md) defines the `CONTOUR_*` scalar
  attributes.
- [Editing API](editing-api.md) defines the lifetime of derived data after tree
  mutation.
