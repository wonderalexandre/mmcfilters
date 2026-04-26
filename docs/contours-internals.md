# Incremental Contours Internals

This document describes the internal design of
`mmcfilters/contours/ContoursComputedIncrementally.hpp` and its helper storage
types under `mmcfilters/contours/detail/`.

It is intended for maintainers. For the public API, see
[contours-api.md](contours-api.md).

## Design Goal

The contour computation has two separate responsibilities:

- extract local contour events for each tree node;
- materialize final node contours only when a caller actually iterates them.

The implementation avoids keeping two public access models. All contour reads go
through `getContour(node)` or `contoursByNode()`, and both use the same
incremental materialization path. `materializeAll()` is only an explicit prefetch
for workloads that know they will visit many or all contours.

## Files

- `ContoursComputedIncrementally.hpp`
  Public API, extraction traversal, lazy materialization, cache state, and
  diagnostics.
- `detail/PendingPixelLists.hpp`
  Transient per-node linked lists backed by one reusable contiguous buffer.
- `detail/ContourDeltaStore.hpp`
  Persistent compact delta store with CSR-like spans for additions and removals.

The `detail/` headers are implementation details. They should not be treated as
stable public API.

## Data Model

After `extractCompactContours(tree)`, each node has two local delta lists:

- additions: pixels that enter the node contour locally;
- removals: pixels that must be removed from the contour accumulated from
  descendants.

These deltas are stored in `ContourDeltaStore`:

```text
addValues      = all compacted addition pixels
addSpans[node] = offset/size into addValues

removeValues      = all compacted removal pixels
removeSpans[node] = offset/size into removeValues
```

The final materialized contour cache is stored separately:

```text
cachedContourValues_       = all materialized contour pixels
cachedContourOffset_[node] = offset into cachedContourValues_
cachedContourSize_[node]   = number of pixels for node
cachedContourReady_[node]  = whether offset/size are valid
```

This keeps extraction compact and lets broad iteration become incremental:
materializing one node also materializes any missing descendants in that
subtree.

## Extraction Pipeline

`extractCompactContours(tree)` performs one post-order traversal of the tree.

For each processed node:

1. It consumes pending removals that have reached the node.
2. It forwards removals to a stricter ancestor or to the lowest common ancestor
   when neighbor pixels belong to incomparable components.
3. It scans the node proper parts in the image domain using 4-neighbour side
   adjacency.
4. It records local contour additions and local removals in `PendingPixelLists`.
5. At the end of traversal, it compacts the transient lists into
   `ContourDeltaStore`.

The contour side adjacency is intentionally fixed to radius `1.0`, independently
from the adjacency used to build the morphological tree. This follows the
side-contour definition where exposed pixel sides are counted in the original
image domain.

## Materialization Pipeline

`getContour(node)` returns a `ContourRange`. The range materializes its node on
first `begin()`/`end()` use.

Materialization is handled by `ensureSubtreeMaterialized(root)`:

1. It walks the requested subtree in post-order with an explicit stack.
2. Already materialized children are reused and not traversed again.
3. For each missing node, it accumulates materialized child contours.
4. It adds the node local additions.
5. It applies the node local removals.
6. It commits the resulting compact slice into `cachedContourValues_`.

`contoursByNode()` is only an iterator over live nodes. It does not implement a
second contour model; each yielded contour range still uses `getContour`'s
materialization mechanism.

`materializeAll()` calls the same path from the root. It is a prefetch command,
not a different representation.

## Scratch Marking

Both delta compaction and contour materialization remove duplicates with a
generation-marked pixel array:

```text
pixelMark[pixel] == currentGeneration means "present in the current temporary set"
```

The generation counter is `uint16_t` to reduce memory. When it wraps to zero, the
whole mark buffer is cleared and the generation restarts from one. This makes the
common reset path O(1) and keeps the rare wraparound path correct.

The mark buffer is scratch state. It is not part of the materialized result.

## Invariants

The implementation relies on these invariants:

- `localDeltas_` is immutable after extraction.
- `cachedContourReady_[node] == 1` means `cachedContourOffset_[node]` and
  `cachedContourSize_[node]` identify a valid slice in `cachedContourValues_`.
- `getContour(node)` and `isContourMaterialized(node)` only accept live internal
  node ids.
- `PendingPixelLists` is transient. It is consumed or compacted before the
  returned `IncrementalContours` object is created.
- `ContourDeltaStore` stores compacted per-node deltas; duplicates inside one
  node span are removed during compaction.
- Pixel ids are image-domain linear indices. Invalid negative or out-of-domain
  pixels are ignored when compacting or materializing.

## Complexity

Let:

- `P` be the number of image pixels;
- `N` be the number of internal node slots;
- `E4` be the number of 4-neighbour pixel-side checks, bounded by `O(P)`;
- `D` be the total number of compacted local delta values;
- `C(S)` be the total number of cached contour values committed while
  materializing a subtree `S`.

Extraction:

```text
O(P + E4 * A + D)
```

`A` is the cost of ancestry/comparability/LCA queries on the current
`MorphologicalTree` implementation. With the Euler/RMQ LCA support, LCA is
constant time after preprocessing; ancestry/comparability still depend on the
tree query implementation.

First materialization of a subtree:

```text
O(C(S) + local_delta_values_in_S)
```

Repeated access to an already materialized node:

```text
O(contour_size(node))
```

`materializeAll()`:

```text
O(total_cached_contour_values + D)
```

The cost can be high for trees where many large contours are materialized,
because the cache stores each node contour as its own contiguous slice.

## Space Usage

Before materialization, dominant storage is:

```text
addValues + removeValues + addSpans + removeSpans + pixelMark_
```

After broad materialization, dominant storage is usually:

```text
cachedContourValues_
```

For large images, changing pixel ids from `int` to `uint16_t` is not generally
valid because image-domain indices exceed 65535. `uint32_t` would have the same
size as `int` on the supported platforms, so it does not materially reduce the
main cache. The current memory reduction comes mainly from using `uint16_t` for
the scratch mark buffer.

## Benchmark Cases

`examples/contour_benchmark.cpp` reports these workloads:

- `extract only [extractCompactContours]`
  Measures extraction and delta compaction without materializing final contours.
- `extract + root via getContour(root)`
  Measures extraction plus materializing the root subtree through normal access.
- `extract + iterate all via getContour(node)`
  Iterates every live node in tree order through `getContour`.
- `extract + iterate all via contoursByNode()`
  Uses the all-node range wrapper. It should be close to tree-order
  `getContour`.
- `extract + iterate all via getContour(node) random order`
  Measures broad access when cache reuse is less aligned with tree order.
- `extract + materializeAll() + iterate via getContour(node)`
  Separates prefetch/materialization from later iteration.

When the workload really needs all contours, `materializeAll()` is the clearest
way to express that intent. When the workload may stop early or visit only a
subtree, plain `getContour(node)` keeps the computation incremental.

## Maintenance Notes

Prefer preserving the single public access model:

```cpp
contours.getContour(node)
contours.contoursByNode()
contours.materializeAll()
```

Avoid reintroducing separate "on-demand" APIs unless they expose a genuinely
different semantic contract. The current implementation already computes on
demand and caches incrementally.

If memory becomes the dominant problem, benchmark before changing storage. The
largest post-materialization cost is usually the repeated materialized contour
cache, not the delta store. Delta compression such as varint encoding may reduce
memory, but it would add decoding cost and complicate random access; it should be
guarded by benchmarks on Component Tree and Tree of Shapes.
