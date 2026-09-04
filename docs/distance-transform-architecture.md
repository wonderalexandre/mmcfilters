# Distance-transform architecture

This contributor guide documents the implementation behind the public scalar
contract in [Distance-transform attributes](distance-transform.md). The
distance-transform computers use the common scheduling and result-assembly path
described in [Attribute computer architecture](attribute-computer-architecture.md);
they do not define an alternate public API.

Both families operate on node supports and their foreground 4-connected
contours in the original image grid. The unsuffixed family uses the adaptive-A8
DIFT approximation, while the `_EXACT` family uses an exact separable squared
Euclidean distance transform. Neither backend dispatches from the descriptive
tree kind.

## Provider and reducer boundary

`NodeDistanceFieldProviderFor<Provider, Reducer>` separates numerical field
computation from scalar projection. A provider declares
`DistanceFieldAccuracy::Exact` or `DistanceFieldAccuracy::Approximate` and feeds
the `NodeDistanceTransformReducer` protocol. Selection is compile-time; the
public attribute path has no virtual dispatch or runtime backend branch.

Reducers materialize only requested scalars. Per-pixel distance fields and
public profiles are outside this API. An isolated maximum request therefore
uses a smaller policy than a request for distributional, localization, plateau,
or distance-weighted spatial descriptors.

## Exact backend

`ExactNodeDistanceFieldProvider` currently resolves to
`MorphologicalTreeDistanceTransform`. It evaluates every live node from its
support and foreground contour, using one implementation for max-trees,
min-trees, trees of shapes, residual trees, and generic trees.

### Support and contour indexes

`NodeSupportIndex` appends every proper-part pixel exactly once in
pre-order. The proper parts of every subtree consequently form a contiguous
support interval. Dense intervals remain indexed by internal `NodeId`,
including trees with empty proper parts or dead slots. The temporary
index costs `O(P + N)` time and memory and is invalidated by tree mutation.

`ContourLifetimeIndex` derives contour lifetimes without rescanning every
support. Let `o(p)` be the smallest node of pixel `p`, and let
`h(p)` be the LCA of `o(p)` and the smallest nodes of its four side neighbors.
An interior pixel is a contour site on the path from `o(p)` to the root, up to
but excluding `h(p)`. A pixel touching the image boundary remains active through
the root.

The index builds pre-order DFS intervals and submits the required LCA pairs as
one batch to `MorphologicalTree`. Construction takes `O(P + N)` time when the
tree already has an Euler/RMQ cache. Without that cache, the tree chooses
between building it and running `O((P + N) alpha(N))` iterative offline Tarjan.
DFS intervals resolve equal nodes and ancestor pairs directly. Tarjan receives
only the remaining incomparable pairs. None of these traversals uses recursion.

### Heavy-path transform schedule

`NodeSupportBoxIndex` adds the bounding boxes required by distance transforms.
`forEachContourUpdate` adapts the shared `contours::detail::ContourTraversal`,
which also implements public contour iteration. The traversal selects the child
with the largest support as the heavy child. It emits light subtrees first and
the heavy child immediately before its parent. The parent reuses the heavy
child's contour, transfers the other child contours, and applies the lifetime
events. Each transferred pixel enters a support at least twice as large, so it
is transferred at most `O(log P)` times.

Only the distance-transform adapter collects transition lists and heavy-path
tops; public pixel-contour iteration does not allocate those buffers.

Each heavy path uses the bounding box of its top node as one translated
transform domain. The leaf initializes the exact separable transform. A
heavy-child-to-parent update toggles contour sites, recomputes only affected
columns vertically, and recomputes only rows whose vertical-cost vector
changed. An unchanged support performs no line transform.

For a tree of shapes, support and contour are projected onto the original image
grid. The implementation does not define a distance transform on an immersed or
Khalimsky grid.

### Exact sample stream

The backend emits a borrowed `NodeDistanceTransformFrame` whose support span and
workspace remain valid only during synchronous iteration. One frame can feed
several reducers from the same squared-distance sample stream.

The maximum-only reducer keeps the isolated `MAX_DIST_EXACT` and
`MAX_SQUARED_DIST_EXACT` paths lightweight. Summary reducers can additionally
accumulate sample count, sums, sums of squares, `sum(sqrt(z))`, a sparse
histogram, and distance-weighted spatial raw moments. Separate public calls
still perform separate traversals.

Localization stores `(squaredDistance, pixel)` and resolves equal distances by
the smallest row-major `PixelId`. Plateau projection also stores the maximizer
count and coordinate sums. For `S` emitted samples, this projection costs
`O(S)` time and `O(1)` reducer state without changing the transform schedule.

The exact workspace stores signed 64-bit integer squared costs. Public buffers
use the requested floating-point result type, so `float32` can round squared
integers above `2^24`; use `double` in C++ or `dtype=np.float64` in Python when
integer-level preservation matters.

## Approximate backend

The unsuffixed family uses
`distance_transform_approx::EdtDIFT2D`. Foreground contour events remain exact
A4 events, while the field is propagated with the adaptive A8 stencil. One
global workspace processes topology-derived levels and inserts each proper-part
pixel once.

A structural preflight checks whether the smallest nodes of the pixels on every
A8 edge are comparable. If they are, all edges remain active and the hot kernel
has the unrestricted Opt3 form. Otherwise, an edge whose pixels have
incomparable smallest nodes is activated when their LCA is processed. This
generalizes the component-tree schedule to all supported hierarchies without
dispatching on a tree-kind label.

### Queue and validation policy

The production queue uses the PQueue32-style `GFT_FAST` intrusive-node layout
with FIFO ordering inside each integer squared-distance bucket. The DIFT
boundary establishes dense element and cost domains once, so `insert`,
`update`, `erase`, `contains`, and `pop` do not repeat caller-contract checks.
The constructor still rejects invalid allocation domains. Audited policies are
available to tests but are not selected through CMake or the public API.

### Approximate observers

`MorphologicalTreeApproximateDistanceTransform::forEachNodeMaximum` owns the
DIFT state and emits one squared maximum per live node. Attribute projection
either returns it as `MAX_SQUARED_DIST` or applies one square root for
`MAX_DIST`.

Other entry points select independent compile-time observers:

- union-find component moments for sums, means, RMS, and variances;
- ordered small-to-large sparse histograms for quantiles, entropy, positive
  area, and level count;
- additive spatial raw moments for distance-weighted geometry;
- `MaximumPixelTracker2D` for the canonical maximum pixel;
- `MaximumPlateauTracker2D` for maximum, canonical pixel, maximizer count, and
  coordinate sums.

Unrequested observers compile to no-op policies. Histogram-only requests do not
maintain spatial moments, spatial-only requests do not maintain histograms, and
maximum-only requests do not allocate localization or plateau state.

### Correspondence with the JMIV algorithm

The production implementation uses descriptive C++ names while preserving the
following correspondence with the JMIV 2025 pseudocode:

| JMIV state | C++ state | Role |
| --- | --- | --- |
| `bin` | `support_` | Pixels in an active binary component. |
| `root` | `root_` | Contour seed assigned to a support pixel. |
| `cost` | `cost_` | Squared distance to the assigned seed. |
| `open` | `open_` | Pixels whose labels may still change. |
| `adjmap` | `stencil_` | Adaptive A8 propagation/removal stencil. |
| `Bedt` | `rootMaximum_` | Maximum processed squared cost for a DIFT root. |
| `Q` | `queue_` | FIFO integer bucket queue. |
| removal stack `T` | `invalidationStack_` | Work stack for invalidating a removed seed's forest. |

`EdtDIFT2D::run` implements `EDTDiff`, `EdtDIFT2D::removeSeeds` implements
`treeRemoval`, `insertSupportNeighbours` performs A4 frontier insertion, and
`maximumRootDistance` performs the final `max Bedt` reduction. Conditional edge
activation and its LCA index are the topology-driven generalization and do not
come from the paper's altitude-level schedule.

## Correctness references

Two implementations are retained only under `unit-tests`:

- an independent brute-force Euclidean oracle reconstructs every small node
  support and contour;
- the former altitude-level DIFT sweep is an optimized component-tree oracle.

The brute-force oracle is the correctness authority for the exact family. The
component-tree oracle preserves the paper-style approximation for regression
and performance comparison but cannot overrule an exact discrepancy. Neither
oracle is installed, linked, or callable from production code.

## Complexity

For the approximate backend, level construction and the comparability preflight
cost `O(N + P)`. Trees whose adjacent pixels have comparable smallest nodes need
no activation index. The general path reuses the tree's Euler/RMQ LCA cache,
built in `O(N log N)` for a stable topology, performs each A8 smallest-node pair
query in `O(1)`, and stores `O(P)` edge-schedule state. Propagation remains
data-dependent: a local change may invalidate a large forest.

If `U` is the number of finite-cost assignments and `E` is the number of
activated domain edges, approximate summary observers add
`O((P + E + U) alpha(P))` time and `O(P)` memory. Plateau tracking adds
`O(Q + C)` work for `Q` popped labels and `C` reduced contour-root summaries.

For the exact backend, let `B(v)` be the area of node `v`'s minimal support box.
For heavy path `pi`, let its fixed box have height `H_pi` and width `W_pi`; let
`c_e` be the number of dirty columns and `r_e` the number of rows whose vertical
costs changed. Exact EDT line work is

`sum_pi 2*H_pi*W_pi + sum_heavy_edges (H_pi*c_e + W_pi*r_e)`.

The structural schedule adds `O(P log P + P + N)` when the tree already has an
Euler/RMQ cache. Without that cache, construction adds either the selected RMQ
index or `O((P + N) alpha(N))` offline Tarjan. Emitting all samples additionally
costs `sum_v |support(v)|`. Local working memory is `O(P + N)`. When the tree
selects RMQ, its estimated representation does not exceed the linear Tarjan
buffers for this batch and remains available to later tree operations.
There is no universal subquadratic claim: a chain of distinct full-span
supports has the output lower bound `Omega(P*N)`.

## Validation

Focused validation targets are:

```bash
cmake --build build --target \
  unit_maxdist_support \
  unit_distance_transform_study

ctest --test-dir build --output-on-failure -R \
  "unit_(maxdist_support|distance_transform_study)"
```

The scalar timing and sensitivity protocols are documented in the
[benchmark guide](https://github.com/wonderalexandre/mmcfilters/blob/main/benchmarks/README.md).
Runtime attribution must use the external native profiling workflow;
production kernels must not acquire timers, counters, or diagnostic branches
solely for profiling.

## Related guides

- [Distance-transform attributes](distance-transform.md): public names,
  mathematical definitions, units, and result conventions.
- [Attribute computer architecture](attribute-computer-architecture.md): common
  family registration, scheduling, numeric policy, and extension contract.
- [Pixel contours](contours.md): public contour semantics and extraction API.
