# Attribute computer architecture

This guide documents the attribute-computer architecture and the extension path
for adding or changing attributes. For public usage, see [Attributes](attributes.md);
for the attribute table, see [Attribute catalog](attribute-catalog.md).

The public API calls each scalar an attribute. In this guide, descriptor refers
only to the mathematical quantity represented by a public scalar attribute.

In this subsystem, incremental means bottom-up and finite-window-oriented
computation over the current tree. It does not mean that every public attribute
buffer stays live after arbitrary topology edits.

## Public boundary

Ordinary application code should include `mmcfilters/attributes/Attributes.hpp`
and call `AttributeComputation`. Concrete computers are advanced extension
components, not alternate public orchestration APIs.

The C++ library is header-only, so installed packages include some `detail`
headers as transitive implementation dependencies. Those headers are shipped so
public headers compile downstream; they are not compatibility-contract headers.
The explicit installation manifest is
`cmake/mmcfiltersPublicHeaders.cmake`; adding a repository header does not
implicitly publish or install it.

## Computer contract

An attribute computer owns one coherent attribute family. Every computer must
provide:

- `inline static constexpr familyName` for diagnostics;
- `inline static constexpr family` for scheduler grouping;
- `inline static constexpr domain` for execution routing;
- `inline static constexpr producedAttributes` as the canonical attribute list;
- `static compute(context)` for internal-node rows;
- `static computeUnitRows(unitContext)` for compact exported-Higra unit rows.

Computers are stateless static kernels. The produced-attribute list has a single
source of truth: the computer class. `runtimeProducedAttributes<Computer>()`
materializes it only for call sites that need runtime storage.

Unit-row support is mandatory. If an attribute has a degenerate one-pixel
meaning, the computer defines that value explicitly. Otherwise it still defines
the exported unit-row convention.

## Registry and metadata

`AttributeRegistry.hpp` stores public attribute metadata:

- public name;
- description;
- group membership;
- altitude requirement;
- topology/support eligibility (`topologyOnly` in the registry).

Group membership is metadata. Public requests may mix scalar attributes and
groups; the pipeline expands groups, deduplicates scalars, and returns only the
requested public attributes.

`AttributeComputerRegistry.hpp` defines the computer protocol and
`RegisteredAttributeComputers`. Produced attributes are declared only by
`Computer::producedAttributes`, and scheduler grouping is declared only by
`Computer::family`.

## Execution model

At a high level, a request follows this path:

```text
request -> expand groups -> validate support -> materialize dependencies
        -> compute buffers -> assemble requested result -> project if needed
```

`AttributeFamilyScheduler` adds hidden dependencies, groups attributes by
family, and preserves dependency order. The central executors are:

- `executeAttributeComputationPlan(...)` for altitude-aware requests;
- `executeTopologyAttributeComputationPlan(...)` for topology/support requests.

The internal orchestration path is `detail::AttributePipeline`;
topology/support families are delegated to `TopologyAttributeBackend`. New code
should extend this path instead of adding another top-level execution pipeline.

Dependencies are ordinary attribute results consumed by another computer. They
are passed as `DependencySourceT<Real>`, a non-owning pair of `AttributeNames`
and `const Real*`. Dependency buffers are reusable only when they contain the
requested attributes and use `NodeIdSpace::MorphologicalTree`.

Several computers use bottom-up accumulation: preprocess the current node, merge
children into the parent, then finalize the current node. Increment-augmented public
calls compute the base attribute first, then materialize ancestor/descendant
sample offsets from a typed positive altitude step, sampling radius,
representative-descendant policy, and missing-sample policy.

## Hierarchical distance-transform backend

`MAX_DIST_EXACT` is a topology/support family. Its production path names the
`ExactNodeDistanceFieldProvider` explicitly. That provider is currently
`MorphologicalTreeDistanceTransform`, which evaluates every live node from its support
and foreground 4-connected contour in the original 2D domain. Its reusable
workspace applies the exact separable squared Euclidean transform along a
heavy-path schedule. Max-trees, min-trees, trees of shapes, residual trees, and
generic trees take this same path. There is no runtime tree-kind selection and
no fallback to another algorithm.

At the start of the traversal, `MorphologicalTreeRegionIndex` appends every
proper-part pixel exactly once in pre-order. The proper parts of every subtree
therefore form one contiguous support interval. Dense interval and bounding-box
metadata remain indexed by internal `NodeId`, including trees with empty proper
parts or dead slots. This temporary index costs `O(P + N)` time and memory and
is invalidated by a tree mutation.

`MorphologicalTreeBoundaryLifetimeIndex` derives the complete contour schedule
without rescanning every support. Let `o(p)` be the inclusion-smallest owner of
pixel `p`. For an interior-domain pixel, let `h(p)` be the LCA of `o(p)` and the
owners of its four side-neighbours. Then `p` is a contour site exactly on the
owner-to-root path from `o(p)` up to, but excluding, `h(p)`. A pixel touching the
global image boundary stays active through the root because its exterior
neighbour never enters a support. Thus every pixel contributes at most one
contour-addition event and one contour-removal event. This statement depends
only on support inclusion, not altitude or tree kind.

The index does not use the tree's persistent Euler/RMQ LCA cache or lazy DFS
interval cache. It owns iterative DFS-preorder intervals, reduces each
five-owner LCA to the LCA of the first and last owner in preorder, and resolves
all resulting pairs in one iterative offline Tarjan pass. Boundary-lifetime
construction therefore uses `O((P + N) alpha(N))` time and `O(P + N)` memory,
without recursive descent.

`MorphologicalTreeContourScheduler` chooses the support-largest child as the
heavy child. Light subtrees are emitted first and the heavy child immediately
precedes its parent. The parent destructively reuses the heavy contour, transfers
light contours, and applies the lifetime events. A transferred pixel enters a
support at least twice as large, so it is transferred at most `O(log P)` times.
Resident contour storage remains `O(P + N)`.

Each heavy path uses the bounding box of its top node as one fixed translated
transform domain. Translation is an isometry and that box contains every
support and contour on the path. The path's leaf initializes the exact separable
transform once. Moving from the heavy child to its parent applies the exact site
insertions/removals: only affected columns are recomputed vertically, and only
rows whose vertical-cost vector changed are recomputed horizontally. If a
parent and child have identical supports, the empty update performs no line
transform. No approximate propagation, DIFT state, or alternative production
kernel is involved.

For a tree of shapes, the region is the node support projected onto the original
image grid and the boundary is its foreground 4-connected contour there. This
contract does not claim a distance transform on an immersed or Khalimsky grid.

The transform keeps squared costs internally. If `d2(N)` is the maximum emitted
cost, `MAX_SQUARED_DIST` and `MAX_SQUARED_DIST_EXACT` expose `d2(N)` in squared
pixels, while `MAX_DIST` and `MAX_DIST_EXACT` expose `sqrt(d2(N))` in pixels.
The maximizer, canonical center, and maximum plateau are identical under both
projections.
The transform workspace stores exact signed 64-bit integer costs. Public
attribute buffers use the requested floating-point result type, so the default
`float32` output can round squared integers above `2^24`; request `double` in C++
or `dtype=np.float64` in Python when integer-level preservation matters.

`NodeDistanceFieldProviderFor<Provider, Reducer>` is the policy boundary between
numerical field computation and attribute projection. A provider declares
`DistanceFieldAccuracy::Exact` or `Approximate` and feeds the same
`NodeDistanceTransformReducer` protocol. This is compile-time selection; there
is no virtual dispatch or runtime backend branch in the public attribute path.

The exact backend emits a borrowed `NodeDistanceTransformFrame` while the workspace is
valid. A frame exposes per-pixel squared-distance labels through synchronous
iteration and can fan one sample stream out to several
`NodeDistanceTransformReducer` objects. The lightweight max-only reducer remains
the path for an isolated `MAX_DIST_EXACT` or `MAX_SQUARED_DIST_EXACT` request;
the former performs one square root per node and the latter exposes the integer
maximum directly. If any distributional summary is
requested, one summary reducer accumulates the maximum, sample count, sum, and
sum of squares from the same stream and materializes `DIST_SQUARED_SUM_EXACT`,
`DIST_SQUARED_MEAN_EXACT`, `DIST_RMS_EXACT`, and `DIST_SQUARED_VARIANCE_EXACT` as requested. No
additional contour reconstruction, transform, or support walk is introduced by
requesting several summaries together. Separate public calls still compute
separate traversals by design. The frame is deliberately non-copyable because
its support span and workspace are reused for the next node.

The exact summary reducer can additionally accumulate `sum(sqrt(z))`, a sparse
histogram, and distance-weighted spatial raw moments from that same sample
stream. These internal summaries materialize only the requested scalar
attributes; per-pixel fields and public profiles are outside this API.

Exact maximum localization uses the same stream. A localized reducer stores
`(squaredDistance, pixel)` and resolves equal distances by the smallest
row-major `PixelId`, then materializes zero-based row and column attributes.
A plateau reducer additionally stores the maximizer count and row/column sums.
It resets that state on a strictly larger sample and accumulates equal samples,
which yields plateau area and centroid without another support traversal. The
original max-only reducer remains selected when neither localization, plateau
geometry, nor distributional moments are requested.
For `S` samples emitted by the exact provider, plateau projection is `O(S)`
time with `O(1)` reducer state; it does not change the provider's transform or
contour-scheduling complexity.

The correctness references are test-only:

- an independent brute-force Euclidean oracle reconstructs every small node
  support and contour;
- the former altitude-level DIFT sweep is retained as an optimized oracle for
  max-trees and min-trees only.

The brute-force implementation is the correctness authority. The optimized
component-tree oracle preserves the paper's adaptive DIFT for regression and
performance comparison; because that IFT is approximate in rare pixel cases,
it cannot overrule a discrepancy with the exact oracle.

Neither oracle is installed, linked, or callable from production code. Further
hierarchical scheduling improvements must remain internal optimizations of the
same uniform exact backend and preserve results against the brute-force oracle.

The unsuffixed `MAX_DIST` family preserves the historical paper-style
approximate contract. Its production engine is the adaptive-A8
`distance_transform_approx::EdtDIFT2D`; foreground contour events remain exact
A4 events. The separate `_EXACT` family uses the exact provider described
above. One global approximate workspace processes incomparable nodes in
topology-derived levels. Each proper-part pixel is inserted once. A structural
preflight checks whether the owners of every A8 edge are comparable. If so,
all edges remain active and production uses the unrestricted Opt3-style hot
kernel. Otherwise, an edge between incomparable pixel owners is activated only
when their lowest common ancestor is processed. This is a topology-derived
dispatch shared by max-trees, min-trees, trees of shapes, residual trees, and
generic trees; it does not select behavior from the declared tree kind.

The approximation also has an explicit provider/consumer boundary.
`MorphologicalTreeApproximateDistanceTransform::forEachNodeMaximum` owns all
DIFT state and emits one squared local maximum per live node. The attribute
projection either exposes that cost as `MAX_SQUARED_DIST` or applies one square
root for `MAX_DIST`; this remains the maximum-only path. `forEachNodeSummary` additionally exposes additive moments of
the current approximate cost field. A union-find observer tracks active binary
components; inserting a support pixel or activating an A8 edge merges component
moments, while every finite DIFT cost assignment removes the old contribution
and adds the new one. The attribute computer therefore obtains sum, mean, RMS,
and population variance without reconstructing or rescanning node supports.
Real-distance sums and variance extend the same additive moment observer.
Requests for quantiles, entropy, positive area, or level count select a
deterministic sparse squared-distance histogram. Distance-weighted geometry
selects additive spatial raw moments. `forEachNodeSelectedStatistics` combines
scalar moments, histogram, and spatial moments as independent compile-time
policies, so each request retains exactly the sufficient statistics it needs.
Histogram merges use ordered small-to-large component storage. Histogram-only
scalar requests therefore perform no spatial-moment updates, while spatial-only
requests perform no histogram updates. The max-only and
moment-only entry points remain specialized policies as well.
The standalone `computeMaxSquaredDistance` facade is only a materializing
consumer of the same provider. Thus storage policy and attribute projection are
isolated from propagation and contour scheduling.

`forEachNodeExtremum` and `forEachNodeSummaryAndExtremum` expose the same Bedt
maximum together with its attaining support pixel. `MaximumPixelTracker2D`
stores one deterministic pixel per active Bedt root and is updated whenever the
root maximum is processed. The tracker is a compile-time policy; max-only and
moment-only traversals instantiate `NoopMaximumPixelTracker2D`, so existing
requests allocate no maximum-pixel vector and perform no tie comparisons.

`forEachNodePlateau` and `forEachNodeSummaryAndPlateau` select
`MaximumPlateauTracker2D`. For each active `Bedt` root, the tracker stores the
maximum, canonical pixel, maximizer count, and coordinate sums. Node output
combines only root summaries whose `Bedt` value equals the node maximum. This
adds constant work per processed DIFT label and one fixed-size summary per
pixel only when a plateau attribute is requested; max-only, moment-only, and
localization-only calls retain their previous policies and storage.
If `Q` labels are popped and `C` contour-root summaries are reduced, the added
plateau work is `O(Q + C)` and the dense auxiliary state is `O(P)`. These terms
do not change the approximate traversal's asymptotic bound, but they do add a
measurable constant factor when plateau output is enabled.

The DIFT implementation keeps descriptive C++ names while documenting one
authoritative correspondence with the JMIV 2025 pseudocode:

| JMIV state | C++ state | Role |
| --- | --- | --- |
| `bin` | `support_` | Pixels currently present in an active binary component. The global topology-driven workspace may contain several incomparable components separated by inactive edges. |
| `root` | `root_` | Contour seed assigned to each support pixel. |
| `cost` | `cost_` | Squared Euclidean cost to the assigned seed. |
| `open` | `open_` | Membership bitmap for pixels whose labels may still change. |
| `adjmap` | `stencil_` | Index of the adaptive A8 propagation/removal stencil. |
| `Bedt` | `rootMaximum_` | Maximum processed squared cost associated with each DIFT root. |
| `Q` | `queue_` | Production PQueue32-style FIFO integer bucket queue. |
| removal stack `T` | `invalidationStack_` | Work stack used while invalidating a removed seed's forest. |

Likewise, `EdtDIFT2D::run` implements `EDTDiff` from Algorithm 1,
`EdtDIFT2D::removeSeeds` implements `treeRemoval` from Algorithm 2,
`insertSupportNeighbours` performs the A4 frontier insertion from Algorithm 3,
and `maximumRootDistance` performs its final `max Bedt` contour reduction.
The conditional `activeEdges_` state and LCA activation index are intentionally absent from
this table's JMIV column: they implement the topology-driven generalization,
not the paper's altitude-level schedule.

Level construction and the A8 comparability preflight cost `O(N + P)`. Every
support pixel is inserted once. Comparable-owner trees require no activation
index. On the general path, the tree's shared Euler/RMQ LCA cache costs
`O(N log N)` once per stable topology, each A8 owner query is `O(1)`, and the
edge schedule occupies `O(P)` additional memory. Remaining work depends on
contour-root invalidations and adaptive IFT propagation; a local change can
still invalidate a large forest in the worst case. The adaptive stencil may
select a non-global nearest contour root in rare discrete Voronoi
configurations, so equality with exact `MAX_DIST_EXACT` is not a public
guarantee.

Let `U` be the number of finite-cost assignments performed by that DIFT work
and `E` the number of activated domain edges. Enabling approximate summaries
adds `O((P + E + U) alpha(P))` time for union-find/moment events and `O(P)`
memory. This does not change the DIFT's propagation-dependent worst case. A
compile-time no-op observer is used for either maximum-only request, so it does not pay
for component moment maintenance. Exact summaries add constant work per
already-emitted support sample and three accumulator scalars per active
reducer; their asymptotic cost is unchanged from the complete exact sample
stream.

Maximum localization adds `O(1)` work per processed exact EDT sample or DIFT
queue pop and `O(1)` exact reducer state. The approximate tracker additionally
uses `O(P)` dense storage. Contour reduction already scans the active Bedt
roots, so returning the pixel with the maximum does not add a support scan.

Let `B(v)` denote the area of node `v`'s minimal support box. For a heavy path
`pi`, let its fixed top box have height `H_pi` and width `W_pi`. Let `c_e` be the
number of dirty columns in a heavy update and `r_e` the number of rows whose
vertical costs actually change. The exact EDT line work is

`sum_pi 2*H_pi*W_pi + sum_heavy_edges (H_pi*c_e + W_pi*r_e)`.

The structural schedule adds `O(P log P + (P + N) alpha(N))`. The shared reducer
stream additionally requires `sum_v |support(v)|` samples. Working memory is
`O(P + N)`. Deterministic statistics distinguish logical support samples from
actually consumed samples and expose path initializations, site toggles, dirty
columns/rows, contour transfers, and the former `sum B(v)` / `P*N` baselines.

There is no universal subquadratic claim. The complete sample stream has the
output lower bound `Omega(sum_v |support(v)|)`, which is `Omega(P*N)` for a
chain of distinct full-span supports. A contour-site change may also alter an
entire column and every horizontal row. Heavy scheduling therefore improves
reuse and practical work while preserving the same worst-case bound.

## Contexts and concepts

The context types in `AttributeKernelSupport.hpp` are the adapter boundary:

- topology/support node rows use `AttributeComputeContext<Real>`;
- altitude-aware node rows use `AltitudeAttributeComputeContext<Real, T>`;
- topology/support unit rows use `UnitAttributeComputeContext<Real>`;
- altitude-aware unit rows use `AltitudeUnitAttributeComputeContext<Real, T>`.

`TopologyAttributeComputer` and `AltitudeAttributeComputer` enforce the standard
computer protocol. A new family should not add public family-specific method
names. Private helpers and `detail` kernels may keep narrower signatures when
that makes implementation or testing clearer.

The generic finite-window extension API exposes role-typed `EventDelta`,
`LocalAttributeIncrement`, and `NodeAttribute` records so its mathematical
pipeline can be inspected and tested. `ConnectedSubsetTreeLocalizer` owns the
tree-dependent join `LCA(P(anchor), P(sample))`; `LocalDecision` owns only the
state-to-value map; `EventAlgebra` owns identity, addition, and subtraction;
and `FiniteWindowLocalEventCompiler` turns those policies into sparse node
events. Attribute-family storage such as
`detail::BitquadFamilyCounts` and `detail::ContourSideCounts` remains an
implementation detail and is not part of the public attribute-computer
contract.

Bitquad-family counting is hierarchy-independent. Connectivity-dependent scalar
formulas are a later materialization step and receive an explicit
`detail::BitquadConnectivityPolicy`. For a tree of shapes with unequal
complementary connectivity, the policy consumes `ShapePolarity::Lower` or
`ShapePolarity::Upper` derived from exact node-versus-parent altitudes. The root
has no polarity and is handled by the policy's separate root entry.

## Numeric policy

Computers use `AttributeNumericPolicy.hpp` for degenerate divisions, square
roots, non-negative clamping, finite fallbacks, and ratio bounds. Attribute
buffers should not expose accidental `NaN` or infinite values for ordinary
finite inputs.

## Adding or changing attributes

Start by deciding whether the attribute belongs to an existing family or
requires a new family. Prefer an existing family when traversal, dependencies,
or intermediate state are shared.

Common metadata steps:

1. Add the scalar enum in `AttributeTypes.hpp` and one matching row in
   `AttributeRegistry.hpp`.
2. Classify the attribute as topology/support, altitude-aware,
   adjacency-dependent, or tree-kind specific.
3. Add it to a group only when the group semantics still hold.

For an attribute in an existing family:

1. Add it to the family's `producedAttributes`.
2. Extend request selection and `compute(context)`.
3. Use `DependencyResolver<Real>` for semantic dependencies and
   `AttributeNumericPolicy.hpp` for finite fallbacks.
4. Add attribute-level dependencies in `AttributeFamilyScheduler.hpp` only when
   another materialized attribute is consumed.
5. Extend `computeUnitRows(unitContext)`.
6. Add focused value tests, plus plumbing tests when registry, dependencies,
   projection, or public layout changes.

For a new family:

1. Add a new `AttributeComputerFamily` value.
2. Create a computer under `mmcfilters/attributes/computers/`.
3. Declare `familyName`, `family`, `domain`, and `producedAttributes`.
4. Implement `compute(context)` and `computeUnitRows(unitContext)`.
5. Register the computer in `RegisteredAttributeComputers`.
6. Register execution in `AttributePipeline.hpp` or `TopologyAttributeBackend.hpp`.
7. Add contract/plumbing tests and focused value tests.

When the public surface or attribute semantics change, update Python bindings
as needed and keep [Attributes](attributes.md) and the
[Attribute catalog](attribute-catalog.md) synchronized with the registry.

## Validation

Useful checks while changing this subsystem are:

```bash
cmake --build build --target \
  unit_public_attribute_api \
  unit_attribute_plumbing \
  unit_attribute_unit_values \
  unit_attributes_on_morphological_tree \
  unit_finite_window_local_attribute_computations \
  unit_maxdist_support

ctest --test-dir build --output-on-failure -R \
  "unit_(public_attribute_api|attribute_plumbing|attribute_unit_values|attributes_on_morphological_tree|finite_window_local_attribute_computations|maxdist_support|installed_consumer)"
```

Run Python tests when bindings or the Python facade change.

## Non-goals

- Do not introduce virtual attribute-computer classes.
- Do not move away from dense node ID buffers.
- Do not add runtime polymorphism to hot attribute kernels.
- Do not expose attribute-family-specific finite-window storage as public
  computer API; use the generic role-typed pipeline for extension work.
