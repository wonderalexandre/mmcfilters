# Morphological Trees

This guide describes the public tree model used by `mmcfilters`. It is the
orientation document for `MorphologicalTree`, `WeightedMorphologicalTree<T>`,
`WeightedTreeView<T>`, construction factories, node-id domains, edits, and Higra
interoperability.

## Public Model

The public tree API separates three concerns:

- `MorphologicalTree` owns topology only.
- `WeightedMorphologicalTree<T>` owns topology plus a dense altitude buffer.
- `WeightedTreeView<T>` borrows topology plus an external altitude span.

`MorphologicalTree` stores:

- the rooted parent/child topology;
- dense internal `NodeId` slots;
- direct proper-part ownership;
- optional regular 2D proper-part-domain metadata;
- optional immutable regular-grid 2D adjacency capability.

It does not own node altitudes. Any operation that reads or mutates altitude
must use `WeightedMorphologicalTree<T>` or a read-only `WeightedTreeView<T>`.

## Capabilities And Descriptive Kinds

`MorphologicalTree` describes interpretation requirements with
`HierarchySemantics`:

- `AltitudeOrder::INCREASING_FROM_ROOT` declares
  `altitude(parent) < altitude(child)`,
  `AltitudeOrder::DECREASING_FROM_ROOT` declares
  `altitude(parent) > altitude(child)`, and `UNCONSTRAINED` declares no global
  parent-to-child altitude order;
- `AdjacencyMode` declares no adjacency, one uniform adjacency, or directional
  decreasing/increasing adjacencies;
- `MorphologicalTreeKind` is a descriptive label, not the
  algorithm-dispatch contract. `GENERIC` is available when no named family is
  appropriate.

The factory assigns the expected capabilities to max-trees, min-trees, and
Trees of Shapes. Generic algorithms query `getAltitudeOrder()`,
`getAdjacencyMode()`, `hasUniformGridAdjacency2D()`, and
`hasDirectionalGridAdjacency2D()`. Uniform relations are read through the
immutable `getUniformGridAdjacency2D()` pointer. The older
tree-family-specific adjacency spellings are no longer part of the API.

Component trees require an image-domain adjacency relation. On the 2D square
grid, `radius=1.0` corresponds to 4-connectivity and `radius=1.5` corresponds
to 8-connectivity. `RegularGridAdjacency2D` is the explicit capability name;
there is no second compatibility type for the same relation.

The relation stores an immutable stencil. Traversal ranges own their own cursor,
so concurrent, interleaved, and nested iteration over the same relation is
safe and allocation-free:

```cpp
const RegularGridAdjacency2D adjacency(rows, cols, 1.5);
for (int p : adjacency.getNeighborIndices(properPart)) {
    for (int q : adjacency.getForwardNeighborIndices(p)) {
        // The inner traversal does not alter the outer traversal.
    }
}
```

Radius-based disks are not the only supported stencil. Centered rectangles,
digital lines, and arbitrary adjacency-inducing structuring elements use the
same concrete relation and therefore add no type erasure or virtual dispatch to
neighbor loops:

```cpp
auto rectangle = RegularGridAdjacency2D::rectangular(
    rows, cols, 1, 2);
auto horizontal = RegularGridAdjacency2D::horizontalLine(
    rows, cols, 3);

std::array<GridOffset2D, 5> cross{{
    {0, 0}, {-1, 0}, {0, -1}, {0, 1}, {1, 0}
}};
auto custom = RegularGridAdjacency2D::fromStructuringElement(
    rows, cols, cross);

auto maxTree = MorphologicalTreeFactory::createMaxTree(
    image, custom);
```

An adjacency-inducing structuring element must contain the origin exactly once
and be centrally symmetric. This preserves the undirected graph contract used
by `getForwardNeighborIndices()`. A component-tree adjacency must additionally
connect its finite image domain; for example, a horizontal-only line connects a
one-row image but not a multi-row image. Formulas defined specifically for
4-/8-connectivity, such as BitQuads, reject other stencils instead of silently
interpreting them as 8-connectivity.

Traversal names end in `Indices` because their values are row-major indices in
the regular grid; the adjacency is not restricted to an image-valued tree.

Algorithms declare their actual requirements. For example, `MAX_DIST` requires
a globally monotone altitude order and uniform adjacency. It does not need to
identify the builder that produced the hierarchy.

| Operation family | Required capabilities |
| --- | --- |
| topology traversal, ancestry, `AREA` | finite proper-part ownership only |
| `LEVEL`, `GRAY_HEIGHT`, `VOLUME` | altitude buffer |
| image reconstruction and attribute mapping | regular 2D proper-part domain |
| moments, bounding boxes, contours, BitQuads | regular 2D proper-part domain |
| `MAX_DIST` | regular 2D domain, uniform adjacency, monotone altitude order |
| directional contour/saliency operations | regular 2D domain and compatible uniform or directional adjacency |

## Construction

Use `MorphologicalTreeFactory` for public construction:

```cpp
#include <mmcfilters/trees/MorphologicalTreeFactory.hpp>
#include <mmcfilters/utils/Image.hpp>

using namespace mmcfilters;

auto image = ImageUInt8::create(4, 4, std::uint8_t{0});

auto maxTree = MorphologicalTreeFactory::createMaxTree(image, 1.5);
auto minTree = MorphologicalTreeFactory::createMinTree(image, 1.5);
auto residual = MorphologicalTreeFactory::createSelfDualResidualTree(
    image, 1.5);
auto saturatedResidual =
    MorphologicalTreeFactory::createSaturatedSelfDualResidualTree(
        image, NodeId{0}, 1.5);
auto tos = MorphologicalTreeFactory::createTreeOfShapes(image);
```

### Self-dual residual trees

`createSelfDualResidualTree` and
`createSaturatedSelfDualResidualTree` use one mutable max-tree and one mutable
min-tree built with the same symmetric adjacency. Current regional extrema are
selected by increasing support cardinality, with a deterministic common tie
policy. The unrestricted construction accepts every current extremum. The
saturated construction additionally requires its complement to remain
connected to the row-major `infinityPixel` supplied by the caller.

An explicit `RegularGridAdjacency2D` may be supplied instead of a radius. Its
stencil may be any centrally symmetric structuring element whose induced graph
connects the finite image domain:

```cpp
auto adjacency = RegularGridAdjacency2D::fromStructuringElement(
    rows, cols, offsets);
auto residual = MorphologicalTreeFactory::createSelfDualResidualTree(
    image, adjacency);
auto saturated =
    MorphologicalTreeFactory::createSaturatedSelfDualResidualTree(
        image, adjacency, infinityPixel);
```

Both results have descriptive kind `SELF_DUAL_RESIDUAL_TREE`, uniform adjacency,
and `AltitudeOrder::UNCONSTRAINED`. Their gray-level valuation alternates
polarities and is therefore not itself a hierarchy saliency valuation.
Topological levels, support cardinality, or another non-decreasing hierarchy
valuation can be projected with `HierarchySaliencyMap`; normalized-altitude
projection intentionally rejects these trees.

Direct users of `MinMaxResidualTreeBuilder` can inspect only the semantic
diagnostics retained for consistency and regression checks: emitted residual
events, rejected extrema, and exact complement-traversal certificates. The
production builder does not collect phase durations or detailed profiling
counters. Construction time must be measured externally around the complete
build call, as done by the benchmarks under `benchmarks/`.

### Tree-of-Shapes producer policies

Tree-of-Shapes construction variants belong to `TreeOfShapesProducer`, not to
`MorphologicalTree`. `TreeOfShapesProducerOptions` selects:

- `SelfDual`, `Min4cMax8c`, or `Min8cMax4c` interpolation/connectivity;
- an exterior immersion ring or no padding;
- the \f$p_\infty\f$ seed in the selected transient immersion domain.

The default remains the historical exterior-padded self-dual construction. A
convenience factory overload selects that producer policy directly:

```cpp
auto selfDual = MorphologicalTreeFactory::createTreeOfShapes(
    image,
    ToSInterpolation::SelfDual);

TreeOfShapesProducerOptions options{
    ToSInterpolation::Min4cMax8c,
    ToSPaddingPolicy::NoPadding,
    0,
    0};
auto unpadded = MorphologicalTreeFactory::createTreeOfShapes(
    image,
    options);
```

Padding changes only the transient 2D immersion. With exterior padding its
extent is `(2 * rows + 1, 2 * cols + 1)`; without padding it is
`(2 * rows - 1, 2 * cols - 1)`. Both factories publish exactly one proper part
per source pixel and attach the original `GridDomain2D`. Producer options are
not copied into the tree, so generic attribute computation sees only the
published topology, ownership, altitudes, domain, and adjacency capabilities.

The `uint8_t` factory quantizes self-dual half-level structural nodes. C++
callers that need those levels can request an exact owner:

```cpp
WeightedMorphologicalTree<ToSGrayLevel> exact =
    MorphologicalTreeFactory::createTreeOfShapesExact(
        image,
        options);
```

Exact altitudes use doubled gray units: source value `v` is represented by
`2 * v`, and odd values represent half levels. Consequently,
`reconstructionImage()` also returns doubled values. The topology, node ids,
proper-part ownership, and original 2D domain are otherwise the same as in the
quantized result. Python currently exposes the producer options and the
`uint8_t` factory; the typed exact owner remains C++-only.

Max-tree and min-tree construction is implemented by a concrete component-tree
producer. It preserves the union-find ordering while emitting generic
node-parent, proper-part-owner, and altitude buffers. The factory then uses the
same proof-carrying native materialization boundary as ToS; no
component-tree algorithm is embedded in `MorphologicalTree`. The union-find
implementation is a private producer detail; there is no virtual
“morphological-tree builder” base because the producers do not share a runtime
polymorphic contract.

C++ max-tree and min-tree factories are typed:

```cpp
ImagePtr<float> floatImage = Image<float>::create(4, 4, 0.0f);
WeightedMorphologicalTree<float> weighted =
    MorphologicalTreeFactory::createMaxTree(floatImage, 1.5);
```

Python construction is intentionally narrower:

- input images must be two-dimensional, C-contiguous `np.uint8` arrays;
- max-tree and min-tree factories return the canonical Python
  `WeightedMorphologicalTree` binding;
- Python does not currently expose typed `int32` or `float32` tree owners.

## Owner And View Boundary

Use `WeightedMorphologicalTree<T>` when the operation needs ownership:

- construction from images or Higra arrays;
- topology edits;
- altitude mutation;
- image reconstruction;
- Higra export/projection;
- filters, extinction values, UAO, CASF, and paired min/max adjustment;
- Python-facing workflows.

Use `WeightedTreeView<T>` for read-only C++ kernels that should not own or copy
altitude storage:

```cpp
const MorphologicalTree& topology = weighted.topology();
std::vector<float> altitude(topology.getNumInternalNodeSlots(), 0.0f);

WeightedTreeView<float> view(topology, altitude);
auto value = view.getAltitude(topology.getRoot());
```

The caller must keep the topology and altitude storage alive for the lifetime of
the view. The view captures the topology mutation version and rejects reads when
the topology changes.

## NodeId And Proper Parts

The canonical tree domain is the dense internal `NodeId` space:

```text
0 <= node_id < tree.getNumInternalNodeSlots()
```

Some internal slots may be dead after edits. Node-indexed buffers keep the full
slot count so existing row ids remain addressable. Iterate
`tree.getAliveNodeIds()` when consuming live tree nodes.

Image pixels are proper parts. They use row-major linear ids:

```text
proper_part = row * numCols + col
```

The key ownership queries are:

```cpp
NodeId owner = tree.getProperPartOwner(pixel);
auto direct = tree.getProperParts(node);
auto support = tree.getConnectedComponent(node);
```

`getProperParts(node)` returns only the direct proper parts owned by `node`.
`getConnectedComponent(node)` returns the full subtree support.

Internal nodes and proper parts are independent domains. A live node is not
required to own a direct proper part: its support may be contributed entirely
by descendants. This is particularly important for a branching root, but the
contract applies uniformly to every native morphological tree.
`isStructuralNode(node)` derives this role from empty direct ownership.

A committed live node must nevertheless have non-empty full subtree support.
An empty node may exist temporarily inside an edit session, but checked native
import, `validateAndCommit()`, and `commit()` reject it if it remains empty.
This distinguishes a valid structural node from an attached hierarchy artifact
that represents no partial partition.

Builders that already have this representation can import it without encoding
nodes as pixels:

```cpp
std::vector<NodeId> nodeParent{0, 0, 0};
std::vector<NodeId> properPartOwner{1, 2};
std::vector<std::uint8_t> altitude{0, 1, 0};

auto tree = MorphologicalTreeFactory::createFromNativeTopology(
    std::span<const NodeId>(nodeParent),
    std::span<const NodeId>(properPartOwner),
    std::span<const std::uint8_t>(altitude),
    0, // root
    HierarchySemantics{
        MorphologicalTreeKind::GENERIC,
        AltitudeOrder::UNCONSTRAINED,
        NoAdjacency{}});
```

Here the root owns no direct proper part and has two children. `AREA` remains
the number of proper parts in the complete subtree support. This overload
attaches no coordinate interpretation: `hasGridDomain2D()` is false, so
topological/support attributes work while reconstruction, contours, moments,
bounding boxes, and other geometry-dependent operations reject the missing
capability explicitly.

The overload with `rows` and `cols` attaches a `GridDomain2D`; it additionally
requires `properPartOwner.size() == rows * cols`. Tree-of-Shapes construction
uses the same generic representation and always attaches the source image
domain; its internal immersion grid is not materialized as the public
proper-part domain. No 3D layout is implied or implemented.

Native C++ builders can expose the same contract directly through
`NativeHierarchyView<T>`:

```cpp
auto weighted =
    MorphologicalTreeFactory::createFromNativeHierarchy(
        NativeHierarchyView<std::uint8_t>{
            nodeParent,
            properPartOwner,
            altitude,
            root,
            std::nullopt,
            HierarchySemantics{}});
```

The view is consumed synchronously and copied into owning tree storage.
Producer-owned validated buffers are instead moved into that storage. Max/min
and ToS use this boundary. Higra first converts
its leaf-first id space in the interoperability adapter and then crosses the
same native materialization boundary.

`getDescriptiveKind()` is the family-label accessor.

Optional adjacency semantics are also generic. `DirectionalGridAdjacency2D` stores
one relation for decreasing branches and another for increasing branches.
Algorithms that need this distinction query
`hasDirectionalGridAdjacency2D()`, `getDecreasingGridAdjacency2D()`, and
`getIncreasingGridAdjacency2D()` without testing the concrete tree kind.
They select the branch relation by comparing child and parent altitudes. If
the two connectivities differ, equal altitudes are rejected as ambiguous;
there is no hidden Tree-of-Shapes polarity sidecar in the generic topology.

## Altitude Contract

Altitude buffers are indexed by dense internal `NodeId`:

```text
altitude[node_id]
```

`WeightedMorphologicalTree<T>` owns an `AltitudeBuffer<T>`. `WeightedTreeView<T>`
borrows an `AltitudeSpan<T>`. Accepted altitude types satisfy `AltitudeValue`:

- floating-point values are accepted when finite;
- integral values must fit safely inside signed 64-bit differences;
- `bool` is rejected.

Weighted construction and edits validate the declared `AltitudeOrder` with
strict inequalities on every alive parent-child arc.
`UNCONSTRAINED` hierarchies skip only the monotonicity check; they still validate
finite altitude values and buffer shape. Operators declare their real
requirements: altitude-neighborhood stability and `MAX_DIST`, for example,
require a globally monotone order, independently of `MorphologicalTreeKind`.
Use topological depth-stability operators when `delta` should mean a number of
tree edges instead of an altitude difference.

`GRAY_HEIGHT` is defined for every altitude order as the maximum absolute
altitude difference between a node and any node in its subtree. On monotone
max/min hierarchies this is exactly the traditional one-sided grey-level span.

## Hierarchy Edge Saliency

The mathematical definitions, operator comparison, validity contracts, tie
policy, and complexity bounds are centralized in
[Hierarchy Saliency Maps](saliency.md). This section records the tree-facing
projection API. For an executable English tutorial, see
[`Saliency_Maps_Tutorial.ipynb`](../notebooks/Saliency_Maps_Tutorial.ipynb).
The complete bibliographic references and the code-to-paper correspondence are
recorded in the canonical guide's
[primary-reference section](saliency.md#primary-references-and-implementation-correspondence).

`HierarchySaliencyMap` projects a tree hierarchy onto an image adjacency graph.
The formal API is valuation-based: a caller supplies one scalar valuation per
internal node, and that valuation must be compatible with the hierarchy:

```text
valuation(parent) >= valuation(child)
```

For every undirected adjacency edge `(p, q)`, the implementation finds the
direct proper-part owners and computes their lowest common ancestor. Two level
conventions are explicit. With `HierarchyLevelConvention::EdgeSaliencyValue`,
the node valuation is already the value written on a transition edge:

```text
saliency(p, q) = valuation(LCA(owner(p), owner(q)))  if owner(p) != owner(q)
saliency(p, q) = 0                                  otherwise
```

With `HierarchyLevelConvention::PartitionAppearanceLevel`, the valuation is the
positive index of the first partition in which the region appears and Algorithm
1 of Cousty is applied literally:

```text
saliency(p, q) = level(LCA(owner(p), owner(q))) - 1
```

The zero case keeps the map defined on the full graph edge set while respecting
the cut definition: an edge whose endpoints already belong to the same finest
region represented by the tree is not a transition contour at any positive
hierarchy level.

This implements the edge-indexed saliency map `Phi(H)` of a connected hierarchy
as defined by Cousty et al. for quasi-flat zones. It does not implement the full
paper pipeline `Psi(w) = Phi(QFZ(G, w))` from an arbitrary edge-weighted graph;
callers provide the hierarchy and a compatible node valuation.

The formal projection always enforces the paper's non-negative edge-weight
domain because edges internal to a finest represented region use the fixed base
level `0`. To convert an arbitrary compatible finite valuation to dense integer
levels `0..k-1`, use
`HierarchySaliencyMapValidation::rankHierarchyValuation(...)` before projection
when node ranks themselves are required. For the canonical dense edge scale, use
`computeCanonicalRankedSaliencyEdgeMap(...)`: it ranks only levels that actually
occur on graph edges, so unused leaf values cannot introduce gaps.

The output is an `EdgeSaliencyMap<T>` with parallel `sources`, `targets`, and
`values` arrays. This keeps the result edge-indexed, matching the saliency-map
representation used for connected hierarchies and quasi-flat zones. Convert it
to a display image only as a later visualization step.

```cpp
#include <mmcfilters/trees/saliency/HierarchySaliencyMapValidation.hpp>
#include <mmcfilters/trees/saliency/HierarchySaliencyMap.hpp>

std::vector<float> valuation(tree.getNumInternalNodeSlots(), 0.0f);
HierarchySaliencyMapValidation::validateHierarchyValuation(
    tree,
    std::span<const float>(valuation),
    HierarchyValuationPolicy::AllowLevelCollapse,
    HierarchyValuationRangePolicy::RequireNonNegative);

std::vector<int> ranked = HierarchySaliencyMapValidation::rankHierarchyValuation(
    tree,
    std::span<const float>(valuation));

std::vector<double> normalized = HierarchySaliencyMapValidation::computeNormalizedScores(
    tree,
    std::span<const float>(valuation));

auto edgeMap = HierarchySaliencyMap::computeSaliencyEdgeMap(
    tree,
    adjacency,
    std::span<const int>(ranked),
    HierarchyValuationPolicy::AllowLevelCollapse);

auto canonicalEdgeRanks = HierarchySaliencyMap::computeCanonicalRankedSaliencyEdgeMap(
    tree,
    adjacency,
    std::span<const float>(valuation));

std::vector<int> appearance =
    ComponentTreePartitionHierarchyAdapter::computePartitionAppearanceLevels(tree);
auto appearanceMap = ComponentTreePartitionHierarchyAdapter::computeSaliencyEdgeMap(
    tree,
    adjacency,
    std::span<const int>(appearance));

auto levelMap = HierarchySaliencyMap::computeTopologicalLevelEdgeMap(tree);
auto normalizedMap = HierarchySaliencyMap::computeNormalizedAltitudeEdgeMap(weightedTree);
```

`computeTopologicalLevelEdgeMap(...)` starts from level `0` on internal leaves
and one plus the maximum child level on ancestors, then canonically ranks the
effective graph-edge levels. This removes structural levels that do not change
any graph partition. It works for max-trees, min-trees, trees of shapes, and
imported hierarchies whose supports are connected in the selected graph.

`computeNormalizedAltitudeEdgeMap(...)` returns `double` values in `[0, 1]` and
orients max-tree/min-tree altitudes so ancestors receive values greater than or
equal to descendants. Max-tree altitudes are inverted; min-tree altitudes keep
their coarse-to-fine order. Equal altitudes remain equal, so this policy can
collapse distinct tree nodes that live at the same image level.

`HierarchySaliencyMapValidation::computeNormalizedScores(...)` is the generic
valuation normalizer. It validates a caller-supplied compatible valuation and
maps its live-node value range to `[0, 1]` with an increasing affine transform.
It does not infer altitude polarity; callers must supply a valuation that
already grows toward the root.

`computeEdgeMap(...)` remains the low-level LCA projection primitive for internal
experiments, but it does not validate monotonicity. Use `computeSaliencyEdgeMap`
when the result must be a formal quasi-flat-zone saliency map.

Formal projection validates spatial connectedness by default. This check is
distinct from validating one connected parent tree: every node support must be
connected in the graph `G`. It uses one LCA assignment per graph edge followed
by a post-order disjoint-set pass, for `O(m + p + e)` work. Trusted internal
paths may explicitly select `HierarchyConnectivityPolicy::AssumeConnected`.

`ComponentTreePartitionHierarchyAdapter` names the completion required because
a component tree is originally a hierarchy of partial partitions. Partition 0
contains pixel singletons; zero-valued same-owner edges form connected direct
proper-part atoms; component-tree nodes merge those atoms and child supports.
The adapter validates that this completion is connected in `G` and supplies
positive partition-appearance levels for the `level(LCA)-1` convention.

When no explicit adjacency is supplied, the tree must carry the construction
adjacency. A uniform relation is used directly. A directional context can be
used implicitly only when its decreasing and increasing relations have the same
stencil and therefore define one unambiguous projection graph. When those
relations differ—for example, a 4/8 Tree of Shapes—the formal saliency map still
needs one fixed graph `G`, so callers must pass the intended
`RegularGridAdjacency2D` explicitly. Imported or self-dual trees without an
attached adjacency have the same explicit requirement. Python entrypoints that
construct this adjacency from a radius require a finite radius of at least
`1.0`.

## Hierarchy Saliency Projections

`HierarchySaliencyMapProjection` materializes contours from the edge saliency
representation. A threshold cut selects the image-domain adjacency edges whose
saliency is greater than or equal to the requested level:

```text
contour_lambda = { (p, q) | saliency(p, q) >= lambda }
```

This is the cut complement of the quasi-flat-zone convention in which edges
with saliency strictly below the level remain connected.

```cpp
#include <mmcfilters/trees/saliency/HierarchySaliencyMapProjection.hpp>

auto saliency = HierarchySaliencyMap::computeNormalizedAltitudeEdgeMap(weightedTree);
auto contour = HierarchySaliencyMapProjection::thresholdCut(saliency, 0.5);
```

For visualization or notebook inspection, an edge map can be rasterized to a
pixel image by aggregating incident edge values at each endpoint. This helper is
not the formal saliency map; it is only a display projection:

```cpp
auto pixelView = HierarchySaliencyMapProjection::edgeMapToPixelImage(
    saliency,
    EdgeToPixelReducer::Max);
```

The recommended sparse contour representation is the per-node incremental
contour map. For each node `u`, it stores the transition-edge set

```text
E_u = { (p, q) | owner(p) != owner(q) and LCA(owner(p), owner(q)) = u }
```

as a contiguous slice:

```cpp
auto incremental = HierarchySaliencyMapProjection::computeIncrementalNodeContours(weightedTree);

for (NodeId u = 0; u < incremental.numNodeSlots; ++u) {
    for (std::size_t i = incremental.offsets[u]; i < incremental.offsets[u + 1]; ++i) {
        NodeId p = incremental.sources[i];
        NodeId q = incremental.targets[i];
        // (p, q) is a transition contour edge born at node u.
    }
}
```

Dense node-indexed valuations can then be projected onto that sparse transition
support without recomputing LCAs:

```cpp
std::vector<float> nodeValuation(tree.getNumInternalNodeSlots(), 0.0f);
auto edgeMap = HierarchySaliencyMapProjection::projectNodeValuation(
    incremental,
    std::span<const float>(nodeValuation));
auto cut = HierarchySaliencyMapProjection::thresholdByNodeValuation(
    incremental,
    std::span<const float>(nodeValuation),
    0.5f);
```

These incremental helpers are generic projection/cut routines over transition
edges. They do not materialize the full graph edge set; edges internal to finest
regions are omitted and have implicit value `0` in the formal saliency map. Use
`HierarchySaliencyMap::computeSaliencyEdgeMap` when every graph edge must be
present explicitly.

For debugging or simple inspection, `nodeContourEdges(...)` returns a flat
`NodeContourEdgeMap` with parallel `sources`, `targets`, and `nodes` arrays.
The incremental map is the preferred representation when the goal is to treat
node valuations as a hierarchy valuation and materialize saliency cuts.

## Shape-Space Extinction Saliency

`ShapeSpaceSaliency` is a separate construction for Xu-style shaping. It does
not interpret the input attribute as a monotone valuation of the original
hierarchy. Instead, it changes objects: every live node of the original tree is
a vertex of a new graph, and every original parent-child relation is an edge of
that graph. Regional minima or maxima of an arbitrary floating-point node
attribute are then ranked by their extinction in this shape space.

```cpp
#include <mmcfilters/trees/saliency/ShapeSpaceSaliency.hpp>

std::vector<double> attribute(tree.getNumInternalNodeSlots(), 0.0);
auto result = ShapeSpaceSaliency::compute(
    tree,
    std::span<const double>(attribute),
    ShapeSpaceExtremaPolarity::Maxima,
    adjacency);
```

The result contains the regional extrema, a dense node buffer with each
extinction placed on its canonical original-tree representative, and an
`EdgeSaliencyMap<double>`. Equal-valued connected nodes form one plateau. Its
representative is the outermost node in that plateau: the node whose parent is
not in the same plateau. Birth and death levels remain in the input attribute's
coordinate system, and the dominant extremum receives a finite death level at
the opposite global level.

The edge projection is maximum-on-contours. If `score(n)` is the sparse
extinction score of original region `n`, then

```text
saliency(p, q) = max { score(n) | (p, q) belongs to the boundary of region n }
```

Equivalently, let `u=owner(p)`, `v=owner(q)`, and `a=LCA(u,v)`. The relevant
regions are the nodes on the two paths from `u` and `v` toward `a`, excluding
`a` itself. Consequently, this operation is not
`score(LCA(owner(p), owner(q)))` and does not require `score` to be monotone
toward the original root.

The stages can also be called independently:

```cpp
auto extinction = ShapeSpaceSaliency::computeExtinctionValues(
    tree,
    std::span<const double>(attribute),
    ShapeSpaceExtremaPolarity::Maxima);

auto edgeMap = ShapeSpaceSaliency::projectContourScores(
    tree,
    std::span<const double>(extinction.nodeScores),
    adjacency);
```

Use `HierarchySaliencyMap` when a known hierarchy valuation must be projected
by LCA. Use `ShapeSpaceSaliency` when an arbitrary attribute on the original
regions must first induce a second component hierarchy and be converted to
extinction-weighted contours. The latter follows Xu, Carlinet, Géraud, and
Najman, *Hierarchical Segmentation Using Tree-Based Shape Spaces*, IEEE TPAMI
39(3), 457–469, 2017,
[DOI 10.1109/TPAMI.2016.2554550](https://doi.org/10.1109/TPAMI.2016.2554550).

## Editing And Derived State

Safe committed edits are public:

```cpp
tree.pruneNode(node);
tree.mergeNodeIntoParent(node);
```

Multi-step topology changes must go through an editor session:

```cpp
auto editor = tree.edit();
NodeId inserted = editor.createDetachedNode();
editor.reparent(childA, inserted);
editor.reparent(childB, inserted);
editor.attach(parent, inserted);
editor.commit();
```

Weighted edits use `WeightedMorphologicalTree<T>::edit()` and additionally
validate the declared altitude order on commit. Checked commits also reject any
live node whose complete subtree support is empty; detached or empty nodes are
therefore transient edit state only.

Every public editor is recoverable. A lazy copy-on-first-write journal records
only the topology, ownership, free-slot, and altitude entries affected by the
edit. Explicit `rollback()`, destruction, and exception unwinding restore the
previous committed owner without a complete-tree snapshot.

Topology mutations advance the tree mutation version. Objects that cache
node-indexed state against a tree reject reads after the topology changes. This
includes:

- `ContoursComputedIncrementally::IncrementalContours`;
- `AttributeFilters`;
- `ExtinctionValues`;
- `UltimateAttributeOpening`.

Plain attribute result buffers are not versioned; recompute them after topology
mutation.

## Higra Interoperability

`createFromHigraParent(...)` imports a compact Higra-style hierarchy whose
external layout is:

```text
[proper parts | internal nodes]
```

The imported Higra node-id domain is preserved only until topology mutation.
`NodeIdSpace::HIGRA` refers to that preserved imported domain, while
`exportHigraHierarchy()` creates a fresh compact snapshot for the current live
tree. For import/export rules and attribute projection to both domains, see
[Higra Interoperability](higra-interoperability.md).

## Python Notes

Python exposes the stable owner-oriented surface:

```python
max_tree = mmcfilters.MorphologicalTreeFactory.createMaxTree(image, radius=1.5)

root = max_tree.getRoot()
children = max_tree.getChildren(root)
pixels = max_tree.getProperParts(root)
component = max_tree.getConnectedComponent(root)
reconstructed = max_tree.reconstructionImage()
```

Python does not expose `TreeEditor`, `WeightedTreeEditor`, or a mutable topology
handle from `WeightedMorphologicalTree`. Unchecked altitude setters are not
part of either the C++ or Python public API.

Factory and import paths return `WeightedMorphologicalTree`. Python does not
register a separate topology-only class or topology-view return type, avoiding
a second object with ambiguous ownership and lifetime. Topology-only algorithms
read the topology owned by the weighted tree.

## Related Guides

- [Attributes](attributes.md): choosing weighted versus topology-only attribute
  computation.
- [Attribute Filters, Extinction Values, And UAO](filters.md): reconstructing
  image outputs from dense node-indexed buffers.
- [Higra Interoperability](higra-interoperability.md): preserved imported ids,
  exported compact layouts, and attribute projection.
- [Editing API](editing-api.md): safe mutators, staged edits, and derived-state
  lifetime.
- [Python API Guide](python-api.md): Python construction, `NodeId` domains,
  attributes, filters, contours, and Higra interop.
- [Incremental Contours](contours.md): pixel-contour extraction over tree
  supports.
- [Contour Traces](contour-traces.md): definitive geometric contour API with
  oriented boundary edges and external/internal loops.
