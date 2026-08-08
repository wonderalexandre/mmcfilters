# Python API Guide

This document summarizes the stable Python surface of `mmcfilters`. It focuses
on practical workflows and the contracts that are easiest to misuse from Python.
For the lower-level pybind11 declarations, use the generated internal
documentation and Python `help(...)` output.

For the C++/Python tree model shared by these bindings, see
[Morphological Trees](trees.md).
For the filtering operators and their C++/Python buffer contracts, see
[Attribute Filters, Extinction Values, And UAO](filters.md).

## Input Contract

Python construction is intentionally limited to the canonical 8-bit contract:

- image inputs must be 2D, C-contiguous `np.uint8` arrays;
- altitude buffers passed from Python must be 1D `np.uint8` arrays or integer
  sequences in `[0, 255]`;
- max-tree and min-tree factories return `WeightedMorphologicalTree`;
- a separate topology-only `MorphologicalTree` Python type is not exposed.

Use `np.ascontiguousarray` at application boundaries when data may come from a
view, transpose, or external image loader:

```python
import numpy as np
import mmcfilters

image = np.asarray(
    [
        [3, 3, 2, 2],
        [3, 4, 4, 2],
        [1, 4, 5, 2],
        [1, 1, 5, 0],
    ],
    dtype=np.uint8,
)
image = np.ascontiguousarray(image)
```

## Building Trees

Component-tree factories use an adjacency radius. On the 2D square grid,
`radius=1.0` corresponds to 4-connectivity and `radius=1.5` corresponds to
8-connectivity.

```python
max_tree = mmcfilters.MorphologicalTreeFactory.createMaxTree(image, radius=1.5)
min_tree = mmcfilters.MorphologicalTreeFactory.createMinTree(image, radius=1.5)
residual = mmcfilters.MorphologicalTreeFactory.createSelfDualResidualTree(
    image,
    radius=1.5,
)
saturated_residual = (
    mmcfilters.MorphologicalTreeFactory.createSaturatedSelfDualResidualTree(
        image,
        infinityPixel=0,
        radius=1.5,
    )
)
tos = mmcfilters.MorphologicalTreeFactory.createTreeOfShapes(image)
```

The altitude-order capability uses strict parent-child contracts:

- `AltitudeOrder.INCREASING_FROM_ROOT` means
  `altitude(parent) < altitude(child)`;
- `AltitudeOrder.DECREASING_FROM_ROOT` means
  `altitude(parent) > altitude(child)`;
- `AltitudeOrder.UNCONSTRAINED` declares no global direction.

Max-trees and min-trees receive the first and second contracts respectively.
Trees of Shapes and self-dual residual trees remain `UNCONSTRAINED`; the generic
tree abstraction does not infer a single altitude polarity. Residual trees can
still be projected to formal saliency maps with topological levels or another
non-decreasing hierarchy valuation. Their alternating gray-level altitudes are
intentionally rejected by normalized-altitude projection.

The two residual-tree factories use synchronized max-tree/min-tree states with
one shared symmetric adjacency. The saturated form restricts selection to
regional extrema whose complement remains connected to `infinityPixel`.
`SdrtTiePolicy.CONTRAST_INVARIANT_SPATIAL` is the default deterministic tie
policy for both factories.

A generic hierarchy over an abstract finite proper-part set does not need
`rows` or `cols`:

```python
semantics = mmcfilters.HierarchySemantics()
abstract_tree = (
    mmcfilters.MorphologicalTreeFactory.createFromNativeTopology(
        [0, 0, 0],       # node parent
        [1, 2],          # direct proper-part owners
        np.array([10, 3, 20], dtype=np.uint8),
        0,               # root
        semantics=semantics,
    )
)
assert not abstract_tree.hasGridDomain2D
```

Pass `rows`, `cols`, and the same `semantics` argument to attach a regular 2D
layout. This enables reconstruction and geometric attributes; it does not
change topology or proper-part ownership. No 3D domain is implied.

Directional 2D semantics use the explicit
`DirectionalGridAdjacency2D(decreasing, increasing)` type. The shorter
`DirectionalAdjacency` compatibility name is not exposed.

Attribute requirements are inspectable before computation:

```python
requirements = mmcfilters.Attribute.requirements(
    mmcfilters.Attribute.MAX_DIST
)
assert requirements["gridDomain2D"]
assert requirements["adjacency"] == "uniform"
assert requirements["monotoneAltitudeOrder"]
assert not requirements["canonical4Or8Adjacency"]
```

The returned object owns both topology and a dense altitude buffer:

```python
root = max_tree.getRoot()
alive_nodes = max_tree.getAliveNodeIds()
leaf_nodes = max_tree.getLeafNodeIds()
altitude_at_root = max_tree.getAltitude(root)
reconstructed = max_tree.reconstructionImage()
```

All factory and import paths keep returning `WeightedMorphologicalTree` in the
current API. A separate read-only topology-view object is not exposed:
it would duplicate the query surface while introducing lifetime and identity
questions without enabling a current public workflow. Topology/support-only
attribute and contour operations accept the weighted tree and use its owned
topology directly.

## NodeId Domains

The main Python API uses the dense internal `NodeId` domain of the current
`MorphologicalTree`. Attribute rows, altitude buffers, predicates, traversal
results, and mutation methods are indexed by this same internal domain.

After topology edits, some internal node slots may be dead. Dense outputs keep
the full `tree.numInternalNodeSlots` shape so existing row ids remain
addressable; iterate `tree.getAliveNodeIds()` when a consumer needs only live
nodes.

```python
for node_id in max_tree.getAliveNodeIds():
    children = max_tree.getChildren(node_id)
    direct_pixels = max_tree.getProperParts(node_id)
    subtree = max_tree.getNodeSubtree(node_id)
```

Image pixels are proper parts. They use row-major linear ids:

```python
pixel_id = 10
owner_node = max_tree.getProperPartOwner(pixel_id)
component_pixels = list(max_tree.getConnectedComponent(owner_node))
component_mask = max_tree.reconstructNode(owner_node)
```

## Attributes

Use weighted entry points for attributes that may read altitude values, and
topology entry points for support-only descriptors.

Single attributes return a 1D floating-point array indexed by the selected
`NodeIdSpace`. Multi-attribute calls return `(layout, values)`, where `layout`
maps attribute names to columns in a 2D floating-point array. The optional
`dtype` keyword accepts `np.float32` or `np.float64`; the default is
`np.float32`.

The selected `dtype` controls the NumPy output array. Attribute values are
computed internally through the same double-precision facade and cast only when
the output buffer is materialized.

```python
area = mmcfilters.Attribute.computeSingleTopologyAttribute(
    max_tree,
    mmcfilters.Attribute.AREA,
)

level = mmcfilters.Attribute.computeSingleAttribute(
    max_tree,
    mmcfilters.Attribute.LEVEL,
)

layout, values = mmcfilters.Attribute.computeAttributes(
    max_tree,
    [
        mmcfilters.Attribute.AREA,
        mmcfilters.Attribute.VOLUME,
        mmcfilters.Attribute.RELATIVE_VOLUME,
    ],
    dtype=np.float64,
)

area_column = values[:, layout["AREA"]]
volume_column = values[:, layout["VOLUME"]]
```

Filtering helpers documented below consume `np.float32` and `np.float64`
attribute buffers, so arrays can be passed through without downcasting after
selecting either supported dtype.

Attribute groups expand to a stable set of scalar attributes:

```python
gray_layout, gray_values = mmcfilters.Attribute.computeAttributes(
    max_tree,
    [mmcfilters.Attribute.Group.GRAY_LEVEL],
)

boundary_layout, boundary_values = mmcfilters.Attribute.computeTopologyAttributes(
    max_tree,
    [mmcfilters.Attribute.Group.BOUNDARY],
)
```

Delta-augmented attributes sample one scalar attribute along ancestor and
descendant offsets:

```python
delta_layout, delta_values = mmcfilters.Attribute.computeSingleAttributeWithDelta(
    max_tree,
    mmcfilters.Attribute.AREA,
    1,
    "null-padding",
)

center = delta_values[:, delta_layout["AREA"]]
asc_1 = delta_values[:, delta_layout["AREA_ASC_1"]]
desc_1 = delta_values[:, delta_layout["AREA_DESC_1"]]
```

The `padding` argument controls missing ancestor/descendant samples:
`"last-padding"` repeats the nearest available value, `"zero-padding"` writes
zero, and `"nan-padding"` or `"null-padding"` leave missing samples as `NaN`.

To project one node attribute back to the image domain, use
`computeAttributeMapping`:

```python
area_by_pixel = mmcfilters.Attribute.computeAttributeMapping(
    max_tree,
    mmcfilters.Attribute.AREA,
)
```

## Hierarchy Edge Saliency

The canonical mathematical definitions, preconditions, tie policy, and
complexity bounds are in [Hierarchy Saliency Maps](saliency.md). Its
[primary-reference section](saliency.md#primary-references-and-implementation-correspondence)
provides complete citations and identifies the exact equations, algorithms,
sections, and library-specific generalizations. This section lists the Python
signatures and return containers.

Use `HierarchySaliencyMap` when the desired saliency map is the edge-indexed
projection of a hierarchy onto an adjacency graph. The formal API is
valuation-based: the valuation must be a dense node-indexed `int32`, `float32`,
or `float64` array that is non-decreasing toward the root. Under
`HierarchyLevelConvention.EdgeSaliencyValue`, each transition edge receives the
valuation of the lowest common ancestor of its endpoint owners. Under
`PartitionAppearanceLevel`, it receives `level[LCA] - 1`, following Algorithm 1
of Cousty. Same-owner edges receive the base value `0`.

This API implements the edge-indexed saliency map `Phi(H)` of a connected
hierarchy as defined by Cousty et al. for quasi-flat zones. It does not implement
the full paper pipeline `Psi(w) = Phi(QFZ(G, w))` from an arbitrary edge-weighted
graph; callers provide the hierarchy and a compatible node valuation.

Formal projection requires a non-negative valuation because same-owner edges
receive the explicit base value `0`. Use
`HierarchySaliencyMapValidation.rankHierarchyValuation(...)` to convert an
arbitrary finite ordered valuation to dense integer levels `0..k-1` before
projection. This ranks every node. Use
`HierarchySaliencyMap.computeCanonicalRankedSaliencyEdgeMap(...)` when ranks
must be dense only over effective graph transitions.

```python
valuation = np.zeros(max_tree.numInternalNodeSlots, dtype=np.float32)
ranked = mmcfilters.HierarchySaliencyMapValidation.rankHierarchyValuation(
    max_tree,
    valuation,
)
normalized = mmcfilters.HierarchySaliencyMapValidation.computeNormalizedScores(
    max_tree,
    valuation,
)
edge_map = mmcfilters.HierarchySaliencyMap.computeSaliencyEdgeMap(
    max_tree,
    ranked,
)

canonical_edge_map = mmcfilters.HierarchySaliencyMap.computeCanonicalRankedSaliencyEdgeMap(
    max_tree,
    valuation,
)

sources = edge_map["sources"]
targets = edge_map["targets"]
values = edge_map["values"]
```

The returned dictionary also contains `numRows`, `numCols`, and
`adjacencyRadius`.

Validate a valuation directly when the caller needs to fail before projection.
Valuations may be `int32`, `float32`, or `float64`:

```python
mmcfilters.HierarchySaliencyMapValidation.validateHierarchyValuation(
    max_tree,
    valuation,
    strict=False,
    nonnegative=False,
)
```

With `strict=False`, parent values may equal child values, which can collapse
adjacent tree levels in the induced quasi-flat-zone hierarchy. With
`strict=True`, every live parent-child relation must be strictly increasing
toward the root. With `nonnegative=True`, every live-node valuation must be
greater than or equal to zero.

Formal projection validates by default that every hierarchy support is
connected in the selected graph. This is different from checking that the
parent array is one rooted tree. It can also be invoked directly:

```python
mmcfilters.HierarchySaliencyMapValidation.validateHierarchyConnectivity(max_tree)
```

For a component tree, the partial-partition completion is explicit through
`ComponentTreePartitionHierarchyAdapter`: pixel singletons form partition 0,
same-owner zero edges form proper-part atoms, and component nodes merge those
atoms and child supports.

```python
appearance = mmcfilters.ComponentTreePartitionHierarchyAdapter.computePartitionAppearanceLevels(max_tree)
appearance_map = mmcfilters.ComponentTreePartitionHierarchyAdapter.computeSaliencyEdgeMap(
    max_tree,
    appearance,
    strict=True,
)
```

If the tree does not carry a construction adjacency, pass a radius explicitly.
The same applies when a directional Tree of Shapes carries distinct decreasing
and increasing stencils: a formal saliency map is defined on one fixed graph,
so the caller must choose it. A directional context whose two stencils coincide
is unambiguous and can be used without `radius`. An explicit radius must be
finite and at least `1.0`, ensuring that the stencil contains at least one
non-central neighbour.

```python
edge_map = mmcfilters.HierarchySaliencyMap.computeSaliencyEdgeMap(
    imported_tree,
    valuation,
    radius=1.5,
)
```

Convenience policies cover common hierarchy scales:

```python
level_map = mmcfilters.HierarchySaliencyMap.computeTopologicalLevelEdgeMap(
    max_tree,
)
normalized_map = mmcfilters.HierarchySaliencyMap.computeNormalizedAltitudeEdgeMap(
    max_tree,
)
```

Topological levels are structural: internal leaves receive `0`, and ancestors
increase toward the root. Normalized altitude maps return `float64` values in
`[0, 1]` for max-trees and min-trees, with max-tree altitudes inverted so that
ancestors are greater than or equal to descendants. Trees without a single
max/min altitude polarity are rejected by the normalized-altitude policy.
Use `HierarchySaliencyMapValidation.computeNormalizedScores(...)` when the input
is already a compatible generic valuation and only needs to be reparameterized
to `[0, 1]`.

## Shape-Space Extinction Saliency (Xu Shaping)

`ShapeSpaceSaliency` implements the shaping construction used by Xu et al. and
is deliberately separate from `HierarchySaliencyMap`. The input morphological
tree is first treated as a graph whose vertices are the original tree nodes and
whose edges are the original parent-child relations. A second component tree is
formed from an arbitrary dense node score on that graph. Local minima or maxima
of this second-level hierarchy are selected with
`ShapeSpaceExtremaPolarity.Minima` or `ShapeSpaceExtremaPolarity.Maxima`, and
their extinction values are measured between their birth and death levels.

```python
scores = np.asarray(node_attribute, dtype=np.float32)
extinction = mmcfilters.ShapeSpaceSaliency.computeExtinctionValues(
    tree,
    scores,
    mmcfilters.ShapeSpaceExtremaPolarity.Minima,
)

for extremum in extinction["extrema"]:
    representative = extremum["representative"]
    value = extremum["extinction"]

dense_extinction = extinction["nodeScores"]
```

The result preserves the input floating-point dtype. `extrema` contains
`representative`, `birthLevel`, `deathLevel`, and `extinction` for every
selected extremum; `nodeScores` is the dense node-indexed extinction score used
for contour projection.

To materialize the image-domain map, each original region contour receives its
node score, and every adjacency edge keeps the maximum score among the region
contours that contain it:

```python
edge_map = mmcfilters.ShapeSpaceSaliency.projectContourScores(
    tree,
    dense_extinction,
)

result = mmcfilters.ShapeSpaceSaliency.compute(
    tree,
    scores,
    mmcfilters.ShapeSpaceExtremaPolarity.Minima,
)
edge_map = result["edgeMap"]
```

The `compute` result contains the same `extrema` and `nodeScores` entries plus
the usual edge-map dictionary under `edgeMap`. The Shape Space convenience
overload requires a stored uniform construction adjacency; pass `radius=...`
when that uniform relation is unavailable.

This is not the `valuation[LCA(owner(source), owner(target))]` projection of
`HierarchySaliencyMap`. In Xu shaping, extinction is computed on the
second-level node graph, then combined on the contours of the original regions
with a maximum. The distinction matters when several original regions share an
image edge or when the arbitrary input score is not a monotone hierarchy
valuation.

References: Xu, Carlinet, Géraud, and Najman,
[*Hierarchical Segmentation Using Tree-Based Shape Spaces*](https://doi.org/10.1109/TPAMI.2016.2554550),
IEEE TPAMI 39(3), 457–469, 2017; and Xu, Géraud, and Najman,
[*Connected Filtering on Tree-Based Shape-Spaces*](https://doi.org/10.1109/TPAMI.2015.2441070),
IEEE TPAMI 38(6), 1126–1140, 2016. The representative-node and maximum-on-contours
projection follows Section 4.3 of the first article. The
[Higra shaping notebook](https://github.com/higra/Higra-Notebooks/blob/master/Computing%20a%20saliency%20map%20with%20the%20shaping%20framework.ipynb)
is an interoperability example, not the normative definition of the algorithm.

## Hierarchy Saliency Projections

Use `HierarchySaliencyMapProjection` when a saliency map should be turned into
explicit contour edges. A threshold cut keeps every edge whose saliency is
greater than or equal to the threshold:

```python
normalized_map = mmcfilters.HierarchySaliencyMap.computeNormalizedAltitudeEdgeMap(
    max_tree,
)
contour = mmcfilters.HierarchySaliencyMapProjection.thresholdCut(
    normalized_map,
    threshold=0.5,
)
```

For display only, an edge-indexed map can be rasterized to a pixel image. This
is not the formal saliency representation; each edge value is only
aggregated into its two endpoint pixels:

```python
pixel_view = mmcfilters.HierarchySaliencyMapProjection.edgeMapToPixelImage(
    normalized_map,
    reducer=mmcfilters.EdgeToPixelReducer.Max,
)
```

The threshold-cut dictionary contains `numRows`, `numCols`, `adjacencyRadius`,
`sources`, and `targets`. To attribute transition edges back to hierarchy nodes
in a flat debug layout, use the node projection:

```python
node_edges = mmcfilters.HierarchySaliencyMapProjection.nodeContourEdges(max_tree)
nodes = node_edges["nodes"]
```

Here `nodes[i]` is `LCA(owner(sources[i]), owner(targets[i]))`, so it identifies
the hierarchy node whose boundary contains that transition edge. Edges whose
endpoints have the same owner are omitted; in the full formal saliency map they
have implicit value `0`.

For hierarchy-valued experiments, prefer the incremental per-node layout:

```python
contours = mmcfilters.HierarchySaliencyMapProjection.computeIncrementalNodeContours(
    max_tree,
)
```

The returned dictionary has `offsets`, `sources`, and `targets`. For node `u`,
the slice `offsets[u]:offsets[u + 1]` contains the transition contour edges born
at that node. Dense node valuations can be projected to a sparse transition-edge map,
or thresholded directly:

```python
node_valuation = np.zeros(max_tree.numInternalNodeSlots, dtype=np.float32)
edge_map = mmcfilters.HierarchySaliencyMapProjection.projectNodeValuation(
    contours,
    node_valuation,
)
cut = mmcfilters.HierarchySaliencyMapProjection.thresholdByNodeValuation(
    contours,
    node_valuation,
    threshold=0.5,
)
```

The incremental projection helpers do not materialize internal zero-valued
edges. Use `HierarchySaliencyMap.computeSaliencyEdgeMap(...)` when every graph
edge must be present explicitly. Validate or rank `node_valuation` through
`HierarchySaliencyMapValidation` first when the transition-edge values must
satisfy the formal saliency-map contract.

## Filtering

Filtering helpers snapshot the tree topology at construction time. Create
`AttributeFilters`, `ExtinctionValues`, or `UltimateAttributeOpening` after the
structural mutations that should be reflected by the operation.
The same dense internal `NodeId` buffer contract is covered in more detail in
[Attribute Filters, Extinction Values, And UAO](filters.md).

```python
filters = mmcfilters.AttributeFilters(max_tree)

box_height = mmcfilters.Attribute.computeSingleAttribute(
    max_tree,
    mmcfilters.Attribute.BOX_HEIGHT,
)

pruned_min = filters.filteringByPruningMin(box_height, 2.0)
pruned_max = filters.filteringByPruningMax(box_height, 2.0)

keep_large = (area >= 4.0).tolist()
direct = filters.filteringDirectRule(keep_large)
subtractive = filters.filteringSubtractiveRule(keep_large)
```

Topological depth stability is available directly when you need the variation
scores or the selected mask:

```python
depth = mmcfilters.DepthStableRegionComputer(tree)
depth_mask = depth.computeByDepth(depthDelta=2)
depth_variation = depth.getVariations()
```

The result getters (`getVariation`, `getVariations`, `getNumNodes`, and the
window-neighbour getters) raise `RuntimeError` until `computeByDepth(...)`
completes successfully.

Extinction-value filtering and saliency maps use a dense node-indexed attribute
buffer and require a globally monotone altitude-order capability. Standard
max-tree and min-tree producers provide that capability:

```python
extinction = mmcfilters.ExtinctionValues(max_tree, level)
strongest = mmcfilters.ExtinctionSelectionPolicy.byTopK(8)
high_extinction = mmcfilters.ExtinctionSelectionPolicy.byThreshold(10.0)

filtered = extinction.filtering(strongest)
filtered_by_threshold = extinction.filtering(high_extinction)
contours = extinction.contourMap(
    strongest,
    mmcfilters.ExtinctionContourScorePolicy.RankScore,
)
extinction_tuples = extinction.getRegionalExtrema()
```

`ExtinctionSelectionPolicy.byTopK(...)` keeps the strongest extrema by decreasing
extinction ranking. `ExtinctionSelectionPolicy.byThreshold(...)` instead keeps
every regional extremum satisfying `extinction >= threshold`; the threshold must
be finite. The same selection policy is passed to `filtering(...)` and
`contourMap(...)`.

`ExtinctionValues.contourMap(...)` is an image-domain contour visualization. For
the QFZ-compatible hierarchical-watershed representation, use the persistence
path from Section 8.1 of Cousty:

```python
valuation = extinction.computeRankedExtinctionValueAttribute()
raw_edge_map = extinction.computeFormalSaliencyEdgeMap()
edge_map = extinction.computeFormalSaliencyEdgeMap(ranked=True)
```

`getExtinctionValueAttribute()` returns the raw cached node valuation, including
the dominant-extremum sentinel. The valuation extends each regional extremum
value from its leaf to all ancestors by subtree maximum; non-leaf nodes therefore
receive the largest extinction value among the extrema they contain. This is an
intermediate quantity: `computeFormalSaliencyEdgeMap(...)` builds an
altitude-ordered MST/BPTAO, assigns each binary merge the minimum of the two
child maximum extinctions, and extends the persistence-weighted MST to the full
graph. With `ranked=True`, only values present on the final edge map are ranked.

The former direct projection of the max-propagated node valuation remains
available under an explicit name:

```python
legacy_projection = extinction.computeMonotoneExtinctionProjection()
legacy_ranked = extinction.computeMonotoneExtinctionProjection(ranked=True)
```

That projection is monotone and QFZ-compatible, but it is not the Cousty
hierarchical-watershed persistence construction.

The count passed to `ExtinctionSelectionPolicy.byTopK(...)` must be
non-negative. The dominant extremum is reported with the
finite sentinel value `numpy.finfo(dtype).max`/`std::numeric_limits<Real>::max()`.
Standard trees of shapes and self-dual residual trees declare unconstrained
altitude order and are rejected by this API until their regional-extrema set is
computed explicitly. Validation uses the altitude-order capability, not the
descriptive tree kind.

Ultimate Attribute Opening also consumes a dense node-indexed attribute buffer:

```python
uao = mmcfilters.UltimateAttributeOpening(max_tree, box_height)
uao.execute(maxCriterion=image.shape[0])
uao.executeWithMSER(maxCriterion=image.shape[0], deltaMSER=2)
uao.executeWithDepthStability(maxCriterion=image.shape[0], depthDelta=2)

max_contrast = uao.getMaxContrastImage()
associated = uao.getAssociatedImage()
associated_color = uao.getAssociatedColoredImage()
```

## Contours

Contours are extracted once and materialized lazily. Iterating one contour may
cache that node without materializing all contours.

```python
contours = mmcfilters.ContourComputation.extraction(max_tree)

root_contour = list(contours.getContour(max_tree.getRoot()))

for node_id, contour in contours.contoursByNode():
    pixels = list(contour)
```

Call `materializeAll()` only when the workload will revisit many contours:

```python
contours.materializeAll()
```

### Contour Traces

`ContourTraceComputation` is the geometric contour API. It exposes side-level
boundary edges and ordered loops, separating external loops from internal loops.

```python
traces = mmcfilters.ContourTraceComputation.extraction(max_tree)

root_edges = traces.getEdges(max_tree.getRoot())
root_loops = traces.getLoops(max_tree.getRoot())

for loop in root_loops:
    edges = traces.getLoopEdges(loop)
```

Loop metadata includes `kind`, `edge_count`, and `signed_area2`. External loops
have positive signed area under the C++ orientation convention; internal loops
have negative signed area. The returned `ContourTraces` object keeps
`max_tree` alive for subsequent lazy queries.

## Higra Interoperability

`exportHigraHierarchy()` creates a compact Higra-style `(parent, altitude)`
layout for the current live tree. This exported layout is distinct from the
internal `NodeId` domain. The complete import/export contract is documented in
[Higra Interoperability](higra-interoperability.md).

```python
parent, altitude = max_tree.exportHigraHierarchy()

area_exported = max_tree.project_node_values_to_exported_higra(
    area,
    mmcfilters.Attribute.AREA,
)
```

Import preserves the original Higra node-id domain until the tree is edited.
For max-trees and min-trees, pass an explicit radius so adjacency-dependent
attributes remain available:

```python
roundtrip = mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
    parent,
    altitude,
    max_tree.numRows,
    max_tree.numCols,
    mmcfilters.MorphologicalTreeKind.MAX_TREE,
    radius=1.5,
)

area_in_imported_higra_space = mmcfilters.Attribute.computeSingleAttribute(
    roundtrip,
    mmcfilters.Attribute.AREA,
    mmcfilters.NodeIdSpace.HIGRA,
)
```

## CASF And Paired Tree Adjustment

The CASF wrapper builds paired component trees internally and applies thresholds
for the selected attribute:

```python
casf = mmcfilters.CasfComponentTrees(
    image,
    mmcfilters.CasfComponentTreesAttribute.AREA,
    radius=1.5,
)

casf_result = casf.filter([2.0, 4.0])
min_parent, min_altitude = casf.exportMinTree()
max_parent, max_altitude = casf.exportMaxTree()
```

For paired min/max adjustment, build compatible trees on the same image domain:

```python
adjust = mmcfilters.DualMinMaxTreeIncrementalFilter(min_tree, max_tree)
candidate_max_nodes = [
    node_id
    for node_id in max_tree.getAliveNodeIds()
    if node_id != max_tree.getRoot() and max_tree.getNumProperParts(node_id) <= 1
]

adjust.pruneMaxTreeAndUpdateMinTree(candidate_max_nodes[:1])
updated_min_tree = adjust.minTree
updated_max_tree = adjust.maxTree
```

## Failure Modes

Python bindings fail explicitly instead of silently converting unsupported
inputs. Common errors include:

- non-`uint8`, non-contiguous, or non-2D image arrays passed to factories;
- altitude arrays with the wrong dtype, shape, contiguity, length, or value
  range;
- topology-only attribute calls with altitude-dependent attributes;
- stale filtering helper objects used after topology mutation;
- `NodeIdSpace.HIGRA` requested on trees that do not preserve an imported Higra
  domain.
