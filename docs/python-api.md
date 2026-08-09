# Python API

This guide describes the Python-specific input, output, ownership, and naming
contracts of `mmcfilters`. The subsystem guides define the corresponding tree,
attribute, filter, contour, saliency, and interoperability semantics.

## Input contract

Python construction stores 8-bit altitudes in `WeightedMorphologicalTree`:

- image inputs are two-dimensional C-contiguous `np.uint8` arrays;
- altitude arrays are one-dimensional C-contiguous `np.uint8` arrays or integer
  sequences in `[0, 255]`;
- factories return `WeightedMorphologicalTree`;
- a separate unweighted `MorphologicalTree` Python type is not exposed.

Normalize an image at the application boundary:

```python
import numpy as np
import mmcfilters

image = np.ascontiguousarray(
    [
        [3, 3, 2, 2],
        [3, 4, 4, 2],
        [1, 4, 5, 2],
        [1, 1, 5, 0],
    ],
    dtype=np.uint8,
)
```

Methods that accept node attributes require one-dimensional C-contiguous
`np.float32` or `np.float64` arrays with one entry per internal node slot.

## Building trees

On a 2D square grid, `radius=1.0` selects 4-connectivity and `radius=1.5`
selects 8-connectivity:

```python
max_tree = mmcfilters.MorphologicalTreeFactory.createMaxTree(
    image,
    radius=1.5,
)
min_tree = mmcfilters.MorphologicalTreeFactory.createMinTree(
    image,
    radius=1.5,
)
residual_tree = mmcfilters.MorphologicalTreeFactory.createSelfDualResidualTree(
    image,
    radius=1.5,
)
saturated_residual_tree = (
    mmcfilters.MorphologicalTreeFactory.createSaturatedSelfDualResidualTree(
        image,
        infinityPixel=0,
        radius=1.5,
    )
)
tree_of_shapes = mmcfilters.MorphologicalTreeFactory.createTreeOfShapes(image)
```

Max-trees declare `AltitudeOrder.INCREASING_FROM_ROOT`; min-trees declare
`DECREASING_FROM_ROOT`. Tree of Shapes hierarchies and self-dual residual trees
declare `UNCONSTRAINED` because their altitudes do not have one global polarity.

The residual-tree factories use a shared symmetric adjacency.
`SdrtTiePolicy.CONTRAST_INVARIANT_SPATIAL` is their default deterministic tie
policy.

### Native hierarchies

A hierarchy over an abstract finite proper-part set does not require `rows` or
`cols`:

```python
semantics = mmcfilters.HierarchySemantics()
abstract_tree = mmcfilters.MorphologicalTreeFactory.createFromNativeTopology(
    [0, 0, 0],
    [1, 2],
    np.array([10, 3, 20], dtype=np.uint8),
    0,
    semantics=semantics,
)
assert not abstract_tree.hasGridDomain2D
```

Pass `rows`, `cols`, and `semantics` to attach a regular 2D layout. Geometry does
not change the topology or proper-part ownership.

Directional semantics use
`DirectionalGridAdjacency2D(decreasing, increasing)`. Attribute requirements
can be inspected before computation:

```python
requirements = mmcfilters.Attribute.requirements(
    mmcfilters.Attribute.MAX_DIST
)
assert requirements["gridDomain2D"]
assert requirements["adjacency"] == "uniform"
assert requirements["monotoneAltitudeOrder"]
```

## Node IDs and queries

Tree queries, attributes, altitude buffers, predicates, and mutations use the
dense internal `NodeId` domain. Edits may leave dead slots, so dense results keep
`max_tree.numInternalNodeSlots` rows. Iterate live IDs when appropriate:

```python
root = max_tree.getRoot()

for node_id in max_tree.getAliveNodeIds():
    children = max_tree.getChildren(node_id)
    direct_pixels = max_tree.getProperParts(node_id)
    subtree = max_tree.getNodeSubtree(node_id)
```

Pixels are row-major proper parts:

```python
pixel_id = 10
owner = max_tree.getProperPartOwner(pixel_id)
component_pixels = list(max_tree.getConnectedComponent(owner))
component_mask = max_tree.reconstructNode(owner)
```

Python exposes safe local mutations and checked altitude setters:

```python
max_tree.pruneNode(node_id)
max_tree.mergeNodeIntoParent(node_id)
max_tree.setAltitude(node_id, value)
max_tree.setAltitudeBuffer(altitude)
```

`TreeEditor`, `WeightedTreeEditor`, unchecked setters, and a mutable topology
handle are not exposed in Python. See [Editing API](editing-api.md) for mutation
and derived-state lifetime.

## Attributes

Use weighted entry points when an attribute may read altitude and
topology/support entry points otherwise. Single-attribute methods return a 1D
array; multi-attribute methods return `(layout, values)`, where `layout` maps
names to columns. `dtype` accepts `np.float32` or `np.float64` and defaults to
`np.float32`.

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
```

A sliced column is generally not C-contiguous. Copy it before passing it to a
filter helper, or compute the scalar attribute directly as above.

Groups expand to stable scalar sets:

```python
boundary_layout, boundary_values = (
    mmcfilters.Attribute.computeTopologyAttributes(
        max_tree,
        [mmcfilters.Attribute.Group.BOUNDARY],
    )
)
```

Delta sampling returns center, ancestor, and descendant columns. Missing values
use the selected padding policy:

```python
delta_layout, delta_values = (
    mmcfilters.Attribute.computeSingleAttributeWithDelta(
        max_tree,
        mmcfilters.Attribute.AREA,
        1,
        "null-padding",
    )
)
```

Use `computeAttributeMapping` to project a node attribute to the image domain.
The complete attribute and layout contracts are in [Attributes](attributes.md)
and the [Attribute catalog](attribute-catalog.md).

## Filters and stability

Filter objects capture the topology mutation version at construction. Create
them after the edits that should affect the result and do not reuse them after a
topology mutation.

```python
filters = mmcfilters.AttributeFilters(max_tree)

box_height = mmcfilters.Attribute.computeSingleTopologyAttribute(
    max_tree,
    mmcfilters.Attribute.BOX_HEIGHT,
)
keep_large = (area >= 4.0).tolist()

direct = filters.filteringDirectRule(keep_large)
subtractive = filters.filteringSubtractiveRule(keep_large)
pruned_min = filters.filteringByPruningMin(box_height, 2.0)
pruned_max = filters.filteringByPruningMax(box_height, 2.0)
adaptive = filters.getAdaptiveCriterion(keep_large, delta=2)
depth_adaptive = filters.getAdaptiveCriterionByDepth(
    keep_large,
    depthDelta=2,
)
```

Use `DepthStableRegionComputer` when the variation scores or selected mask are
needed directly:

```python
depth = mmcfilters.DepthStableRegionComputer(max_tree)
depth_mask = depth.computeByDepth(depthDelta=2)
depth_variation = depth.getVariations()
```

Its result getters raise `RuntimeError` until `computeByDepth` succeeds.

Extinction selection and Ultimate Attribute Opening (UAO) use the same dense
attribute-buffer convention:

```python
extinction = mmcfilters.ExtinctionValues(max_tree, area)
strongest = mmcfilters.ExtinctionSelectionPolicy.byTopK(8)
filtered = extinction.filtering(strongest)

uao = mmcfilters.UltimateAttributeOpening(max_tree, box_height)
uao.execute(maxCriterion=image.shape[0])
max_contrast = uao.getMaxContrastImage()
associated = uao.getAssociatedImage()
```

See [Filters](filters.md) for rule selection, extinction contracts, stability,
and UAO outputs.

## Saliency maps

The Python surface exposes three distinct operations:

| Need | API |
| --- | --- |
| project a monotone hierarchy valuation | `HierarchySaliencyMap` |
| compute a persistence hierarchy from extinction values | `ExtinctionValues.computeFormalSaliencyEdgeMap` |
| compute extinction in tree-node shape space | `ShapeSpaceSaliency` |

```python
edge_map = mmcfilters.HierarchySaliencyMap.computeNormalizedAltitudeEdgeMap(
    max_tree,
)

shape_space = mmcfilters.ShapeSpaceSaliency.compute(
    max_tree,
    area,
    mmcfilters.ShapeSpaceExtremaPolarity.Minima,
)
```

Edge-map dictionaries contain `sources`, `targets`, `values`, `numRows`,
`numCols`, and `adjacencyRadius`. Cuts and display projections are provided by
`HierarchySaliencyMapProjection`. See [Saliency maps](saliency.md) for operator
definitions and preconditions.

## Contours

Pixel contours are materialized lazily:

```python
contours = mmcfilters.ContourComputation.extraction(max_tree)
root_contour = list(contours.getContour(max_tree.getRoot()))

for node_id, contour in contours.contoursByNode():
    pixels = list(contour)
```

Use contour traces for oriented sides and ordered external/internal loops:

```python
traces = mmcfilters.ContourTraceComputation.extraction(max_tree)
root_edges = traces.getEdges(max_tree.getRoot())
root_loops = traces.getLoops(max_tree.getRoot())

for loop in root_loops:
    edges = traces.getLoopEdges(loop)
```

See [Pixel contours](contours.md) and [Contour traces](contour-traces.md).

## Higra interoperability

`exportHigraHierarchy()` returns a compact `(parent, altitude)` snapshot of the
live tree. Its node domain differs from internal `NodeId` values:

```python
parent, altitude = max_tree.exportHigraHierarchy()
area_exported = max_tree.project_node_values_to_exported_higra(
    area,
    mmcfilters.Attribute.AREA,
)
```

Import preserves the supplied Higra domain until the topology is edited:

```python
roundtrip = mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
    parent,
    altitude,
    max_tree.numRows,
    max_tree.numCols,
    mmcfilters.MorphologicalTreeKind.MAX_TREE,
    radius=1.5,
)
area_in_higra_space = mmcfilters.Attribute.computeSingleAttribute(
    roundtrip,
    mmcfilters.Attribute.AREA,
    mmcfilters.NodeIdSpace.HIGRA,
)
```

See [Higra interoperability](higra-interoperability.md) for the two external
node domains and projection rules.

## Connected alternating sequential filter

The connected alternating sequential filter (CASF) API is exposed as
`CasfComponentTrees`. It constructs paired component trees and applies
thresholds for a selected attribute:

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

Paired adjustment requires `WeightedMorphologicalTree` instances for a min-tree
and max-tree over the same image domain:

```python
adjust = mmcfilters.DualMinMaxTreeIncrementalFilter(min_tree, max_tree)
candidate = next(
    node_id
    for node_id in max_tree.getAliveNodeIds()
    if node_id != max_tree.getRoot()
)
adjust.pruneMaxTreeAndUpdateMinTree([candidate])
updated_min_tree = adjust.minTree
updated_max_tree = adjust.maxTree
```

## Failure modes

Python bindings reject:

- non-`np.uint8`, non-contiguous, or non-2D factory images;
- altitude arrays with an invalid dtype, shape, contiguity, length, or range;
- topology/support calls that request altitude-dependent attributes;
- non-contiguous attribute arrays passed to filter helpers;
- cached helper objects used after topology mutation;
- `NodeIdSpace.HIGRA` on trees without a valid imported Higra domain.
