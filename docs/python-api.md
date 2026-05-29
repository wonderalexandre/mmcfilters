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
- direct image construction through `MorphologicalTree(...)` is not public API.

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
tos = mmcfilters.MorphologicalTreeFactory.createTreeOfShapes(image)
```

The returned object owns both topology and a dense altitude buffer:

```python
root = max_tree.getRoot()
alive_nodes = max_tree.getAliveNodeIds()
leaf_nodes = max_tree.getLeafNodeIds()
altitude_at_root = max_tree.getAltitude(root)
reconstructed = max_tree.reconstructionImage()
```

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

Extinction-value filtering and saliency maps use a dense node-indexed attribute
buffer:

```python
extinction = mmcfilters.ExtinctionValues(max_tree, level)
filtered = extinction.filtering(leafToKeep=8)
saliency = extinction.saliencyMap(leafToKeep=8)
extinction_tuples = extinction.getExtinctionValues()
```

Ultimate Attribute Opening also consumes a dense node-indexed attribute buffer:

```python
uao = mmcfilters.UltimateAttributeOpening(max_tree, box_height)
uao.execute(maxCriterion=image.shape[0])

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
