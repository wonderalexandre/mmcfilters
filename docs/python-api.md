# Python API

This guide describes the Python-specific input, output, ownership, and naming
contracts of `mmcfilters`. The subsystem guides define the corresponding tree,
attribute, filter, contour, saliency, and interoperability semantics.

## Input contract

Python exposes one `ValuedMorphologicalTree` abstraction while retaining the
native altitude representation required by each hierarchy:

- image inputs are two-dimensional C-contiguous `np.uint8` arrays;
- component-tree and residual-tree altitude arrays are one-dimensional
  C-contiguous `np.uint8` arrays or integer sequences in `[0, 255]`;
- tree-of-shapes node altitudes follow the `altitude_encoding` declared by the
  topographic convention. The default `TopographicAltitudeEncoding.UINT8`
  publishes unchanged source gray levels as `np.uint8`;
  `TopographicAltitudeEncoding.EXACT_DOUBLED` publishes exact doubled gray units
  as `np.uint16`, where an original gray level `g` is represented by `2*g` and
  odd values represent exact half levels;
- factories return `ValuedMorphologicalTree`;
- a separate topology-only `MorphologicalTree` Python type is not exposed.

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
max_tree = mmcfilters.MorphologicalTreeFactory.create_max_tree(
    image,
    radius=1.5,
)
min_tree = mmcfilters.MorphologicalTreeFactory.create_min_tree(
    image,
    radius=1.5,
)
residual_tree = mmcfilters.MorphologicalTreeFactory.create_unrestricted_residual_tree(
    image,
    radius=1.5,
)
saturated_residual_tree = (
    mmcfilters.MorphologicalTreeFactory.create_saturated_residual_tree(
        image,
        infinity_pixel=0,
        radius=1.5,
    )
)
tree_of_shapes = mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(image)
assert tree_of_shapes.node_altitudes.dtype == np.uint8
assert np.array_equal(tree_of_shapes.reconstruct_from_node_altitudes(), image)

self_dual_tree_of_shapes = (
    mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(
        image,
        mmcfilters.self_dual_span_convention(),
    )
)
assert self_dual_tree_of_shapes.node_altitudes.dtype == np.uint16
assert np.array_equal(
    self_dual_tree_of_shapes.reconstruct_from_node_altitudes(),
    2 * image.astype(np.uint16),
)
```

The default convention selects the canonical minimum-4/maximum-8
complementary-grid immersion without domain padding, uses infinity pixel zero,
and publishes 8-bit altitudes. Its construction levels never leave the source
level set, so those altitudes are exact rather than quantized.

The self-dual span immersion admits `UINT8` only without an exterior ring. Its
boundary reference level is the single source of half levels, and
`TopographicDomainExtension.NONE` crops that level away before any interior cell
reads it. The convention fields may be passed straight to the factory:

```python
unpadded_self_dual = mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(
    image,
    immersion=mmcfilters.SelfDualSpanImmersion(),
    domain_extension=mmcfilters.TopographicDomainExtension.NONE,
    altitude_encoding=mmcfilters.TopographicAltitudeEncoding.UINT8,
)
assert unpadded_self_dual.node_altitudes.dtype == np.uint8
assert np.array_equal(unpadded_self_dual.reconstruct_from_node_altitudes(), image)
```

`mmcfilters.self_dual_span_convention(...)` builds the same convention as a
value when one has to be stored or reused. Passing both a complete convention
and individual fields to the factory is rejected, so a call always has one
source of truth.

Combining the self-dual span with `EXTERIOR_RING` and `UINT8` is rejected.

Max-trees declare `NodeAltitudeOrder.INCREASING`; min-trees declare
`NodeAltitudeOrder.DECREASING`. Trees of shapes and self-dual residual trees
declare `NodeAltitudeOrder.UNCONSTRAINED` because their altitudes do not have one global polarity.

The residual-tree factories use a shared symmetric adjacency and one canonical
`SelfDualResidualSchedule`. Candidate keys are
`SelfDualResidualKey(support_cardinality, spatial_minimum)`, with
`RowMajorSpatialOrder` by default. A custom `SpatialOrder` may be supplied; no
polarity-first tie policy is exposed.

### Native hierarchies

A hierarchy over an abstract finite pixel set does not require `rows` or
`columns`:

```python
semantics = mmcfilters.MorphologicalTreeSemantics()
abstract_tree = mmcfilters.MorphologicalTreeFactory.create_from_native_topology(
    parent=[0, 0, 0],
    smallest_node_map=[1, 2],
    node_altitudes=np.array([10, 3, 20], dtype=np.uint8),
    root=0,
    semantics=semantics,
)
assert not abstract_tree.has_grid_domain_2d
```

Pass `rows`, `columns`, and `semantics` to attach a regular 2D layout. Geometry does
not change the topology or smallest-node mapping.

Tree-of-shapes semantics retain a `TopographicConvention`. A complementary
immersion carries `ComplementaryAdjacencies(min_adjacency, max_adjacency)`;
self-dual span immersion is represented explicitly by `SelfDualSpanImmersion`.
Attribute requirements can be inspected before computation:

```python
requirements = mmcfilters.Attribute.requirements(
    mmcfilters.Attribute.CENTRAL_MOMENT_20
)
assert requirements["grid_domain_2d"]
assert not requirements["altitude"]
assert requirements["adjacency"] == "none"
assert not requirements["monotone_altitude_order"]
```

## Node IDs and queries

Tree queries, attributes, altitude buffers, predicates, and mutations use the
dense internal `NodeId` domain. Edits may leave dead slots, so dense results keep
`max_tree.num_internal_node_slots` rows. Iterate live IDs when appropriate:

```python
root = max_tree.root

for node_id in max_tree.alive_node_ids:
    children = max_tree.children(node_id)
    direct_pixels = max_tree.proper_part(node_id)
    subtree = max_tree.subtree_nodes(node_id)
    entry_index = max_tree.dfs_entry_index(node_id)
    exit_index = max_tree.dfs_exit_index(node_id)
```

The DFS entry and exit indices belong to one interleaved event sequence. They
characterize the node interval, so
`num_descendants(node_id) == (exit_index - entry_index - 1) // 2`.

Pixels use row-major indices. The canonical structural queries are:

```python
pixel_id = 10
smallest = max_tree.smallest_node(pixel_id)
support_pixels = list(max_tree.node_support(smallest))
component_mask = max_tree.reconstruct_node(smallest)
```

Python exposes safe local mutations and checked altitude setters:

```python
max_tree.prune_node(node_id)
max_tree.merge_node_into_parent(node_id)
max_tree.set_node_altitude(node_id, value)
max_tree.node_altitudes = altitude
```

`TreeEditor`, `ValuedMorphologicalTreeEditor`, unchecked setters, and a mutable topology
handle are not exposed in Python. See [Editing API](editing-api.md) for mutation
and derived-state lifetime.

## Attributes

Use valued-tree entry points when an attribute may read node altitude and
topology/support entry points otherwise. Single-attribute methods return a 1D
array; multi-attribute methods return `(layout, values)`, where `layout` maps
names to columns. `dtype` accepts `np.float32` or `np.float64` and defaults to
`np.float32`.

Every method that accepts an `Attribute` or `Attribute.Group` also accepts its
stable symbolic name as a string. The accepted names are exactly the ones that
`layout` already uses as keys, so the request and the result share one
vocabulary:

```python
area = mmcfilters.Attribute.compute_single_attribute(max_tree, "AREA")
layout, values = mmcfilters.Attribute.compute_attributes(
    max_tree,
    ["AREA", "VOLUME", "BOUNDARY"],
)
```

Matching is exact and case-sensitive. An unknown name raises `ValueError` and
reports the closest known names.

```python
area = mmcfilters.Attribute.compute_single_topology_attribute(
    max_tree,
    mmcfilters.Attribute.AREA,
)
altitude = max_tree.node_altitude(node_id)
gray_level_height = mmcfilters.Attribute.compute_single_attribute(
    max_tree,
    mmcfilters.Attribute.GRAY_LEVEL_HEIGHT,
)

layout, values = mmcfilters.Attribute.compute_attributes(
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
    mmcfilters.Attribute.compute_topology_attributes(
        max_tree,
        [mmcfilters.Attribute.Group.BOUNDARY],
    )
)
```

Stable scalar and group membership, input contracts, and family-specific guides
are indexed by the [Attribute catalog](attribute-catalog.md).

Altitude-based node-attribute sampling returns current-node, ancestor, and
representative-descendant columns. Missing values use a typed policy:

```python
sample_layout, sampled_values = (
    mmcfilters.Attribute.compute_sampled_node_attribute(
        max_tree,
        mmcfilters.Attribute.AREA,
        1,
        1,
        missing_sample_policy=(
            mmcfilters.MissingNodeAttributeSamplePolicy.NOT_A_NUMBER
        ),
    )
)
```

The two integer arguments are `altitude_step` and `sampling_radius`.
Serialized layout labels use `_ANCESTOR_n` and `_DESCENDANT_n`.

Use `compute_attribute_mapping` to project a node attribute to the image domain.
The complete attribute and layout contracts are in [Attributes](attributes.md)
and the [Attribute catalog](attribute-catalog.md).

## Filters and stability

Filter objects capture the topology mutation version at construction. Create
them after the edits that should affect the result and do not reuse them after a
topology mutation.

```python
filters = mmcfilters.AttributeFilters(max_tree)

box_height = mmcfilters.Attribute.compute_single_topology_attribute(
    max_tree,
    mmcfilters.Attribute.BOUNDING_BOX_HEIGHT,
)
keep_large = mmcfilters.compute_node_preservation_mask(area, 4.0)

direct = mmcfilters.DirectAttributeFilter(max_tree).apply_direct_attribute_filter(
    keep_large,
)
subtractive = mmcfilters.SubtractiveAttributeFilter(
    max_tree,
).apply_subtractive_attribute_filter(keep_large)
pruned_min = filters.filtering_by_pruning_min(box_height, 2.0)
pruned_max = filters.filtering_by_pruning_max(box_height, 2.0)
altitude_adjusted = mmcfilters.adjust_node_preservation_mask_by_altitude_stability(
    max_tree,
    keep_large,
    altitude_window_radius=2,
)
depth_adjusted = mmcfilters.adjust_node_preservation_mask_by_depth_stability(
    max_tree,
    keep_large,
    depth_window_radius=2,
)
```

Use `DepthStableRegionComputer` when the variation scores or selected mask are
needed directly:

```python
depth = mmcfilters.DepthStableRegionComputer(max_tree)
depth_mask = depth.compute_by_depth(depth_window_radius=2)
depth_variation = depth.get_variations()
```

Its result getters raise `RuntimeError` until `compute_by_depth` succeeds.

Extinction selection and Ultimate Attribute Opening (UAO) use the same dense
attribute-buffer convention:

```python
extinction = mmcfilters.ExtinctionValues(max_tree, area)
strongest = mmcfilters.ExtinctionSelectionPolicy.by_top_k(8)
filtered = extinction.filtering(strongest)

uao = mmcfilters.UltimateAttributeOpening(max_tree, box_height)
uao.execute(maximum_attribute_threshold=image.shape[0])
max_contrast = uao.get_max_contrast_image()
associated = uao.get_associated_image()
```

See [Filters](filters.md) for rule selection, extinction contracts, stability,
and UAO outputs.

## Saliency maps

The Python surface exposes three distinct operations:

| Need | API |
| --- | --- |
| project a monotone hierarchy valuation | `HierarchySaliencyMap` |
| compute a persistence hierarchy from extinction values | `ExtinctionValues.compute_formal_saliency_edge_map` |
| compute extinction in tree-node shape space | `ShapeSpaceSaliency` |

```python
edge_map = mmcfilters.HierarchySaliencyMap.compute_normalized_altitude_edge_map(
    max_tree,
)

shape_space = mmcfilters.ShapeSpaceSaliency.compute(
    max_tree,
    area,
    mmcfilters.ShapeSpaceExtremaPolarity.MINIMA,
)
```

Edge-map dictionaries contain `sources`, `targets`, `values`, `num_rows`,
`num_columns`, and `adjacency_radius`. Cuts and display projections are provided by
`HierarchySaliencyMapProjection`. See [Saliency maps](saliency.md) for operator
definitions and preconditions.

## Contours

Pixel contours are materialized lazily:

```python
contours = mmcfilters.ContourComputation.extraction(max_tree)
root_contour = list(contours.get_contour(max_tree.root))

for node_id, contour in contours.contours_by_node():
    pixels = list(contour)
```

Use contour traces for oriented sides and ordered external/internal loops:

```python
traces = mmcfilters.ContourTraceComputation.extraction(max_tree)
root_edges = traces.get_edges(max_tree.root)
root_loops = traces.get_loops(max_tree.root)

for loop in root_loops:
    edges = traces.get_loop_edges(loop)
```

See [Pixel contours](contours.md) and [Contour traces](contour-traces.md).

## Higra interoperability

`export_higra_hierarchy()` returns a compact `(parent, altitude)` snapshot of the
live tree. Its node domain differs from internal `NodeId` values:

```python
parent, altitude = max_tree.export_higra_hierarchy()
area_exported = max_tree.project_node_values_to_exported_higra(
    area,
    mmcfilters.Attribute.AREA,
)
```

Import preserves the supplied Higra domain until the topology is edited:

```python
roundtrip = mmcfilters.MorphologicalTreeFactory.create_from_higra_parent(
    parent,
    altitude,
    max_tree.num_rows,
    max_tree.num_columns,
    mmcfilters.MorphologicalTreeKind.MAX_TREE,
    radius=1.5,
)
area_in_higra_space = mmcfilters.Attribute.compute_single_attribute(
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
min_parent, min_altitude = casf.export_min_tree()
max_parent, max_altitude = casf.export_max_tree()
```

Paired adjustment requires `ValuedMorphologicalTree` instances for a min-tree
and max-tree over the same image domain:

```python
adjust = mmcfilters.DualMinMaxTreeIncrementalFilter(min_tree, max_tree)
candidate = next(
    node_id
    for node_id in max_tree.alive_node_ids
    if node_id != max_tree.root
)
adjust.prune_max_tree_and_update_min_tree([candidate])
updated_min_tree = adjust.min_tree
updated_max_tree = adjust.max_tree
```

## Failure modes

Python bindings reject:

- non-`np.uint8`, non-contiguous, or non-2D factory images;
- altitude arrays with an invalid dtype, shape, contiguity, length, or range;
- topology/support calls that request altitude-dependent attributes;
- non-contiguous attribute arrays passed to filter helpers;
- cached helper objects used after topology mutation;
- `NodeIdSpace.HIGRA` on trees without a valid imported Higra domain.
