# Higra interoperability

This guide documents the Higra-facing boundary of `mmcfilters`: importing a
static hierarchy, preserving imported Higra node IDs, exporting the current
tree, and projecting node attributes to exported layouts.

For the underlying tree model, see [Morphological trees](trees.md). For
attribute-buffer layout, see [Attributes](attributes.md).

## Two Higra domains

The library exposes two distinct Higra-related domains:

| Domain | Created by | Edit-stable | Main use |
| --- | --- | --- | --- |
| Preserved imported Higra domain | `createFromHigraParent(...)` | No | Original imported node IDs |
| Exported compact Higra domain | `exportHigraHierarchy()` | Snapshot only | Current live tree export |

`NodeIdSpace::Higra` selects the original imported node IDs while that domain
remains valid. Normal tree queries, filters, and contours use internal
`NodeIdSpace::MorphologicalTree` indexing.

`exportHigraHierarchy()` always computes a new compact layout for the current
live rooted tree. Use `projectNodeValuesToExportedHigra(...)` or Python
`project_node_values_to_exported_higra(...)` when attributes must be aligned
with that exported snapshot.

## Compact layout

The compact Higra layout used by import and export is:

```text
[pixel leaves | internal nodes]
```

For an image domain with `rows * columns` pixels:

- pixel leaf IDs are `[0, rows * columns)`;
- internal node IDs are `[rows * columns, parent.size())`;
- every pixel leaf points to its smallest node in the external layout;
- every internal node points to another internal node or to itself;
- exactly one internal node is self-parented and is the root.

Each pixel leaf represents one pixel. A tree node's proper part is the set of
pixels mapped to that node; it may contain several pixels or be empty.

When exporting, pixel leaves are emitted in row-major order. Internal nodes are
assigned compact node IDs from the live rooted tree. For max-trees and min-trees,
the export order follows the tree altitude polarity with a deterministic
post-order tie-breaker. Trees of shapes and other
`NodeAltitudeOrder::Unconstrained` hierarchies use deterministic post-order directly. Consequently,
every non-root internal node appears before its parent even when one branch
increases in altitude and another decreases.

The exported altitude array has the same length as the exported parent array.
Each pixel leaf receives the altitude of its smallest node.

This layout policy lives at the interoperability boundary. Import converts it
to separate dense buffers for node parents and smallest nodes before generic
tree materialization; `MorphologicalTree` does not parse Higra parent arrays.
While the topology is unchanged, it retains only the affine external ID offset
needed by the imported-ID queries below.

## Importing a static hierarchy

Use `MorphologicalTreeFactory::createFromHigraParent(...)` in C++:

```cpp
#include <cstdint>
#include <span>
#include <vector>

#include <mmcfilters/trees/MorphologicalTreeFactory.hpp>

using namespace mmcfilters;

std::vector<NodeId> parent = /* compact [pixel leaves | internal nodes] */;
std::vector<std::uint8_t> altitude = /* same size as parent */;

auto tree = MorphologicalTreeFactory::createFromHigraParent(
    std::span<const NodeId>(parent),
    std::span<const std::uint8_t>(altitude),
    rows,
    columns,
    MorphologicalTreeKind::MaxTree,
    RegularGridAdjacency2D(rows, columns, 1.5));
```

C++ import is generic over the altitude type:

```cpp
std::vector<float> floatAltitude(parent.size(), 0.0f);

auto floatTree = MorphologicalTreeFactory::createFromHigraParent<float>(
    std::span<const NodeId>(parent),
    std::span<const float>(floatAltitude),
    rows,
    columns,
    MorphologicalTreeKind::MinTree,
    RegularGridAdjacency2D(rows, columns, 1.5));
```

Python exposes the canonical 8-bit path:

```python
import mmcfilters

tree = mmcfilters.MorphologicalTreeFactory.create_from_higra_parent(
    parent=parent,
    node_altitudes=node_altitudes,
    rows=rows,
    columns=columns,
    kind=mmcfilters.MorphologicalTreeKind.MAX_TREE,
    radius=1.5,
)
```

Python altitude inputs must be integer sequences in `[0, 255]` or 1D
C-contiguous `np.uint8` arrays. C++ accepts any type satisfying the public
`AltitudeValue` contract.

Max-tree and min-tree imports require adjacency metadata. In Python, pass
`radius`; in C++, pass a `RegularGridAdjacency2D`. Tree-of-shapes imports can
omit component-tree adjacency.

## Preserved imported node IDs

After import, the tree still uses the internal dense `NodeId` domain for normal
tree operations. A preserved mapping from internal live nodes to the original
Higra node IDs is available until the topology is edited:

```cpp
NodeId internal = tree.topology().root();
NodeId originalHigraId = tree.topology().getHigraNodeId(internal);
int higraDomainSize = tree.topology().getNumHigraNodes();
```

```python
internal = tree.root
original_higra_id = tree.higra_node_id(internal)
higra_domain_size = tree.num_higra_nodes
```

For trees imported from the compact layout above, the internal slot associated
with Higra internal node ID `h` starts as:

```text
node_id = h - rows * columns
```

Do not rely on that arithmetic after edits. Safe public code should use
`getHigraNodeId(...)` while the preserved domain is still valid.

Any topology mutation invalidates the preserved imported Higra domain. This
includes safe mutators such as `pruneNode(...)` and `mergeNodeIntoParent(...)`
and staged edit commits. After invalidation:

- `getNumHigraNodes()` fails;
- `NodeIdSpace::Higra` attribute requests fail;
- `getHigraNodeId(node)` returns `InvalidNode`.

Export still works after edits because it creates a new compact domain.

## Attributes in preserved Higra space

Attribute computation always runs internally in `NodeIdSpace::MorphologicalTree`.
Projection to `NodeIdSpace::Higra` is an API-boundary step.

Use preserved Higra output space when a consumer needs rows indexed by the
original imported node IDs:

```cpp
auto computed = AttributeComputation::computeSingleAttribute(
    tree,
    Area,
    NodeIdSpace::Higra);

const AttributeNames& names = computed.attributeNames();
const std::vector<float>& values = computed.values();
```

```python
area_in_imported_space = mmcfilters.Attribute.compute_single_attribute(
    tree,
    mmcfilters.Attribute.AREA,
    mmcfilters.NodeIdSpace.HIGRA,
)
```

Live internal-node rows receive the values computed in the internal
`MorphologicalTree` node ID space. Rows for pixel leaves in the preserved
imported Higra domain receive the same unit-component values used by compact
Higra export.

## Exporting the current tree

Use `exportHigraHierarchy()` when a consumer needs the current live tree:

```cpp
auto [exportedParent, exportedAltitude] = tree.exportHigraHierarchy();
```

```python
exported_parent, exported_altitude = tree.export_higra_hierarchy()
```

The exported layout is a fresh snapshot. It is valid for image-built trees,
imported trees, and edited trees as long as the current topology is one rooted
live component, no edit session is open, and the altitude buffer covers the
internal node slots.

Dead internal slots are not exported. Export size is:

```text
num_pixels + num_live_internal_nodes
```

The exported parent and altitude arrays should be treated as a pair. If the tree
is edited again, export again and reproject any attributes that must align with
the new compact node IDs.

## Projecting attributes to exported layout

To align attributes with `exportHigraHierarchy()`, project a dense internal-node
buffer:

```cpp
auto [names, values] = AttributeComputation::computeAttributes(
    tree,
    std::vector<AttributeOrGroup>{Area, MaxDist});

std::vector<float> exportedValues =
    AttributeComputation::projectNodeValuesToExportedHigra(
        tree,
        names,
        values);
```

Python exposes the same operation on `ValuedMorphologicalTree`:

```python
layout, values = mmcfilters.Attribute.compute_attributes(
    tree,
    [mmcfilters.Attribute.AREA, mmcfilters.Attribute.MAX_DIST],
)

exported_values = tree.project_node_values_to_exported_higra(
    values,
    [mmcfilters.Attribute.AREA, mmcfilters.Attribute.MAX_DIST],
)
```

For a single attribute:

```python
area = mmcfilters.Attribute.compute_single_attribute(
    tree,
    mmcfilters.Attribute.AREA,
)

exported_area = tree.project_node_values_to_exported_higra(
    area,
    mmcfilters.Attribute.AREA,
)
```

The projection output follows the same `[pixel leaves | live internal nodes]`
layout as the exported hierarchy. Internal-node rows are copied from the
node-indexed input. Pixel rows are filled with unit-component values for
the requested attributes.

Examples of values for a pixel leaf:

- `AREA`: `1`;
- `MeanGrayLevel`: the altitude of the smallest node;
- `VOLUME`: one pixel at the smallest node altitude;
- `GrayLevelVariance`, `GrayLevelHeight`, and `MAX_DIST`: `0`;
- bounding-box coordinates: the pixel's row and column; width and height: `1`.

Projection fails if:

- the node-value buffer shape does not match the dense internal-node domain;
- the attribute list does not match the number of columns;
- an attribute has no registered unit-component projection;
- the tree is currently inside an edit session;
- a borrowed `ValuedMorphologicalTreeView<T>` became stale after topology mutation.

## Round-trip pattern

A common interoperability round trip is:

```python
tree = mmcfilters.MorphologicalTreeFactory.create_max_tree(image, radius=1.5)

parent, altitude = tree.export_higra_hierarchy()

roundtrip = mmcfilters.MorphologicalTreeFactory.create_from_higra_parent(
    parent,
    altitude,
    image.shape[0],
    image.shape[1],
    mmcfilters.MorphologicalTreeKind.MAX_TREE,
    radius=1.5,
)

area_exported = tree.project_node_values_to_exported_higra(
    mmcfilters.Attribute.compute_single_attribute(tree, mmcfilters.Attribute.AREA),
    mmcfilters.Attribute.AREA,
)

area_imported_space = mmcfilters.Attribute.compute_single_attribute(
    roundtrip,
    mmcfilters.Attribute.AREA,
    mmcfilters.NodeIdSpace.HIGRA,
)
```

`area_exported` and `area_imported_space` use the same compact node ID layout
after the round trip. Both paths fill pixel rows with unit-component values.

## Related guides

- [Morphological trees](trees.md): owning tree/view boundary, `NodeId`, proper
  parts, altitude, and edits.
- [Attributes](attributes.md): attribute layouts, output spaces, and unit
  export projections.
- [Python API](python-api.md): Python construction and wrapper names.
- [Editing API](editing-api.md): edit-session lifetime and derived-state
  invalidation.
