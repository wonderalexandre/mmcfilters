# Higra Interoperability

This guide documents the Higra-facing boundary of `mmcfilters`: importing a
static hierarchy, preserving imported Higra node ids, exporting the current
tree, and projecting node attributes to exported layouts.

For the underlying tree model, see [Morphological Trees](trees.md). For
attribute buffer layout, see [Attributes](attributes.md).

## Two Higra Domains

The library exposes two distinct Higra-related domains:

| Domain | Created by | Edit-stable | Main use |
| --- | --- | --- | --- |
| Preserved imported Higra domain | `createFromHigraParent(...)` | No | Original imported ids |
| Exported compact Higra domain | `exportHigraHierarchy()` | Snapshot only | Current live tree export |

`NodeIdSpace::HIGRA` refers only to the preserved imported domain. It does not
mean "whatever `exportHigraHierarchy()` would produce now".

`exportHigraHierarchy()` always computes a new compact layout for the current
live rooted tree. Use `projectNodeValuesToExportedHigra(...)` or Python
`project_node_values_to_exported_higra(...)` when attributes must be aligned
with that exported snapshot.

## Compact Layout

The compact Higra layout used by import and export is:

```text
[proper parts | internal nodes]
```

For an image domain with `rows * cols` proper parts:

- proper-part leaf ids are `[0, rows * cols)`;
- internal ids are `[rows * cols, parent.size())`;
- every leaf points to an internal node;
- every internal node points to another internal node or to itself;
- exactly one internal node is self-parented and is the root.

When exporting, proper parts are emitted in row-major order. Internal nodes are
assigned compact ids from the live rooted tree. For max-trees and min-trees, the
export order follows the tree altitude polarity with a deterministic post-order
tie-breaker. Trees of shapes and other
`UNCONSTRAINED` hierarchies use deterministic post-order directly. Consequently,
every non-root internal node appears before its parent even when one branch
increases in altitude and another decreases.

The exported altitude array has the same length as the exported parent array.
Leaf/proper-part altitudes are filled with the altitude of their owner node.

This layout policy lives at the interoperability boundary. Import converts it
to independent dense node-parent and proper-part-owner buffers before generic
tree materialization; `MorphologicalTree` does not parse Higra parent arrays.
While the topology is unchanged, it retains only the affine external-id offset
needed by the compatibility queries below.

## Importing A Static Hierarchy

Use `MorphologicalTreeFactory::createFromHigraParent(...)` in C++:

```cpp
#include <cstdint>
#include <span>
#include <vector>

#include <mmcfilters/trees/MorphologicalTreeFactory.hpp>

using namespace mmcfilters;

std::vector<NodeId> parent = /* compact [proper parts | internal nodes] */;
std::vector<std::uint8_t> altitude = /* same size as parent */;

auto tree = MorphologicalTreeFactory::createFromHigraParent(
    std::span<const NodeId>(parent),
    std::span<const std::uint8_t>(altitude),
    rows,
    cols,
    MorphologicalTreeKind::MAX_TREE,
    RegularGridAdjacency2D(rows, cols, 1.5));
```

The altitude type is typed in C++:

```cpp
std::vector<float> floatAltitude(parent.size(), 0.0f);

auto floatTree = MorphologicalTreeFactory::createFromHigraParent<float>(
    std::span<const NodeId>(parent),
    std::span<const float>(floatAltitude),
    rows,
    cols,
    MorphologicalTreeKind::MIN_TREE,
    RegularGridAdjacency2D(rows, cols, 1.5));
```

Python exposes the canonical 8-bit path:

```python
import mmcfilters

tree = mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
    parent,
    altitude,
    rows,
    cols,
    mmcfilters.MorphologicalTreeKind.MAX_TREE,
    radius=1.5,
)
```

Python altitude inputs must be integer sequences in `[0, 255]` or 1D
C-contiguous `np.uint8` arrays. C++ accepts any type satisfying the public
`AltitudeValue` contract.

Max-tree and min-tree imports require adjacency metadata. In Python, pass
`radius`; in C++, pass a `RegularGridAdjacency2D`. Tree-of-shapes imports can omit
component-tree adjacency.

## Preserved Imported Ids

After import, the tree still uses the internal dense `NodeId` domain for normal
tree operations. A preserved mapping from internal live nodes to the original
Higra ids is available until the topology is edited:

```cpp
NodeId internal = tree.topology().getRoot();
NodeId originalHigraId = tree.topology().getHigraNodeId(internal);
int higraDomainSize = tree.topology().getNumHigraNodes();
```

```python
internal = tree.getRoot()
original_higra_id = tree.getHigraNodeId(internal)
higra_domain_size = tree.numHigraNodes
```

For trees imported from the compact layout above, the internal slot associated
with Higra internal id `h` starts as:

```text
node_id = h - rows * cols
```

Do not rely on that arithmetic after edits. Safe public code should use
`getHigraNodeId(...)` while the preserved domain is still valid.

Any topology mutation invalidates the preserved imported Higra domain. This
includes safe mutators such as `pruneNode(...)` and `mergeNodeIntoParent(...)`
and staged edit commits. After invalidation:

- `getNumHigraNodes()` fails;
- `NodeIdSpace::HIGRA` attribute requests fail;
- `getHigraNodeId(node)` returns `InvalidNode`.

Export still works after edits because it creates a new compact domain.

## Attributes In Preserved Higra Space

Attribute computation always runs internally in `NodeIdSpace::MORPHOLOGICAL_TREE`.
Projection to `NodeIdSpace::HIGRA` is an API-boundary step.

Use preserved Higra output space when a consumer needs rows indexed by the
original imported ids:

```cpp
auto computed = AttributeComputation::computeSingleAttribute(
    tree,
    AREA,
    NodeIdSpace::HIGRA);

const AttributeNames& names = computed.attributeNames();
const std::vector<float>& values = computed.values();
```

```python
area_in_imported_space = mmcfilters.Attribute.computeSingleAttribute(
    tree,
    mmcfilters.Attribute.AREA,
    mmcfilters.NodeIdSpace.HIGRA,
)
```

Live internal-node rows receive the values computed in the internal
`MorphologicalTree` node-id space. Proper-part/leaf rows in the preserved
imported Higra domain receive the same unit-component values used by compact
Higra export.

## Exporting The Current Tree

Use `exportHigraHierarchy()` when a consumer needs the current live tree:

```cpp
auto [exportedParent, exportedAltitude] = tree.exportHigraHierarchy();
```

```python
exported_parent, exported_altitude = tree.exportHigraHierarchy()
```

The exported layout is a fresh snapshot. It is valid for image-built trees,
imported trees, and edited trees as long as the current topology is one rooted
live component, no edit session is open, and the altitude buffer covers the
internal node slots.

Dead internal slots are not exported. Export size is:

```text
num_total_proper_parts + num_live_internal_nodes
```

The exported parent and altitude arrays should be treated as a pair. If the tree
is edited again, export again and reproject any attributes that must align with
the new compact ids.

## Projecting Attributes To Exported Layout

`NodeIdSpace::HIGRA` is not the right tool for exported snapshots. To align
attributes with `exportHigraHierarchy()`, project a dense internal-node buffer:

```cpp
auto [names, values] = AttributeComputation::computeAttributes(
    tree,
    std::vector<AttributeOrGroup>{AREA, MAX_DIST});

std::vector<float> exportedValues =
    AttributeComputation::projectNodeValuesToExportedHigra(
        tree,
        names,
        values);
```

Python exposes the same operation on `WeightedMorphologicalTree`:

```python
layout, values = mmcfilters.Attribute.computeAttributes(
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
area = mmcfilters.Attribute.computeSingleAttribute(
    tree,
    mmcfilters.Attribute.AREA,
)

exported_area = tree.project_node_values_to_exported_higra(
    area,
    mmcfilters.Attribute.AREA,
)
```

The projection output follows the same `[proper parts | live internal nodes]`
layout as the exported hierarchy. Internal-node rows are copied from the
node-indexed input. Proper-part rows are filled with unit-component values for
the requested attributes.

Examples of unit proper-part values:

- `AREA`: `1`;
- `LEVEL` and `MEAN_LEVEL`: the altitude of the proper-part owner;
- `VOLUME`: one pixel at the owner altitude;
- `VARIANCE_LEVEL`, `GRAY_HEIGHT`, and `MAX_DIST`: `0`;
- bounding-box attributes: the proper-part row/column coordinates.

Projection fails if:

- the node-value buffer shape does not match the dense internal-node domain;
- the attribute list does not match the number of columns;
- an attribute has no registered unit-component projection;
- the tree is currently inside an edit session;
- a borrowed `WeightedTreeView<T>` became stale after topology mutation.

## Round Trip Pattern

A common interoperability round trip is:

```python
tree = mmcfilters.MorphologicalTreeFactory.createMaxTree(image, radius=1.5)

parent, altitude = tree.exportHigraHierarchy()

roundtrip = mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
    parent,
    altitude,
    image.shape[0],
    image.shape[1],
    mmcfilters.MorphologicalTreeKind.MAX_TREE,
    radius=1.5,
)

area_exported = tree.project_node_values_to_exported_higra(
    mmcfilters.Attribute.computeSingleAttribute(tree, mmcfilters.Attribute.AREA),
    mmcfilters.Attribute.AREA,
)

area_imported_space = mmcfilters.Attribute.computeSingleAttribute(
    roundtrip,
    mmcfilters.Attribute.AREA,
    mmcfilters.NodeIdSpace.HIGRA,
)
```

`area_exported` and `area_imported_space` use the same compact id layout after
the round trip. Both paths fill proper-part rows with unit-component values.

## Choosing The API

Use preserved `NodeIdSpace::HIGRA` when:

- the tree was imported from Higra;
- the topology has not been edited;
- the downstream consumer wants the original imported node ids;
- proper-part rows should follow the same unit-component convention as export.

Use exported Higra projection when:

- the tree was image-built;
- the topology may have been edited;
- the downstream consumer will use `exportHigraHierarchy()` output.

Use internal `NodeIdSpace::MORPHOLOGICAL_TREE` when:

- the data stays inside `mmcfilters`;
- filters, UAO, contours, or topology queries will consume it;
- dead internal slots and live-node iteration are part of the workflow.

## Related Guides

- [Morphological Trees](trees.md): owner/view boundary, `NodeId`, proper parts,
  altitude, and edits.
- [Attributes](attributes.md): attribute layouts, output spaces, and unit
  export projections.
- [Python API Guide](python-api.md): Python construction and wrapper names.
- [Editing API](editing-api.md): edit-session lifetime and derived-state
  invalidation.
