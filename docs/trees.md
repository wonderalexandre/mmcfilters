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
- image-domain shape metadata;
- optional adjacency metadata.

It does not own node altitudes. Any operation that reads or mutates altitude
must use `WeightedMorphologicalTree<T>` or a read-only `WeightedTreeView<T>`.

## Tree Kinds

The implemented high-level tree kinds are:

- `MAX_TREE`: component tree whose altitude is non-decreasing from parent to
  child.
- `MIN_TREE`: component tree whose altitude is non-increasing from parent to
  child.
- `TREE_OF_SHAPES`: tree of shapes built from the 8-bit ToS construction path;
  monotone altitude validation is not applied globally.
- `SELF_DUAL_RESIDUAL_TREE`: self-dual residual tree built from paired min/max
  component-tree seeds.

Component trees require an image-domain adjacency relation. On the 2D square
grid, `radius=1.0` corresponds to 4-connectivity and `radius=1.5` corresponds
to 8-connectivity.

Some attributes are tree-kind dependent. For example, `MAX_DIST` is currently
defined for max-trees and min-trees with adjacency metadata and rejects
tree-of-shapes and self-dual residual trees explicitly.

## Construction

Use `MorphologicalTreeFactory` for public construction:

```cpp
#include <mmcfilters/trees/MorphologicalTreeFactory.hpp>
#include <mmcfilters/utils/Image.hpp>

using namespace mmcfilters;

auto image = ImageUInt8::create(4, 4, std::uint8_t{0});

auto maxTree = MorphologicalTreeFactory::createMaxTree(image, 1.5);
auto minTree = MorphologicalTreeFactory::createMinTree(image, 1.5);
auto tos = MorphologicalTreeFactory::createTreeOfShapes(image);
```

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

Max-trees and min-trees validate monotone altitude order. Trees of shapes and
self-dual residual trees currently skip global monotone-altitude validation
because their semantics are not a single max/min component-tree order. Operators
that need classical altitude polarity, such as MSER, therefore reject those tree
kinds; use topological depth-stability operators when `delta` should mean a
number of tree edges instead of an altitude difference. `ExtinctionValues` is
also restricted to max-trees and min-trees because it currently identifies
regional extrema with component-tree leaves.

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
validate altitude order on commit.

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

Python does not expose `TreeEditor`, `WeightedTreeEditor`, unchecked altitude
setters, or a mutable topology handle from `WeightedMorphologicalTree`.

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
