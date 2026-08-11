# Morphological trees

This guide describes the public hierarchy model, construction factories,
ownership rules, node domains, and altitude contracts used by `mmcfilters`.

## Public model

The public tree API separates topology from node altitude:

- `MorphologicalTree` owns rooted topology and proper-part ownership;
- `WeightedMorphologicalTree<T>` owns a topology and a dense altitude buffer;
- `WeightedTreeView<T>` borrows a topology and an external altitude span.

`MorphologicalTree` stores:

- parent/child relations and one live root;
- dense internal `NodeId` slots;
- one direct proper-part owner for every proper part;
- optional regular 2D geometry;
- optional regular-grid adjacency semantics.

Operations that read altitude require `WeightedMorphologicalTree<T>` or
`WeightedTreeView<T>`. Topology/support operations can use the topology alone.

## Capabilities and descriptive kinds

`HierarchySemantics` declares the capabilities of a hierarchy:

- `AltitudeOrder::INCREASING_FROM_ROOT` requires
  `altitude(parent) < altitude(child)`;
- `AltitudeOrder::DECREASING_FROM_ROOT` requires
  `altitude(parent) > altitude(child)`;
- `AltitudeOrder::UNCONSTRAINED` declares no global altitude direction;
- `AdjacencyMode` declares no adjacency, one uniform adjacency, or distinct
  decreasing/increasing adjacencies;
- `MorphologicalTreeKind` is a descriptive label, not an algorithm-dispatch
  contract.

Algorithms validate the capabilities they use. For example, `MAX_DIST` requires
an altitude buffer, a regular 2D domain, uniform adjacency, and a globally
monotone altitude order. It does not require a particular descriptive kind.

| Operation | Required capabilities |
| --- | --- |
| topology traversal and `AREA` | finite proper-part ownership |
| `LEVEL`, `GRAY_HEIGHT`, and `VOLUME` | altitude buffer |
| reconstruction and pixel mapping | regular 2D domain |
| moments, bounding boxes, and contours | regular 2D domain |
| bitquad attributes | regular 2D domain and canonical 4/8 adjacency when adjacency is needed |
| `MAX_DIST` | regular 2D domain, uniform adjacency, and monotone altitude order |
| hierarchy saliency projection | regular 2D domain and one compatible projection adjacency |

### Regular-grid adjacency

On a 2D square grid, `radius=1.0` selects 4-connectivity and `radius=1.5`
selects 8-connectivity. `RegularGridAdjacency2D` also supports centered
rectangles, digital lines, and symmetric structuring elements:

```cpp
const RegularGridAdjacency2D disk(rows, cols, 1.5);
auto rectangle = RegularGridAdjacency2D::rectangular(rows, cols, 1, 2);
auto horizontal = RegularGridAdjacency2D::horizontalLine(rows, cols, 3);

std::array<GridOffset2D, 5> cross{{
    {0, 0}, {-1, 0}, {0, -1}, {0, 1}, {1, 0}
}};
auto custom = RegularGridAdjacency2D::fromStructuringElement(
    rows,
    cols,
    cross);
```

An adjacency-inducing structuring element must contain the origin exactly once
and be centrally symmetric. A component-tree adjacency must also connect its
finite image domain. Traversal ranges own their cursors, so nested and
interleaved iteration over the same relation is safe.

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
    image,
    1.5);
auto saturatedResidual =
    MorphologicalTreeFactory::createSaturatedSelfDualResidualTree(
        image,
        NodeId{0},
        1.5);
auto treeOfShapes = MorphologicalTreeFactory::createTreeOfShapes(image);
```

C++ max-tree and min-tree factories are typed:

```cpp
ImagePtr<float> floatImage = Image<float>::create(4, 4, 0.0f);
WeightedMorphologicalTree<float> floatMaxTree =
    MorphologicalTreeFactory::createMaxTree(floatImage, 1.5);
```

Python image factories accept two-dimensional C-contiguous `np.uint8` arrays
and return `WeightedMorphologicalTree`. See [Python API](python-api.md) for the
binding-specific contract.

### Self-dual residual trees

`createSelfDualResidualTree` constructs a hierarchy from synchronized max-tree
and min-tree states using one symmetric adjacency. The saturated variant also
requires the complement of a selected extremum to remain connected to the
row-major `infinityPixel` supplied by the caller.

Both factories publish:

- descriptive kind `SELF_DUAL_RESIDUAL_TREE`;
- uniform adjacency;
- `AltitudeOrder::UNCONSTRAINED`.

Their alternating gray-level altitudes are not a monotone hierarchy valuation.
Use a structural or otherwise non-decreasing node valuation for operations that
require values to increase toward the root.

An explicit adjacency can replace the radius:

```cpp
auto adjacency = RegularGridAdjacency2D::fromStructuringElement(
    rows,
    cols,
    offsets);
auto residual = MorphologicalTreeFactory::createSelfDualResidualTree(
    image,
    adjacency);
```

Construction policies are grouped by mode, so unrestricted callers cannot
accidentally configure saturation-only mechanisms:

```cpp
#include <mmcfilters/trees/sdrt/ResidualTreePolicies.hpp>

sdrt::UnrestrictedResidualTreeOptions unrestrictedOptions;
unrestrictedOptions.tiePolicy = sdrt::SdrtTiePolicy::MaxBeforeMinThenSpatial;
auto unrestricted = MorphologicalTreeFactory::createSelfDualResidualTree(
    image,
    adjacency,
    unrestrictedOptions);

sdrt::SaturatedResidualTreeOptions saturatedOptions;
saturatedOptions.lcaPolicy = sdrt::SaturatedMinMaxLcaPolicy::LinkCut;
saturatedOptions.fallbackPolicy =
    sdrt::SaturatedMinMaxFallbackPolicy::BoundaryMultiSource;
auto saturated =
    MorphologicalTreeFactory::createSaturatedSelfDualResidualTree(
        image,
        adjacency,
        NodeId{0},
        saturatedOptions);
```

By default, construction uses contrast-invariant spatial tie-breaking,
parent-climb LCA, and multi-source boundary fallback. Both modes use
incremental small-to-large flat-zone boundary maintenance and
region-contraction event assembly. Configuration is supplied through
mode-specific C++ option objects.

### Tree-of-shapes policies

`TreeOfShapesProducerOptions` controls:

- `SelfDual`, `Min4cMax8c`, or `Min8cMax4c` interpolation/connectivity;
- exterior immersion padding or no padding;
- the \f$p_\infty\f$ seed in the selected immersion domain.

```cpp
auto selfDual = MorphologicalTreeFactory::createTreeOfShapes(
    image,
    ToSInterpolation::SelfDual);

TreeOfShapesProducerOptions options{
    ToSInterpolation::Min4cMax8c,
    ToSPaddingPolicy::NoPadding,
    0,
    0};
auto unpadded = MorphologicalTreeFactory::createTreeOfShapes(image, options);
```

Both padding policies publish exactly one proper part per source pixel and
attach the original `GridDomain2D`. The immersion grid is not exposed as the
proper-part domain.

The `std::uint8_t` overload quantizes half-level structural nodes. C++ callers
that need exact levels can use doubled gray units:

```cpp
WeightedMorphologicalTree<ToSGrayLevel> exact =
    MorphologicalTreeFactory::createTreeOfShapesExact(image, options);
```

A source value `v` is represented by `2 * v`; odd values represent half levels.
`reconstructionImage()` therefore also returns doubled values for an exact
weighted tree. Python exposes the producer options with `np.uint8` altitudes.

## Owning tree and view boundary

Use `WeightedMorphologicalTree<T>` when an operation must own topology and
altitude state, including:

- construction and import;
- topology or altitude mutation;
- image reconstruction;
- Higra export and projection;
- filters, extinction values, `UltimateAttributeOpening`,
  `adjust::CasfComponentTrees`, and paired-tree adjustment;
- Python workflows.

Use `WeightedTreeView<T>` for read-only C++ operations over caller-owned altitude
storage:

```cpp
const MorphologicalTree& topology = weighted.topology();
std::vector<float> altitude(topology.getNumInternalNodeSlots(), 0.0f);

WeightedTreeView<float> view(topology, altitude);
auto rootAltitude = view.getAltitude(topology.getRoot());
```

The caller must keep the topology and altitude storage alive for the lifetime of
the view. A view captures the topology mutation version and rejects reads after
the topology changes.

## Node IDs and proper parts

The internal node domain is dense:

```text
0 <= node_id < tree.getNumInternalNodeSlots()
```

Edits may leave dead internal slots. Node-indexed buffers retain the complete
slot count so existing IDs remain addressable. Iterate `tree.getAliveNodeIds()`
when consuming live nodes.

Image pixels are proper parts with row-major IDs:

```text
proper_part = row * numCols + col
```

The main ownership queries are:

```cpp
NodeId owner = tree.getProperPartOwner(pixel);
auto direct = tree.getProperParts(node);
auto support = tree.getConnectedComponent(node);
```

`getProperParts(node)` returns only the proper parts directly owned by `node`.
`getConnectedComponent(node)` returns the complete subtree support.

Internal nodes and proper parts are independent domains. A live structural node
may own no proper part directly when its support is supplied by descendants.
Every committed live node must nevertheless have non-empty subtree support.

### Native topology import

Builders that already have node parents and proper-part ownership can use the
native representation directly:

```cpp
std::vector<NodeId> nodeParent{0, 0, 0};
std::vector<NodeId> properPartOwner{1, 2};
std::vector<std::uint8_t> altitude{0, 1, 0};

auto tree = MorphologicalTreeFactory::createFromNativeTopology(
    std::span<const NodeId>(nodeParent),
    std::span<const NodeId>(properPartOwner),
    std::span<const std::uint8_t>(altitude),
    0,
    HierarchySemantics{
        MorphologicalTreeKind::GENERIC,
        AltitudeOrder::UNCONSTRAINED,
        NoAdjacency{}});
```

This overload attaches no coordinate interpretation. Topology/support
attributes remain available, while reconstruction and geometry-dependent
attributes reject the missing `GridDomain2D`. An overload with `rows` and
`cols` attaches a regular 2D domain and requires one proper-part owner entry per
grid position.

`NativeHierarchyView<T>` provides the equivalent synchronous C++ import from
borrowed spans. The factory copies or moves the validated buffers into owning
tree storage.

## Altitude contract

Altitude buffers are indexed by internal `NodeId`:

```text
altitude[node_id]
```

Accepted altitude types satisfy `AltitudeValue`:

- floating-point values must be finite;
- integral values must support differences representable in signed 64-bit
  arithmetic;
- `bool` is rejected.

Construction and committed edits validate strict parent/child inequalities for
ordered hierarchies. `UNCONSTRAINED` skips only monotonicity validation; buffer
shape and finite values remain required.

`GRAY_HEIGHT` is defined for every altitude order as the maximum absolute
altitude difference between a node and any node in its subtree. On monotone
max-tree and min-tree hierarchies this is the usual one-sided gray-level span.

## Edits and derived state

`pruneNode(node)` and `mergeNodeIntoParent(node)` are complete local edits.
Multi-step structural changes use an editor session and cross a checked commit
boundary. Topology mutations invalidate cached objects and node-indexed results
computed against the previous topology. See [Editing API](editing-api.md) for
rollback, weighted edits, and Python exposure.

## Related guides

- [Attributes](attributes.md): topology/support and altitude-aware computation.
- [Filters](filters.md): reconstruction from dense node buffers.
- [Saliency maps](saliency.md): hierarchy and shape-space projections.
- [Higra interoperability](higra-interoperability.md): imported and exported
  node domains.
- [Editing API](editing-api.md): committed and staged mutation.
- [Python API](python-api.md): Python construction and buffer contracts.
