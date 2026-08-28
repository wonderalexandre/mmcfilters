# Morphological trees

This guide describes the public hierarchy model, construction factories,
smallest-node rules, node domains, and altitude contracts used by `mmcfilters`.

## Public model

The public tree API separates topology from node altitude:

- `MorphologicalTree` owns rooted topology and smallest-node mapping;
- `ValuedMorphologicalTree<T>` owns a topology and a dense altitude buffer;
- `ValuedMorphologicalTreeView<T>` borrows a topology and an external altitude span.

`MorphologicalTree` stores:

- parent/child relations and one live root;
- dense internal `NodeId` slots;
- one smallest node for every pixel;
- optional regular 2D geometry;
- optional regular-grid adjacency semantics.

Pixel-domain objects use `PixelId`: this includes the infinity pixel,
proper-part elements, support and boundary pixels, and flat-zone
representatives. Tree slots and parent/child links use `NodeId`. Their current
integer representations are equal, but code must use `InvalidPixel` for the
pixel domain and `InvalidNode` for the node domain. Counts use ordinary
integral count types rather than either identifier alias.

Operations that read altitude require `ValuedMorphologicalTree<T>` or
`ValuedMorphologicalTreeView<T>`. Topology/support operations can use the topology alone.

## Capabilities and descriptive kinds

`MorphologicalTreeSemantics` declares the capabilities of a hierarchy:

- `NodeAltitudeOrder::Increasing` requires
  `altitude(parent) < altitude(child)`;
- `NodeAltitudeOrder::Decreasing` requires
  `altitude(parent) > altitude(child)`;
- `NodeAltitudeOrder::Unconstrained` declares no global altitude direction;
- `MorphologicalTreeConstructionContext` records explicit absence of
  construction provenance, a shared adjacency, saturated-residual parameters,
  or a complete topographic convention;
- `MorphologicalTreeKind` is a descriptive label, not an algorithm-dispatch
  contract.

Algorithms validate the capabilities they use. For example, `MAX_DIST` requires
an altitude buffer, a regular 2D domain, uniform adjacency, and a globally
monotone altitude order. It does not require a particular descriptive kind.

| Operation | Required capabilities |
| --- | --- |
| topology traversal and `AREA` | finite smallest-node mapping |
| `nodeAltitude`, `nodeAltitudes`, and reconstruction from node altitudes | altitude buffer |
| `GrayLevelHeight`, `MeanGrayLevel`, `GrayLevelVariance`, and `VOLUME` | altitude buffer |
| reconstruction and pixel mapping | regular 2D domain |
| moments, bounding boxes, and contours | regular 2D domain |
| bitquad attributes | regular 2D domain, canonical 4/8 projection connectivity, and exact node altitudes when lower/upper shape connectivity differs |
| `MAX_DIST` | regular 2D domain, uniform adjacency, and monotone altitude order |
| hierarchy saliency projection | regular 2D domain and one compatible projection adjacency |

### Regular-grid adjacency

On a 2D square grid, `radius=1.0` selects 4-connectivity and `radius=1.5`
selects 8-connectivity. `RegularGridAdjacency2D` also supports centered
rectangles, digital lines, and symmetric structuring elements:

```cpp
const RegularGridAdjacency2D disk(rows, columns, 1.5);
auto rectangle = RegularGridAdjacency2D::rectangular(rows, columns, 1, 2);
auto horizontal = RegularGridAdjacency2D::horizontalLine(rows, columns, 3);

std::array<GridOffset2D, 5> cross{{
    {0, 0}, {-1, 0}, {0, -1}, {0, 1}, {1, 0}
}};
auto custom = RegularGridAdjacency2D::fromStructuringElement(
    rows,
    columns,
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
auto residual = MorphologicalTreeFactory::createUnrestrictedResidualTree(
    image,
    1.5);
auto saturatedResidual =
    MorphologicalTreeFactory::createSaturatedResidualTree(
        image,
        PixelId{0},
        1.5);
auto treeOfShapes = MorphologicalTreeFactory::createTreeOfShapes(image);
```

C++ max-tree and min-tree factories are typed:

```cpp
ImagePtr<float> floatImage = Image<float>::create(4, 4, 0.0f);
ValuedMorphologicalTree<float> floatMaxTree =
    MorphologicalTreeFactory::createMaxTree(floatImage, 1.5);
```

Python image factories accept two-dimensional C-contiguous `np.uint8` arrays
and return `ValuedMorphologicalTree`. See [Python API](python-api.md) for the
binding-specific contract.

### Self-dual residual trees

`createUnrestrictedResidualTree` constructs a hierarchy from synchronized max-tree
and min-tree states using one symmetric adjacency. The saturated variant also
requires the complement of a selected extremum to remain connected to the
row-major `infinityPixel` supplied by the caller.

The unrestricted and saturated factories publish distinct
`UnrestrictedResidualTree` and `SaturatedResidualTree` kinds. The former
retains a `SharedAdjacencyContext`; the latter retains a
`SaturatedResidualContext`, including its infinity pixel. Both publish:

- `NodeAltitudeOrder::Unconstrained`.

Their alternating gray-level altitudes are not a monotone hierarchy valuation.
Use a structural or otherwise non-decreasing node valuation for operations that
require values to increase toward the root.

An explicit adjacency can replace the radius:

```cpp
auto adjacency = RegularGridAdjacency2D::fromStructuringElement(
    rows,
    columns,
    offsets);
auto residual = MorphologicalTreeFactory::createUnrestrictedResidualTree(
    image,
    adjacency);
```

Construction policies are grouped by mode, so unrestricted callers cannot
accidentally configure saturation-only mechanisms:

```cpp
#include <mmcfilters/trees/sdrt/ResidualTreePolicies.hpp>

sdrt::UnrestrictedResidualTreeOptions unrestrictedOptions;
unrestrictedOptions.spatialOrder =
    sdrt::SpatialOrder(std::vector<PixelId>{3, 2, 1, 0});
auto unrestricted = MorphologicalTreeFactory::createUnrestrictedResidualTree(
    image,
    adjacency,
    unrestrictedOptions);

sdrt::SaturatedResidualTreeOptions saturatedOptions;
saturatedOptions.lcaPolicy = sdrt::SaturatedMinMaxLcaPolicy::LinkCut;
saturatedOptions.fallbackPolicy =
    sdrt::SaturatedMinMaxFallbackPolicy::BoundaryMultiSource;
auto saturated =
    MorphologicalTreeFactory::createSaturatedResidualTree(
        image,
        adjacency,
        PixelId{0},
        saturatedOptions);
```

By default, construction uses `RowMajorSpatialOrder`, parent-climb LCA, and
multi-source boundary fallback. The unique self-dual schedule orders candidates
by `(supportCardinality, spatialMinimum)` and never by polarity. Both modes use
incremental small-to-large flat-zone boundary maintenance and
region-contraction event assembly. Configuration is supplied through
mode-specific C++ option objects.

### Tree-of-shapes policies

`TopographicConvention` records and controls:

- `CanonicalComplementaryGridImmersion` carrying only its `ComplementaryPairing`,
  a `ComplementaryGridImmersion` carrying explicit minimum/maximum adjacencies,
  or `SelfDualSpanImmersion`;
- `TopographicDomainExtension::ExteriorRing` or `None`;
- the row-major `infinityPixel` in the active topographic domain;
- the published altitude scale, `TopographicAltitudeEncoding::UInt8` or
  `ExactDoubled`.

The default convention selects the canonical minimum-4/maximum-8
complementary-grid immersion without domain padding, uses infinity pixel zero,
and publishes unchanged 8-bit source levels:

```cpp
ValuedMorphologicalTree<std::uint8_t> treeOfShapes =
    MorphologicalTreeFactory::createTreeOfShapes(image);
```

A canonical immersion declares only the ordered pairing, so it needs no domain.
It is resolved against the source domain during construction, and the retained
convention always exposes the concrete adjacencies:

```cpp
TopographicConvention convention{
    CanonicalComplementaryGridImmersion{ComplementaryPairing::Min8Max4},
    TopographicDomainExtension::None,
    PixelId{0}};
auto unpadded = MorphologicalTreeFactory::createTreeOfShapes(image, convention);
```

Both padding policies publish the original source pixels as the tree domain and
attach the original `GridDomain2D`. The immersion grid is not exposed as the
pixel domain.

### Tree-of-shapes altitude encodings

A complementary-grid immersion floods the interpolated domain over the source
level set itself, so every construction level is already a source gray level.
`TopographicAltitudeEncoding::UInt8` publishes those levels unchanged and
`reconstructFromNodeAltitudes()` returns the source image. The encoding is exact,
not a quantization: the published hierarchy is the doubled hierarchy with every
altitude halved, so no parent-child altitude distinction is lost.

A self-dual span immersion propagates from a boundary reference level that, on
an even boundary, is the mean of the two central boundary values and may
therefore fall between two source levels. That reference level lives on the
exterior ring and is the only source of half levels, so the admissible encodings
depend on the domain extension:

| Immersion | `TopographicDomainExtension` | Admissible encodings |
| --- | --- | --- |
| Complementary grid | either | `UInt8`, `ExactDoubled` |
| Self-dual span | `None` | `UInt8`, `ExactDoubled` |
| Self-dual span | `ExteriorRing` | `ExactDoubled` only |

Without the exterior ring the reference level is cropped away and no interior
cell reads it, so every construction level stays on the source lattice and the
8-bit encoding is exact:

```cpp
ValuedMorphologicalTree<std::uint8_t> unpaddedSelfDual =
    MorphologicalTreeFactory::createTreeOfShapes(
        image,
        TopographicConvention{SelfDualSpanImmersion{},
                              TopographicDomainExtension::None,
                              PixelId{0},
                              TopographicAltitudeEncoding::UInt8});

ValuedMorphologicalTree<ToSGrayLevel> selfDual =
    MorphologicalTreeFactory::createTreeOfShapes<ToSGrayLevel>(
        image, selfDualSpanConvention());
```

Under `ExactDoubled` a source value `v` is represented by `2 * v`; odd values
represent half levels, and `reconstructFromNodeAltitudes()` also returns doubled
values. Half levels are never quantized into `std::uint8_t`, because doing so can
create equal-altitude parent-child edges and destroy shape polarity.

The published altitude type is a compile-time choice in C++ and must match the
encoding declared by the convention; `createTreeOfShapes` throws
`std::invalid_argument` when they disagree. Python selects the instantiation from
the `altitude_encoding` field at run time.

Every non-root tree-of-shapes node has a polarity derived from these exact
altitudes:

- `ShapePolarity::Upper` when the node altitude is greater than its parent altitude;
- `ShapePolarity::Lower` when the node altitude is less than its parent altitude.

The root has no shape polarity. Scalar bitquad materialization uses an explicit
`BitquadConnectivityPolicy` after hierarchy-independent bitquad-family counts
have been computed. Under a complementary convention, lower shapes use
`minAdjacency`, upper shapes use `maxAdjacency`, and the unpolarized root follows
the policy's explicit root connectivity. The current root projection remains
8-connected when the two complementary choices differ; this is not inferred by
classifying the root as an upper or lower shape.

## Owning tree and view boundary

Use `ValuedMorphologicalTree<T>` when an operation must own topology and
altitude state, including:

- construction and import;
- topology or altitude mutation;
- image reconstruction;
- Higra export and projection;
- filters, extinction values, `UltimateAttributeOpening`,
  `adjust::CasfComponentTrees`, and paired-tree adjustment;
- Python workflows.

Use `ValuedMorphologicalTreeView<T>` for read-only C++ operations over caller-owned altitude
storage:

```cpp
const MorphologicalTree& topology = valuedTree.topology();
std::vector<float> altitude(topology.numInternalNodeSlots(), 0.0f);

ValuedMorphologicalTreeView<float> view(topology, altitude);
auto rootAltitude = view.nodeAltitude(topology.root());
```

The caller must keep the topology and altitude storage alive for the lifetime of
the view. A view captures the topology mutation version and rejects reads after
the topology changes.

## Nodes, supports, proper parts, and pixels

The internal node domain is dense:

```text
0 <= node_id < tree.numInternalNodeSlots()
```

Edits may leave dead internal slots. Node-indexed buffers retain the complete
slot count so existing IDs remain addressable. In C++, iterate
`tree.aliveNodeIds()` when consuming live nodes; in Python, use the
`tree.alive_node_ids` property.

Image pixels have row-major IDs:

```text
pixel = row * numColumns + column
```

The main structural queries are:

```cpp
NodeId smallest = tree.smallestNode(pixel);
auto proper = tree.properPart(node);
auto support = tree.nodeSupport(node);
```

`properPart(node)` returns the pixels in
\f$\rho(X)=X\setminus\bigcup_{Y\in children(X)}Y\f$.
`nodeSupport(node)` returns every pixel in the complete node support. The
smallest-node map satisfies `pixel in properPart(smallestNode(pixel))`.

### Topology queries and traversal contracts

Topology query names describe the relation they return:

| Meaning | C++ | Python |
| --- | --- | --- |
| root node | `root()` | `root` |
| parent node | `parent(node)` | `parent(node)` |
| direct children | `children(node)` | `children(node)` |
| inclusive ancestor chain | `ancestors(node)` | `ancestors(node)` |
| strict descendants | `descendants(node)` | `descendants(node)` |
| inclusive subtree | `subtreeNodes(node)` | `subtree_nodes(node)` |
| lowest common ancestor | `lowestCommonAncestor(a, b)` | `lowest_common_ancestor(a, b)` |

`ancestors(node)` includes `node` itself and ends at the root. Conversely,
`descendants(node)` excludes `node`, whereas `subtreeNodes(node)` includes it.
The subtree range uses pre-order, but callers must rely on the stated relation
rather than on incidental sibling order unless an algorithm explicitly defines
that order.

`postOrder()` and `breadthFirstTraversal()` (Python: `post_order` and
`breadth_first_traversal`) return tree-wide traversal schedules. Post-order
visits every child before its parent. It can schedule bottom-up aggregation,
but traversal and aggregation are distinct operations.

`dfsEntryIndex(node)` and `dfsExitIndex(node)` (Python:
`dfs_entry_index(node)` and `dfs_exit_index(node)`) are indices in one shared,
zero-based DFS event sequence. The counter advances when DFS enters a node and
again when it exits after visiting every descendant. For a tree with
`numNodes()` live nodes, the combined entry and exit indices form a permutation
of `[0, 2 * numNodes())`.

The inclusive ancestor relation is characterized by interval containment:
$u$ is an ancestor of $v$ exactly when
`dfsEntryIndex(u) <= dfsEntryIndex(v)` and
`dfsExitIndex(u) >= dfsExitIndex(v)`. The number of strict descendants is
`(dfsExitIndex(node) - dfsEntryIndex(node) - 1) / 2`. These values are event
indices, not wall-clock timestamps or independent positions in two traversals.

Node IDs and pixel IDs are independent domains. A live node may have an empty
proper part when its support is supplied by descendants. Every committed live
node must nevertheless have non-empty node support. Such a tree remains a valid
`MorphologicalTree`, but it is not a morphological tree of partial partitions;
use `isTreeOfPartialPartitions()` or `validateTreeOfPartialPartitions()` when
that stricter specialization is required.

### Native topology import

Builders that already have node parents and smallest-node mapping can use the
native representation directly:

```cpp
std::vector<NodeId> nodeParent{0, 0, 0};
std::vector<NodeId> smallestNodeMap{1, 2};
std::vector<std::uint8_t> altitude{0, 1, 0};

auto tree = MorphologicalTreeFactory::createFromNativeTopology(
    std::span<const NodeId>(nodeParent),
    std::span<const NodeId>(smallestNodeMap),
    std::span<const std::uint8_t>(altitude),
    0,
    MorphologicalTreeSemantics{
        MorphologicalTreeKind::Generic,
        NodeAltitudeOrder::Unconstrained,
        NoAdjacency{}});
```

This overload attaches no coordinate interpretation. Topology/support
attributes remain available, while reconstruction and geometry-dependent
attributes reject the missing `GridDomain2D`. An overload with `rows` and
`columns` attaches a regular 2D domain and requires one smallest-node-map entry
per grid position.

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
ordered hierarchies. `NodeAltitudeOrder::Unconstrained` skips only monotonicity validation; buffer
shape and finite values remain required.

`GrayLevelHeight` is defined for every altitude order as the maximum absolute
altitude difference between a node and any node in its subtree. On monotone
max-tree and min-tree hierarchies this is the usual one-sided gray-level span.

### Residue and reconstruction baseline

The project adopts the fixed reconstruction baseline zero. For every non-root
node `n`, `nodeResidue(n)` is
`nodeAltitude(n) - nodeAltitude(parent(n))`; for the root it is simply
`nodeAltitude(root)`. No baseline parameter or baseline state is part of the
public API.

With this convention, summing node residues along the path from the root to a
node telescopes to that node's altitude. Assigning the resulting altitude to
each pixel through its smallest node therefore reconstructs the image represented
by the valued tree. `reconstructFromNodeAltitudes()` performs the equivalent
direct reconstruction from the stored node altitudes.

## Edits and derived state

`pruneNode(node)` and `mergeNodeIntoParent(node)` are complete local edits.
Multi-step structural changes use an editor session and cross a checked commit
boundary. Topology mutations invalidate cached objects and node-indexed results
computed against the previous topology. See [Editing API](editing-api.md) for
rollback, valued-tree edits, and Python exposure.

## Related guides

- [Attributes](attributes.md): topology/support and altitude-aware computation.
- [Filters](filters.md): reconstruction from dense node buffers.
- [Saliency maps](saliency.md): hierarchy and shape-space projections.
- [Higra interoperability](higra-interoperability.md): imported and exported
  node domains.
- [Editing API](editing-api.md): committed and staged mutation.
- [Python API](python-api.md): Python construction and buffer contracts.
