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
  construction provenance, a shared adjacency,
  or a complete topographic convention;
- `MorphologicalTreeKind` is a descriptive label, not an algorithm-dispatch
  contract.

Algorithms validate the capabilities they use. A topology/support operation may
require regular 2D geometry without requiring altitude, construction adjacency,
or a particular descriptive tree kind.

`MorphologicalTree::constructionAdjacency()` returns the retained shared
adjacency. `complementaryAdjacencies()` returns the retained minimum/maximum
adjacencies of a complementary grid immersion. Both methods return pointers to
const data in the stored construction context, or `nullptr` when the requested
context is absent.

| Operation | Required capabilities |
| --- | --- |
| topology traversal and `AREA` | finite smallest-node mapping |
| `nodeAltitude`, `nodeAltitudes`, and reconstruction from node altitudes | altitude buffer |
| `GrayLevelHeight`, `MeanGrayLevel`, `GrayLevelVariance`, and `VOLUME` | altitude buffer |
| reconstruction and pixel mapping | regular 2D domain |
| moments, bounding boxes, and contours | regular 2D domain |
| bitquad attributes | regular 2D domain, canonical 4/8 projection connectivity, and exact node altitudes when lower/upper shape connectivity differs |
| distance-transform attributes (`DIST_TRANSF`, `DIST_TRANSF_EXACT`) | regular 2D domain |
| hierarchy saliency projection | regular 2D domain and one compatible projection adjacency |

### Regular-grid adjacency

On a 2D square grid, `radius=1.0` selects 4-connectivity and `radius=1.5`
selects 8-connectivity. Radius construction requires a finite radius of at least
`1.0`. `RegularGridAdjacency2D` also supports centered
rectangles, digital lines, and symmetric structuring elements:

```cpp
const RegularGridAdjacency2D disk(rows, columns, 1.5);
auto rectangle = RegularGridAdjacency2D::rectangular(rows, columns, 1, 2);
auto horizontal = RegularGridAdjacency2D::horizontalLine(1, columns, 3);

std::array<GridOffset2D, 5> cross{{
    {0, 0}, {-1, 0}, {0, -1}, {0, 1}, {1, 0}
}};
auto custom = RegularGridAdjacency2D::fromStructuringElement(
    rows,
    columns,
    cross);
```

An adjacency-inducing structuring element must contain the origin exactly once
and be centrally symmetric. It must connect its complete nonempty finite grid
domain. A horizontal line therefore requires a single-row domain, and a vertical
line requires a single-column domain. Adjacency construction establishes these
properties once in both contract modes; tree builders and factories reuse them.
Traversal ranges own their cursors, so nested and
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
auto extTree = MorphologicalTreeFactory::createExtremumResidualTree(image, 1.5);
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
std::vector<float> altitude(valuedTree.nodeAltitudes().begin(),
                            valuedTree.nodeAltitudes().end());

ValuedMorphologicalTreeView<float> view(topology, altitude);
auto rootAltitude = view.nodeAltitude(topology.root());
```

Public view constructors validate the external buffer shape, finite values, and
strict parent-child inequalities required by the declared altitude order.
`Unconstrained` valuations do not require monotonicity. The caller must keep the
topology and altitude storage alive and unchanged while using the external view.
A view captures the topology mutation version so consumers can reject stale topology.

`valuedTree.asView()` borrows the owned state through a private construction path
without scanning the buffer again. The owner maintains its altitude guarantees
through the altitude setters and the valued-tree editor.

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
| batch of lowest common ancestors | `lowestCommonAncestors(pairs)` | `lowest_common_ancestors(pairs)` |

### Native hierarchy buffers

The public hierarchy import entry point is `createFromHigraParent`. Builders and
internal hierarchy-producing algorithms transfer their native buffers directly
through internal materialization. There is no public native-buffer import API.

Both native hierarchy representations are public types in `mmcfilters`:

- [`NativeHierarchyView<T>`](../mmcfilters/trees/NativeHierarchyView.hpp) groups
  borrowed spans for internal conversion and validation. The buffers must remain
  valid while the view is used. The descriptor is not a public import API.
- [`NativeHierarchyOwned<T>`](../mmcfilters/trees/NativeHierarchyOwned.hpp) owns
  the topology and altitude vectors returned by builders. Copies are disabled;
  moves transfer the buffers and semantic metadata. Materialization consumes
  these buffers without copying them again.

`makeNativeHierarchyOwned` moves existing vectors into the owning representation
without structural validation. Callers must provide valid hierarchy data.

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

Builders establish their declared altitude order by construction; internal
materialization transfers their buffers without repeating validation. Higra imports
validate external structure and altitudes before materialization.
Altitude setters and committed edits preserve the declared order, using local
checks where possible. `NodeAltitudeOrder::Unconstrained` skips only monotonicity
validation; buffer shape and finite values remain required.

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
node telescopes to the altitude of that node. Assigning the resulting altitude to
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
