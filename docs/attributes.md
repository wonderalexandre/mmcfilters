# Attributes

This document describes the public attribute-computation API and the result
contracts for node-indexed attribute buffers.

## Purpose

The attribute layer turns public attribute requests into dense per-node buffers.
It handles:

- scalar attributes and attribute groups;
- topology-only and altitude-aware computation;
- shared intermediate quantities such as area and volume;
- optional projection from internal `NodeId` space to public export layouts.

The canonical execution layout is always the dense internal `NodeId` space of
`MorphologicalTree`. Projection is a boundary operation, not the internal
representation used by computers.

The subsystem is designed for incremental, tree-structured attribute
computation: computers accumulate family-specific state over the current
morphological tree and materialize dense public buffers at API boundaries.
After topology edits, recompute public attribute buffers unless a higher-level
operator documents its own edit-aware update path.

## Public API

For ordinary C++ use, include:

```cpp
#include <mmcfilters/attributes/Attributes.hpp>
```

Main entry points:

- `AttributeComputation`: public facade for computing one attribute,
  one group, or a heterogeneous request.
- `AttributeNames`: describes column layout in flat attribute buffers.
- `ComputedAttributeData<Real>` and `ComputedAttributeDataWithDelta<Real>`:
  owning result types returned by the facade.

Use `AttributeComputation` for normal application code. Concrete computers are
advanced extension components, not an alternate public orchestration path.

For the descriptor list, see [Attribute Catalog](attribute-catalog.md). For the
tree ownership, altitude, and `NodeId` model that attribute buffers use, see
[Morphological Trees](trees.md). For reconstruction operators that consume
node-indexed attribute buffers, see
[Attribute Filters, Extinction Values, And UAO](filters.md).

## Tree Contracts

Topology-only requests may run on `MorphologicalTree`. Requests that read
altitude must use `WeightedMorphologicalTree<T>` or `WeightedTreeView<T>`.

The attribute layer relies on the tree contracts documented in
[Morphological Trees](trees.md):

- `WeightedMorphologicalTree<T>` owns topology plus a dense altitude buffer.
- `WeightedTreeView<T>` borrows topology plus an external altitude span.
- public computation checks that the tree is not inside an editing session;
- altitude buffers must match the number of internal node slots;
- attributes requiring adjacency must check that adjacency metadata is present.

Individual attributes may add stricter checks. For example, `MAX_DIST` requires
max-tree or min-tree topology with adjacency metadata and rejects unsupported
tree kinds explicitly.

## Common C++ Usage Patterns

Choose the public entry point from the input contract:

- use weighted computation when any requested attribute may read altitude;
- use topology computation only when the request is explicitly support/topology
  only;
- choose non-default output spaces only at API boundaries.

Single altitude-aware attribute:

```cpp
auto [names, values] =
    AttributeComputation::computeSingleAttribute(weightedTree, LEVEL);
```

Attribute value buffers default to `float`, but the public facade can also
materialize `double` buffers by selecting the `Real` template argument:

```cpp
auto single32 = AttributeComputation::computeSingleAttribute(weightedTree, LEVEL);
auto single64 = AttributeComputation::computeSingleAttribute<double>(weightedTree, LEVEL);
auto mapped64 = AttributeComputation::computeAttributeMapping<double>(weightedTree, LEVEL);
```

The `Real` argument selects both the public result storage and the internal
floating-point arithmetic used by the typed attribute facade. Integer-valued
support descriptors are still counted discretely and then materialized in the
requested real type.

Several scalar attributes or groups in one coordinated request:

```cpp
auto [names, values] = AttributeComputation::computeAttributes(
    weightedTree,
    std::vector<AttributeOrGroup>{AREA, LEVEL, AttributeGroup::GRAY_LEVEL});
```

Topology-only attributes without requiring an altitude-bearing owner:

```cpp
auto [names, values] = AttributeComputation::computeTopologyAttributes(
    tree,
    std::vector<AttributeOrGroup>{AREA, BOX_WIDTH, BALANCE_NODE});
```

Returning values in a preserved public node-id space, for a tree created by
`createFromHigraParent(...)` and not edited since import:

```cpp
auto [names, values] = AttributeComputation::computeAttributes(
    importedTree,
    std::vector<AttributeOrGroup>{AREA, LEVEL},
    NodeIdSpace::HIGRA);
```

Delta-augmented sampling around one scalar attribute:

```cpp
auto [names, values] =
    AttributeComputation::computeSingleAttributeWithDelta(
        weightedTree,
        LEVEL,
        AltitudeDiff<std::uint8_t>{1},
        2);
```

Projecting node attributes to pixels or to an exported Higra layout:

```cpp
auto mapped = AttributeComputation::computeAttributeMapping(weightedTree, AREA);

auto internal = AttributeComputation::computeAttributes(
    weightedTree,
    std::vector<AttributeOrGroup>{AREA, LEVEL});
auto exported = AttributeComputation::projectNodeValuesToExportedHigra(
    weightedTree,
    internal.attributeNames(),
    internal.values());
```

## Python Surface

Python keeps a smaller public surface than C++:

- image factories currently expose the canonical `np.uint8`
  `WeightedMorphologicalTree` path;
- `Attribute.computeSingleAttribute(...)` and
  `Attribute.computeAttributes(...)` are the weighted attribute entry points;
- `Attribute.computeSingleTopologyAttribute(...)` and
  `Attribute.computeTopologyAttributes(...)` are the explicit
  topology/support-only entry points;
- `NodeIdSpace` can be passed to the Python attribute methods when a preserved
  output node-id space is needed;
- attribute methods accept `dtype=np.float32` or `dtype=np.float64`; the default
  remains `np.float32`;
- `AttributePipeline`, concrete C++ computers, local-event deltas, and bitquad
  delta buffers are not Python API.

Typical Python calls use the same weighted versus topology-only split:

```python
level_by_node = mmcfilters.Attribute.computeSingleAttribute(
    weighted_tree,
    mmcfilters.Attribute.LEVEL,
)
names, values = mmcfilters.Attribute.computeAttributes(
    weighted_tree,
    [mmcfilters.Attribute.AREA, mmcfilters.Attribute.Group.GRAY_LEVEL],
)
topology_names, topology_values = mmcfilters.Attribute.computeTopologyAttributes(
    weighted_tree,
    [mmcfilters.Attribute.AREA, mmcfilters.Attribute.Group.BOUNDARY],
)
```

The `dtype` keyword selects the public result storage and maps directly to the
C++ attribute facade: `np.float32` uses `float` and `np.float64` uses `double`.

```python
area32 = mmcfilters.Attribute.computeSingleAttribute(
    weighted_tree,
    mmcfilters.Attribute.AREA,
)
area64 = mmcfilters.Attribute.computeSingleAttribute(
    weighted_tree,
    mmcfilters.Attribute.AREA,
    dtype=np.float64,
)
names, values64 = mmcfilters.Attribute.computeAttributes(
    weighted_tree,
    [mmcfilters.Attribute.AREA, mmcfilters.Attribute.Group.GRAY_LEVEL],
    dtype=np.float64,
)
```

Python filtering helpers such as `AttributeFilters`, `ExtinctionValues`, and
`UltimateAttributeOpening` consume both `np.float32` and `np.float64` attribute
buffers, matching the dtype selected by the attribute computation call.

## Result Layout And Output Spaces

Attribute results are dense flat buffers interpreted by `AttributeNames`.
The canonical internal layout is:

```text
values[node_id * num_attributes + attribute_column]
```

where `node_id` is an internal dense `MorphologicalTree` node slot. Dead
internal slots keep the default buffer value; consumers that reason about tree
nodes should iterate `tree.getAliveNodeIds()`.

`ComputedAttributeData<Real>` and `ComputedAttributeDataWithDelta<Real>` also
store the `NodeIdSpace` of the returned buffer. The default public computation
type is `ComputedAttributeData<float>` or
`ComputedAttributeDataWithDelta<float>`; selecting `Real=double` returns the
corresponding `double` specialization. Public computation methods can request a
different public node-id space, but projection always happens after the
internal pipeline has computed the result in `NodeIdSpace::MORPHOLOGICAL_TREE`.

`NodeIdSpace::HIGRA` means the preserved imported Higra node-id domain. It is
available only for trees imported from Higra whose original node-id space has
not been invalidated by edits. Direct projection to this space copies live
internal-node rows and fills proper-part rows with unit-component values for
the requested attributes.

For a compact Higra layout exported from the current tree, use
`AttributeComputation::projectNodeValuesToExportedHigra(...)`. That
helper emits the `[proper parts | live internal nodes]` layout produced by
hierarchy export and fills unit proper-part rows through the responsible
attribute computers. `computeAttributeMapping(...)` is the image-domain helper:
each proper part receives the value stored at its owner node.

For the distinction between preserved imported Higra ids and exported compact
Higra snapshots, see [Higra Interoperability](higra-interoperability.md).

## Numeric Stability Contract

Scalar attributes returned by the ordinary public attribute APIs are finite for
valid live nodes and exported proper-part rows. This includes degenerate
supports such as one-pixel components, line-like components, zero continuous
BitQuads perimeter, and BitQuads configurations whose Euler estimate is zero.

The finite fallbacks are part of the public attribute contract:

| Attribute family | Degenerate condition | Returned value |
| --- | --- | --- |
| `BITQUADS_CIRCULARITY` | `BITQUADS_PERIMETER_CONTINUOUS <= eps` | `0` |
| `BITQUADS_PERIMETER_AVERAGE` | estimated Euler component count `<= 0` | `0` |
| `BITQUADS_LENGTH_AVERAGE` | estimated Euler component count `<= 0` | `0` |
| `BITQUADS_WIDTH_AVERAGE` | continuous perimeter `<= eps` | `0` |
| `ECCENTRICITY` | both inertia eigenvalues are numerically zero | `1` |
| `ECCENTRICITY` | smallest inertia eigenvalue `<= eps` or ratio overflows the practical range | `1e6` |

`BITQUADS_WIDTH_AVERAGE` is computed with the algebraically equivalent finite
formula `2 * BITQUADS_AREA / BITQUADS_PERIMETER_CONTINUOUS` when the continuous
perimeter is positive. This avoids the `inf / inf` pattern that appears if the
formula is evaluated through average area and average perimeter separately.

The finite-scalar contract does not mean every buffer cell is finite in every
API mode:

- dead internal node slots are addressable and may retain sentinel `NaN` values;
- `computeSingleAttributeWithDelta(..., "nan-padding")` and
  `"null-padding"` deliberately keep `NaN` for missing ancestor or descendant
  samples;
- callers that provide non-finite floating-point altitudes are rejected by
  weighted-tree and altitude-validation code before attributes are computed.

## Related Guides

- [Attribute Catalog](attribute-catalog.md): public scalar descriptors, groups,
  and input contracts.
- [Attribute Computer Architecture](attribute-computer-architecture.md):
  contributor guide for adding or changing attribute computers.
