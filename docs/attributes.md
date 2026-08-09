# Attributes

This document describes the public attribute-computation API and the result
contracts for node-indexed attribute buffers.

## Purpose

The attribute layer turns public attribute requests into dense per-node buffers.
It handles:

- scalar attributes and attribute groups;
- topology/support and altitude-aware computation;
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

For the attribute list, see [Attribute catalog](attribute-catalog.md). For the
tree ownership, altitude, and `NodeId` model that attribute buffers use, see
[Morphological trees](trees.md). For reconstruction operators that consume
node-indexed attribute buffers, see
[Filters and hierarchy operators](filters.md).

## Tree contracts

Topology/support requests may run on `MorphologicalTree`. Requests that read
altitude must use `WeightedMorphologicalTree<T>` or `WeightedTreeView<T>`.

The attribute layer relies on the tree contracts documented in
[Morphological trees](trees.md):

- `WeightedMorphologicalTree<T>` owns topology plus a dense altitude buffer.
- `WeightedTreeView<T>` borrows topology plus an external altitude span.
- public computation checks that the tree is not inside an editing session;
- altitude buffers must match the number of internal node slots;
- attributes requiring adjacency must check that adjacency metadata is present.

Individual attributes may add stricter checks. For example, `MAX_DIST` requires
an altitude buffer, a regular 2D domain, uniform adjacency, and a globally
monotone altitude order. Its descriptive tree kind is irrelevant.

## Common C++ usage patterns

Choose the public entry point from the input contract:

- use weighted computation when any requested attribute may read altitude;
- use topology/support computation only when the request does not read altitude;
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

The `Real` argument selects the public result storage. The typed attribute
facade computes through the same internal `double` pipeline and casts only when
materializing the returned buffer. Integer-valued support attributes are still
counted discretely and then materialized in the requested real type.

Several scalar attributes or groups in one coordinated request:

```cpp
auto [names, values] = AttributeComputation::computeAttributes(
    weightedTree,
    std::vector<AttributeOrGroup>{AREA, LEVEL, AttributeGroup::GRAY_LEVEL});
```

Topology/support attributes without requiring an altitude-bearing weighted tree:

```cpp
auto [names, values] = AttributeComputation::computeTopologyAttributes(
    tree,
    std::vector<AttributeOrGroup>{AREA, BOX_WIDTH, BALANCE_NODE});
```

Returning values in a preserved public node ID space, for a tree created by
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

## Python surface

Python keeps a smaller public surface than C++:

- image factories currently expose the canonical `np.uint8`
  `WeightedMorphologicalTree` path;
- `Attribute.computeSingleAttribute(...)` and
  `Attribute.computeAttributes(...)` are the weighted attribute entry points;
- `Attribute.computeSingleTopologyAttribute(...)` and
  `Attribute.computeTopologyAttributes(...)` are the explicit
  topology/support entry points;
- `NodeIdSpace` can be passed to the Python attribute methods when a preserved
  output node ID space is needed;
- attribute methods accept `dtype=np.float32` or `dtype=np.float64`; the default
  remains `np.float32`;
- `AttributePipeline`, concrete C++ computers, and local-event storage are not
  part of the Python API.

Typical Python calls use the same weighted versus topology/support split:

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

The `dtype` keyword selects the returned NumPy storage. Both supported dtypes
use the same internal `double` computation pipeline. Filtering helpers accept
either dtype when the array satisfies their one-dimensional, contiguous buffer
contract. See [Python API](python-api.md) for Python-specific examples and
failure modes.

## Result layout and output spaces

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
different public node ID space, but projection always happens after the
internal pipeline has computed the result in `NodeIdSpace::MORPHOLOGICAL_TREE`.

`NodeIdSpace::HIGRA` means the preserved imported Higra node ID domain. It is
available only for trees imported from Higra whose original node ID space has
not been invalidated by edits. Direct projection to this space copies live
internal-node rows and fills proper-part rows with unit-component values for
the requested attributes.

For a compact Higra layout exported from the current tree, use
`AttributeComputation::projectNodeValuesToExportedHigra(...)`. That
helper emits the `[proper parts | live internal nodes]` layout produced by
hierarchy export and fills unit proper-part rows through the responsible
attribute computers. `computeAttributeMapping(...)` is the image-domain helper:
each proper part receives the value stored at its proper-part owner.

For the distinction between preserved imported Higra node IDs and exported
compact Higra snapshots, see [Higra interoperability](higra-interoperability.md).

## Numeric stability contract

Scalar attributes returned by the ordinary public attribute APIs are finite for
valid live nodes and exported proper-part rows. This includes degenerate
supports such as one-pixel components, line-like components, zero continuous
bitquad perimeter, and bitquad configurations whose Euler estimate is zero.

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

## Related guides

- [Attribute catalog](attribute-catalog.md): public scalar attributes, groups,
  and input contracts.
- [Attribute computer architecture](attribute-computer-architecture.md):
  contributor guide for adding or changing attribute computers.
