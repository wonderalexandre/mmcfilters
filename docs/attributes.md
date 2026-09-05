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

Recompute attribute buffers after topology edits. See
[Attribute computer architecture](attribute-computer-architecture.md) for
computation strategies and extension contracts.

## Public API

For ordinary C++ use, include:

```cpp
#include <mmcfilters/attributes/Attributes.hpp>
```

Main entry points:

- `AttributeComputation`: public facade for computing one attribute,
  one group, or a heterogeneous request.
- `AttributeNames`: describes column layout in flat attribute buffers.
- `NodeAttributeSampleLayout`: describes signed sample-offset columns.
- `ComputedAttributeData<Real>` and `SampledNodeAttributeData<Real>`:
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
altitude must use `ValuedMorphologicalTree<T>` or `ValuedMorphologicalTreeView<T>`.

The attribute layer relies on the tree contracts documented in
[Morphological trees](trees.md):

- `ValuedMorphologicalTree<T>` owns topology plus a dense altitude buffer.
- `ValuedMorphologicalTreeView<T>` borrows topology plus an external altitude span.
- public computation checks that the tree is not inside an editing session;
- altitude buffers must match the number of internal node slots;
- attributes requiring adjacency must check that adjacency metadata is present.

Individual attributes may add stricter checks. For example, all distance-field
attributes require a regular 2D domain and use only node supports and their
foreground 4-connected contours. They do not require altitude or
construction-adjacency metadata, and their descriptive tree kind is irrelevant.

Family-specific names, formulas, units, tie rules, and group membership are
documented in the [Attribute catalog](attribute-catalog.md) and the relevant
subsystem guide, such as
[Distance-transform attributes](distance-transform.md).

## Common C++ usage patterns

Choose the public entry point from the input contract:

- use valued-tree computation when any requested attribute may read node altitude;
- use topology/support computation only when the request does not read altitude;
- choose non-default output spaces only at API boundaries.

Single altitude-aware attribute:

```cpp
auto [names, values] =
    AttributeComputation::computeSingleAttribute(valuedTree, GrayLevelHeight);
```

Attribute value buffers default to `float`, but the public facade can also
materialize `double` buffers by selecting the `Real` template argument:

```cpp
auto single32 = AttributeComputation::computeSingleAttribute(valuedTree, GrayLevelHeight);
auto single64 = AttributeComputation::computeSingleAttribute<double>(valuedTree, GrayLevelHeight);
auto mapped64 = AttributeComputation::computeAttributeMapping<double>(valuedTree, GrayLevelHeight);
```

The `Real` argument selects the public result storage. The typed attribute
facade computes through the same internal `double` pipeline and casts only when
materializing the returned buffer. Integer-valued support attributes are still
counted discretely and then materialized in the requested real type.

Several scalar attributes or groups in one coordinated request:

```cpp
auto [names, values] = AttributeComputation::computeAttributes(
    valuedTree,
    std::vector<AttributeOrGroup>{Area, GrayLevelHeight, AttributeGroup::GrayLevel});
```

Topology/support attributes without requiring an altitude-bearing valued tree:

```cpp
auto [names, values] = AttributeComputation::computeTopologyAttributes(
    tree,
    std::vector<AttributeOrGroup>{Area, BoxWidth, BalanceNode});
```

Returning values in a preserved public node ID space, for a tree created by
`createFromHigraParent(...)` and not edited since import:

```cpp
auto [names, values] = AttributeComputation::computeAttributes(
    importedTree,
    std::vector<AttributeOrGroup>{Area, GrayLevelHeight},
    NodeIdSpace::Higra);
```

Altitude-based sampling around one scalar node attribute:

```cpp
auto [names, values] =
    AttributeComputation::computeSampledNodeAttribute(
        valuedTree,
        MeanGrayLevel,
        AltitudeDifference<std::uint8_t>{1},
        2);
```

The altitude step must be positive. Offset `0` is the current node, negative
offsets are ancestor samples, and positive offsets are representative-
descendant samples. The default `LargestSupportDescendant` policy selects the
candidate with greatest node-support cardinality and breaks ties by the
minimum pixel index in the node support. Missing positions use one of
the typed `MissingNodeAttributeSamplePolicy` values: `RepeatNearest` (default),
`NotANumber`, or `Zero`.

Projecting node attributes to pixels or to an exported Higra layout:

```cpp
auto mapped = AttributeComputation::computeAttributeMapping(valuedTree, Area);

auto internal = AttributeComputation::computeAttributes(
    valuedTree,
    std::vector<AttributeOrGroup>{Area, GrayLevelHeight});
auto exported = AttributeComputation::projectNodeValuesToExportedHigra(
    valuedTree,
    internal.attributeNames(),
    internal.values());
```

## Python surface

Python uses the same valued-tree and topology/support distinction:

```python
gray_level_height_by_node = mmcfilters.Attribute.compute_single_attribute(
    valued_tree,
    mmcfilters.Attribute.GRAY_LEVEL_HEIGHT,
)
area_by_node = mmcfilters.Attribute.compute_single_topology_attribute(
    valued_tree,
    mmcfilters.Attribute.AREA,
)
```

See [Python API](python-api.md) for NumPy layouts, dtypes, symbolic names,
sampling, and buffer requirements.

## Result layout and output spaces

C++ attribute results are dense flat buffers interpreted by `AttributeNames`.
Ordinary multi-attribute results use node-major layout, with columns in
scalar-enumerator order:

```text
values[node_id * num_attributes + attribute_column]
```

where `node_id` is an internal dense `MorphologicalTree` node slot. Dead
internal slots keep the default buffer value; consumers that reason about tree
nodes should iterate `tree.aliveNodeIds()`.

`ComputedAttributeData<Real>` and `SampledNodeAttributeData<Real>` also
store the `NodeIdSpace` of the returned buffer. The default public computation
type is `ComputedAttributeData<float>` or
`SampledNodeAttributeData<float>`; selecting `Real=double` returns the
corresponding `double` specialization. Public computation methods can request a
different public node ID space, but projection always happens after the
internal pipeline has computed the result in `NodeIdSpace::MorphologicalTree`.

`NodeIdSpace::Higra` requests the preserved imported node domain.
`projectNodeValuesToExportedHigra(...)` aligns results with a fresh export.
Both layouts contain pixel leaves and internal nodes; pixel rows receive
unit-component values. See [Higra interoperability](higra-interoperability.md)
for domain validity and projection rules.

`computeAttributeMapping(...)` projects to the image domain: each pixel
receives the value stored at its smallest node.

## Numeric stability contract

Scalar attributes returned by the ordinary public attribute APIs are finite for
valid live nodes and exported pixel rows. This includes degenerate
supports such as one-pixel components, line-like components, zero continuous
bitquad perimeter, and bitquad configurations whose Euler estimate is zero.

The finite fallbacks are part of the public attribute contract:

| Attribute family | Degenerate condition | Returned value |
| --- | --- | --- |
| `BITQUAD_CIRCULARITY` | `BITQUAD_PERIMETER_CONTINUOUS <= eps` | `0` |
| `BITQUAD_PERIMETER_AVERAGE` | estimated Euler component count `<= 0` | `0` |
| `BITQUAD_LENGTH_AVERAGE` | estimated Euler component count `<= 0` | `0` |
| `BITQUAD_WIDTH_AVERAGE` | continuous perimeter `<= eps` | `0` |
| `ECCENTRICITY` | both inertia eigenvalues are numerically zero | `1` |
| `ECCENTRICITY` | smallest inertia eigenvalue `<= eps` or ratio overflows the practical range | `1e6` |

`BITQUAD_WIDTH_AVERAGE` is computed with the algebraically equivalent finite
formula `2 * BITQUAD_AREA / BITQUAD_PERIMETER_CONTINUOUS` when the continuous
perimeter is positive. This avoids the `inf / inf` pattern that appears if the
formula is evaluated through average area and average perimeter separately.

The finite-scalar contract does not mean every buffer cell is finite in every
API mode:

- dead internal node slots are addressable and may retain sentinel `NaN` values;
- `computeSampledNodeAttribute(...,
  MissingNodeAttributeSamplePolicy::NotANumber)` deliberately keeps `NaN` for
  missing ancestor or descendant samples;
- callers that provide non-finite floating-point altitudes are rejected by
  valued-tree and altitude-validation code before attributes are computed.

## Related guides

- [Attribute catalog](attribute-catalog.md): public scalar attributes, groups,
  and input contracts.
- [Distance-transform attributes](distance-transform.md): exact/approximate
  scalar contracts and units.
- [Attribute computer architecture](attribute-computer-architecture.md):
  contributor guide for adding or changing attribute computers.
- [Finite-window local-attribute C++ extension](finite-window-local-attributes.md):
  C++ extension contract for observation windows, anchored entries, visibility
  states, local decisions, and event algebras.
