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

Individual attributes may add stricter checks. For example, `MAX_DIST` requires
an altitude buffer, a regular 2D domain, uniform adjacency, and a globally
monotone altitude order. Its descriptive tree kind is irrelevant.

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
    std::vector<AttributeOrGroup>{AREA, BOX_WIDTH, BALANCE_NODE});
```

Returning values in a preserved public node ID space, for a tree created by
`createFromHigraParent(...)` and not edited since import:

```cpp
auto [names, values] = AttributeComputation::computeAttributes(
    importedTree,
    std::vector<AttributeOrGroup>{AREA, GrayLevelHeight},
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
smallest row-major pixel in the candidate support. Missing positions use one of
the typed `MissingNodeAttributeSamplePolicy` values: `RepeatNearest` (default),
`NotANumber`, or `Zero`.

When a radius contains multiple altitude distances, the implementation computes
support cardinalities and smallest-pixel tie-break keys once, then reuses those
metadata and the dense neighbourhood buffers for every distance. This changes
only storage lifetime: each distance still performs the same ancestor search
and `LargestSupportDescendant` selection.

Projecting node attributes to pixels or to an exported Higra layout:

```cpp
auto mapped = AttributeComputation::computeAttributeMapping(valuedTree, AREA);

auto internal = AttributeComputation::computeAttributes(
    valuedTree,
    std::vector<AttributeOrGroup>{AREA, GrayLevelHeight});
auto exported = AttributeComputation::projectNodeValuesToExportedHigra(
    valuedTree,
    internal.attributeNames(),
    internal.values());
```

## Python surface

Python keeps a smaller public surface than C++:

- image factories currently expose the canonical `np.uint8`
  `ValuedMorphologicalTree` path;
- `Attribute.compute_single_attribute(...)` and
  `Attribute.compute_attributes(...)` are the valued-tree attribute entry points;
- `Attribute.compute_single_topology_attribute(...)` and
  `Attribute.compute_topology_attributes(...)` are the explicit
  topology/support entry points;
- `NodeIdSpace` can be passed to the Python attribute methods when a preserved
  output node ID space is needed;
- attribute methods accept `dtype=np.float32` or `dtype=np.float64`; the default
  remains `np.float32`;
- every method taking an `Attribute` or `Attribute.Group` also takes its stable
  symbolic name as a string, matched exactly and case-sensitively;
- `AttributePipeline`, concrete C++ computers, and finite-window intermediate storage are not
  part of the Python API.

Typical Python calls use the same valued-tree versus topology/support split:

```python
gray_level_height_by_node = mmcfilters.Attribute.compute_single_attribute(
    valued_tree,
    mmcfilters.Attribute.GRAY_LEVEL_HEIGHT,
)
names, values = mmcfilters.Attribute.compute_attributes(
    valued_tree,
    [mmcfilters.Attribute.AREA, mmcfilters.Attribute.Group.GRAY_LEVEL],
)
topology_names, topology_values = mmcfilters.Attribute.compute_topology_attributes(
    valued_tree,
    [mmcfilters.Attribute.AREA, mmcfilters.Attribute.Group.BOUNDARY],
)
```

The same calls accept the symbolic names, which are the keys the returned layout
already uses. Enum values and names may be mixed in one request:

```python
gray_level_height_by_node = mmcfilters.Attribute.compute_single_attribute(
    valued_tree,
    "GRAY_LEVEL_HEIGHT",
)
names, values = mmcfilters.Attribute.compute_attributes(valued_tree, ["AREA", "GRAY_LEVEL"])
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
nodes should iterate `tree.aliveNodeIds()`.

`ComputedAttributeData<Real>` and `SampledNodeAttributeData<Real>` also
store the `NodeIdSpace` of the returned buffer. The default public computation
type is `ComputedAttributeData<float>` or
`SampledNodeAttributeData<float>`; selecting `Real=double` returns the
corresponding `double` specialization. Public computation methods can request a
different public node ID space, but projection always happens after the
internal pipeline has computed the result in `NodeIdSpace::MorphologicalTree`.

`NodeIdSpace::Higra` means the preserved imported Higra node ID domain. It is
available only for trees imported from Higra whose original node ID space has
not been invalidated by edits. Direct projection to this space copies live
internal-node rows and fills proper-part rows with unit-component values for
the requested attributes.

For a compact Higra layout exported from the current tree, use
`AttributeComputation::projectNodeValuesToExportedHigra(...)`. That
helper emits the `[proper parts | live internal nodes]` layout produced by
hierarchy export and fills unit proper-part rows through the responsible
attribute computers. `compute_attribute_mapping(...)` is the image-domain helper:
each proper part receives the value stored at its smallest node.

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
- [Attribute computer architecture](attribute-computer-architecture.md):
  contributor guide for adding or changing attribute computers.
- [Finite-window local attributes](finite-window-local-attributes.md):
  observation windows, anchored entries, visibility states, and additive local rules.
