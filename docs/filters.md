# Filters and hierarchy operators

This guide describes the filtering operators in `mmcfilters/filters`, their
node-buffer contracts, and their image-domain outputs. See [Attributes](attributes.md)
for computing the buffers consumed here.

## Choose an operator

MSER denotes maximally stable extremal regions.

| Need | Operator |
| --- | --- |
| reconstruct from explicit keep/reject decisions | `NodePreservationMask` with a direct or subtractive filter |
| remove complete branches by an attribute threshold | pruning-min or pruning-max rule |
| adjust preservation decisions with altitude stability | altitude-stability adjustment |
| adjust preservation decisions by a number of tree edges | depth-stability adjustment |
| rank and select component-tree extrema | `ExtinctionValues` |
| obtain the strongest response over an attribute scale | `UltimateAttributeOpening` |

## Shared contracts

The filtering layer reconstructs images from dense internal-node data:

- buffers use internal `NodeId` indexing;
- their length is `tree.numInternalNodeSlots()`;
- dead slots remain present, while algorithms traverse live nodes;
- image outputs use the tree's regular 2D domain;
- helper objects capture the topology mutation version and reject use after a
  topology change;
- plain attribute arrays are not versioned and must be recomputed after edits.

Expected C++ buffers are:

| Input | Element type | Shape |
| --- | --- | --- |
| attribute | floating point | one value per internal node slot |
| node-preservation mask | `bool` | one decision per internal node slot; `true` preserves |
| score | `float` | one value per internal node slot |
| Ultimate Attribute Opening (UAO) selection | byte/boolean | one value per internal node slot |

Python attributes must be one-dimensional C-contiguous `np.float32` or
`np.float64` arrays of the same length. A column sliced from a multi-attribute
matrix is usually not contiguous; copy it or compute the scalar attribute
directly.

Use the owning `ValuedMorphologicalTree<T>` when an operator needs topology
and altitude state, including Python bindings, image reconstruction, MSER, or
`executeWithMSER`. Read-only C++ filtering rules can also use
`ValuedMorphologicalTreeView<T>` when the caller owns the altitude buffer.

## Attribute filters

`AttributeFilters<T>` groups pruning reconstruction from either an explicit
`NodePreservationMask` or a node-attribute threshold:

```cpp
#include <mmcfilters/attributes/Attributes.hpp>
#include <mmcfilters/filters/AttributeFilters.hpp>

using namespace mmcfilters;

auto boundingBoxHeight = AttributeComputation::computeSingleAttribute(
    valuedTree,
    BoundingBoxHeight);

AttributeFilters<std::uint8_t> filters(valuedTree);
auto prunedMin = filters.filteringByPruningMin(
    boundingBoxHeight.values().data(),
    2.0f);
auto prunedMax = filters.filteringByPruningMax(
    boundingBoxHeight.values().data(),
    2.0f);
```

A single-attribute result is a dense node buffer. In a multi-attribute result,
columns are strided because the layout is node-major; copy the selected column
before passing it to a filter.

### Reconstruction rules

- **Direct:** `DirectAttributeFilter` accepts a `NodePreservationMask`; preserved
  nodes use their altitude and rejected nodes inherit the reconstructed parent
  altitude. The root must be preserved.
- **Hard subtractive:** `SubtractiveAttributeFilter` independently gates
  every zero-baseline node residue, including the root residue, with a
  `NodePreservationMask`. Integral-altitude trees produce a signed image.
- **Soft subtractive:** `SoftSubtractiveAttributeFilter` multiplies every
  zero-baseline node residue by a finite dense score in `[0, 1]`; the output has
  the score dtype.
- **Pruning-min:** accepted branches remain traversable and rejected subtrees
  follow the pruning-min reconstruction convention.
- **Pruning-max:** rejected subtrees are detected bottom-up and follow the
  pruning-max reconstruction convention.

`NodePreservationMask` uses `true` for preservation. `NodePruningMask` uses
`true` for pruning; convert between them only with `toNodePruningMask(...)` or
`toNodePreservationMask(...)`.

Direct reconstruction is separate from subtractive residue modulation:

```cpp
NodePreservationMask keep(
    valuedTree.topology().numInternalNodeSlots(),
    true);

auto direct = DirectAttributeFilter<std::uint8_t>(valuedTree)
                  .applyDirectAttributeFilter(keep);
auto subtractive = SubtractiveAttributeFilter<std::uint8_t>(valuedTree)
                       .applySubtractiveAttributeFilter(keep);
```

## Stability

### Altitude-based MSER stability

`adjustNodePreservationMaskByAltitudeStability(...)` relocates the rejections
in a `NodePreservationMask` using an MSER-style altitude window. It requires a
`ValuedMorphologicalTree<T>` with a globally monotone altitude order. Standard
max-tree and min-tree factories provide this capability.

```cpp
std::vector<float> area = /* one value per internal node slot */;
NodePreservationMask keep = computeNodePreservationMask<float>(area, 4.0f);

NodePreservationMask adjusted = adjustNodePreservationMaskByAltitudeStability(
    valuedTree, keep, AltitudeDifference<std::uint8_t>{2});
```

`MSERComputer<T, Real>` pairs each node with altitude-window ancestor and
descendant samples, computes variation from an increasing attribute, and
selects strict local minima within configured bounds. If no attribute is
provided, it computes `AREA`. Missing windows produce `NaN` and are not
selected.

The variation is

```text
variation(x) = (attr(ancestor(x)) - attr(descendant(x))) / attr(x)
```

Result getters require a successful computation. A `ValuedMorphologicalTreeView<T>` can
run ordinary reconstruction rules but cannot run the MSER path that requires an
owning `ValuedMorphologicalTree<T>`.

### Depth stability

Use depth stability on hierarchies without a global altitude direction or when
the window should represent a number of parent/child edges:

```cpp
NodePreservationMask adjusted = adjustNodePreservationMaskByDepthStability(
    valuedTree, keep, 2);
```

`DepthStableRegionComputer<Real>` climbs exactly `depthWindowRadius` parent links and
selects a descendant at the same depth. When several descendants qualify, it
uses the largest `AREA`, then the smallest `NodeId`. This is a topological
stability operator, not altitude-based MSER stability.

When either window is incomplete, the default
`IncompleteStabilityWindowPolicy::PreserveInputDecision` retains the received
decision at that node. This does not necessarily preserve the node: an input
rejection remains a rejection. Attribute thresholding is intentionally a
separate operation from stability adjustment.

## Extinction values

`ExtinctionValues<T, Real>` ranks component-tree leaves from an increasing node
attribute. It requires a globally monotone altitude order. `Real` defaults to
`float`; use `double` with double-precision attribute buffers.

```cpp
#include <mmcfilters/filters/ExtinctionValues.hpp>

auto area = AttributeComputation::computeSingleAttribute(valuedTree, AREA);
ExtinctionValues<std::uint8_t> extinction(valuedTree, area.values());

auto strongest = ExtinctionSelectionPolicy<float>::byTopK(8);
auto aboveThreshold =
    ExtinctionSelectionPolicy<float>::byThreshold(10.0f);

auto filtered = extinction.filtering(strongest);
auto contourImage = extinction.contourMap(
    strongest,
    ExtinctionContourScorePolicy::RankScore);
```

`byTopK(k)` selects the first `k` extrema in decreasing extinction order.
`byThreshold(t)` selects extinctions greater than or equal to a finite
threshold. The dominant extremum has no stronger merge point and receives the
finite ordering sentinel `std::numeric_limits<Real>::max()`.

`filtering(selection)` reconstructs selected extrema. `contourMap(...)` is a
pixel visualization of selected cutoff-node contours. Edge-indexed hierarchy
maps are provided separately:

```cpp
auto edgeMap = extinction.computeFormalSaliencyEdgeMap();
auto rankedEdgeMap = extinction.computeRankedFormalSaliencyEdgeMap();
auto directProjection = extinction.computeMonotoneExtinctionProjection();
```

See [Saliency maps](saliency.md) for the distinction between persistence and
direct hierarchy projection. `ShapeSpaceSaliency` is the appropriate operator
when extinction must be computed from an arbitrary attribute over the tree-node
graph.

## Ultimate Attribute Opening

`UltimateAttributeOpening<T, Real>` consumes an increasing node attribute and
produces:

- a maximum-contrast image;
- an associated-attribute image;
- an optional colored visualization of the associated image.

```cpp
#include <mmcfilters/filters/UltimateAttributeOpening.hpp>

auto boundingBoxHeight = AttributeComputation::computeSingleAttribute(
    valuedTree,
    BoundingBoxHeight);

UltimateAttributeOpening<std::uint8_t> uao(
    valuedTree,
    boundingBoxHeight.values());
uao.execute(maximumAttributeThreshold);

auto contrast = uao.getMaxContrastImage();
auto associated = uao.getAssociatedImage();
```

`execute(maximumAttributeThreshold)` considers all live nodes.
`execute(maximumAttributeThreshold, selectedForFiltering)` accepts an explicit
dense UAO primitive-selection mask; it is not a preservation/pruning mask.
`executeWithMSER(...)` builds that mask with altitude-based MSER stability and
requires `ValuedMorphologicalTree<T>`.
`executeWithDepthStability(...)` uses topological depth stability.

Contrasts use the altitude type `T`. Choose a wider or floating-point altitude
type in C++ when the contrast range requires it.

## Python usage

Python uses the same dense-buffer contracts:

```python
area = mmcfilters.Attribute.compute_single_topology_attribute(
    tree,
    mmcfilters.Attribute.AREA,
)
box_height = mmcfilters.Attribute.compute_single_topology_attribute(
    tree,
    mmcfilters.Attribute.BOUNDING_BOX_HEIGHT,
)

filters = mmcfilters.AttributeFilters(tree)
decisions = area >= 4.0
decisions[tree.root] = True
keep = mmcfilters.NodePreservationMask(decisions)
direct = mmcfilters.DirectAttributeFilter(tree).apply(keep)
pruned = filters.filtering_by_pruning_min(box_height, 2.0)

extinction = mmcfilters.ExtinctionValues(tree, area)
strongest = mmcfilters.ExtinctionSelectionPolicy.by_top_k(8)
filtered = extinction.filtering(strongest)

uao = mmcfilters.UltimateAttributeOpening(tree, box_height)
uao.execute(maximum_attribute_threshold=image.shape[0])
```

`mmcfilters.compute_node_preservation_mask(area, 4.0)` builds the same mask from
the inclusive `>=` rule when no per-node adjustment is needed.

An operator may also receive the attribute itself, as an `Attribute` value or its
symbolic name, and compute the buffer internally. The declared capability
requirements select the valued or the topology entry point, so the result matches
the buffer the caller would have computed:

```python
pruned = filters.filtering_by_pruning_min("BOUNDING_BOX_HEIGHT", 2.0)
extinction = mmcfilters.ExtinctionValues(tree, "AREA")
uao = mmcfilters.UltimateAttributeOpening(tree, "BOUNDING_BOX_HEIGHT")
```

Pass the buffer instead when it is reused across several thresholds, which avoids
recomputing it per call.

The three reconstruction filters expose `apply(...)` as a short alias of their
explicit method name.

See [Python API](python-api.md) for array and dtype requirements.

## Edits and lifetime

Create helpers after the topology edits that should affect an operation. After a
topology mutation, recompute attributes and construct new `AttributeFilters`,
`ExtinctionValues`, and `UltimateAttributeOpening` objects. Stale objects reject
public reads explicitly.

## Related guides

- [Morphological trees](trees.md): smallest-node mapping, `NodeId`, altitude, and mutation.
- [Attributes](attributes.md): dense node-indexed input buffers.
- [Saliency maps](saliency.md): edge-indexed hierarchy operators.
- [Editing API](editing-api.md): mutation and derived-state lifetime.
