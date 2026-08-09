# Filters and hierarchy operators

This guide describes the filtering operators in `mmcfilters/filters`, their
node-buffer contracts, and their image-domain outputs. See [Attributes](attributes.md)
for computing the buffers consumed here.

## Choose an operator

MSER denotes maximally stable extremal regions.

| Need | Operator |
| --- | --- |
| reconstruct from a keep/reject criterion | `AttributeFilters` direct or subtractive rule |
| remove complete branches by an attribute threshold | pruning-min or pruning-max rule |
| adapt a criterion with altitude stability | MSER-assisted filtering |
| adapt a criterion by a number of tree edges | depth stability |
| rank and select component-tree extrema | `ExtinctionValues` |
| obtain the strongest response over an attribute scale | `UltimateAttributeOpening` |

## Shared contracts

The filtering layer reconstructs images from dense internal-node data:

- buffers use internal `NodeId` indexing;
- their length is `tree.getNumInternalNodeSlots()`;
- dead slots remain present, while algorithms traverse live nodes;
- image outputs use the tree's regular 2D domain;
- helper objects capture the topology mutation version and reject use after a
  topology change;
- plain attribute arrays are not versioned and must be recomputed after edits.

Expected C++ buffers are:

| Input | Element type | Shape |
| --- | --- | --- |
| attribute | floating point | one value per internal node slot |
| criterion | `bool` | one value per internal node slot |
| score | `float` | one value per internal node slot |
| Ultimate Attribute Opening (UAO) selection | byte/boolean | one value per internal node slot |

Python attributes must be one-dimensional C-contiguous `np.float32` or
`np.float64` arrays of the same length. A column sliced from a multi-attribute
matrix is usually not contiguous; copy it or compute the scalar attribute
directly.

Use the owning `WeightedMorphologicalTree<T>` when an operator needs topology
and altitude state, including Python bindings, image reconstruction, MSER, or
`executeWithMSER`. Read-only C++ filtering rules can also use
`WeightedTreeView<T>` when the caller owns the altitude buffer.

## Attribute filters

`AttributeFilters<T>` groups criterion- and threshold-based reconstruction
rules:

```cpp
#include <mmcfilters/attributes/Attributes.hpp>
#include <mmcfilters/filters/AttributeFilters.hpp>

using namespace mmcfilters;

auto boxHeight = AttributeComputation::computeSingleAttribute(
    weightedTree,
    BOX_HEIGHT);

AttributeFilters<std::uint8_t> filters(weightedTree);
auto prunedMin = filters.filteringByPruningMin(
    boxHeight.values().data(),
    2.0f);
auto prunedMax = filters.filteringByPruningMax(
    boxHeight.values().data(),
    2.0f);
```

A single-attribute result is a dense node buffer. In a multi-attribute result,
columns are strided because the layout is node-major; copy the selected column
before passing it to a filter.

### Reconstruction rules

- **Direct:** accepted nodes use their altitude; rejected nodes inherit the
  filtered parent level.
- **Subtractive:** accepted nodes add their altitude residue to the propagated
  parent level; rejected nodes inherit the parent level.
- **Subtractive score:** each residue is multiplied by a dense score; the output
  is a float image.
- **Pruning-min:** accepted branches remain traversable and rejected subtrees
  follow the pruning-min reconstruction convention.
- **Pruning-max:** rejected subtrees are detected bottom-up and follow the
  pruning-max reconstruction convention.

Criterion overloads receive a keep/reject mask. Threshold overloads derive that
mask from one node attribute and a scalar threshold.

Static overloads write into caller-owned output images and are useful when an
application wants to reuse storage:

```cpp
std::vector<bool> keep(tree.getNumInternalNodeSlots(), true);
auto output = ImageUInt8::create(
    tree.getNumRowsOfGridDomain2D(),
    tree.getNumColsOfGridDomain2D());

AttributeFilters<std::uint8_t>::filteringByDirectRule(
    weightedTree,
    keep,
    output);
```

## Stability

### Altitude-based MSER stability

`getAdaptiveCriterion(...)` adjusts a criterion using an MSER-style altitude
window. It requires a `WeightedMorphologicalTree<T>` with a globally monotone
altitude order. Standard max-tree and min-tree factories provide this
capability.

```cpp
AttributeFilters<std::uint8_t> filters(weightedTree);
std::vector<bool> keep = /* one value per internal node slot */;

std::vector<bool> adaptive = filters.getAdaptiveCriterion(
    keep,
    AltitudeDiff<std::uint8_t>{2});
```

`MSERComputer<T, Real>` pairs each node with altitude-delta ancestor and
descendant samples, computes variation from an increasing attribute, and
selects strict local minima within configured bounds. If no attribute is
provided, it computes `AREA`. Missing windows produce `NaN` and are not
selected.

The variation is

```text
variation(x) = (attr(asc_delta(x)) - attr(desc_delta(x))) / attr(x)
```

Result getters require a successful computation. A `WeightedTreeView<T>` can
run ordinary reconstruction rules but cannot run the MSER path that requires an
owning `WeightedMorphologicalTree<T>`.

### Depth stability

Use depth stability on hierarchies without a global altitude direction or when
the window should represent a number of parent/child edges:

```cpp
std::vector<bool> adaptive = filters.getAdaptiveCriterionByDepth(
    keep,
    2);
```

`DepthStableRegionComputer<Real>` climbs exactly `depthDelta` parent links and
selects a descendant at the same depth. When several descendants qualify, it
uses the largest `AREA`, then the smallest `NodeId`. This is a topological
stability operator, not altitude-based MSER stability.

## Extinction values

`ExtinctionValues<T, Real>` ranks component-tree leaves from an increasing node
attribute. It requires a globally monotone altitude order. `Real` defaults to
`float`; use `double` with double-precision attribute buffers.

```cpp
#include <mmcfilters/filters/ExtinctionValues.hpp>

auto area = AttributeComputation::computeSingleAttribute(weightedTree, AREA);
ExtinctionValues<std::uint8_t> extinction(weightedTree, area.values());

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

auto boxHeight = AttributeComputation::computeSingleAttribute(
    weightedTree,
    BOX_HEIGHT);

UltimateAttributeOpening<std::uint8_t> uao(
    weightedTree,
    boxHeight.values());
uao.execute(maxCriterion);

auto contrast = uao.getMaxContrastImage();
auto associated = uao.getAssociatedImage();
```

`execute(maxCriterion)` considers all live nodes.
`execute(maxCriterion, selectedForFiltering)` accepts an explicit dense mask.
`executeWithMSER(...)` builds that mask with altitude-based MSER stability and
requires `WeightedMorphologicalTree<T>`.
`executeWithDepthStability(...)` uses topological depth stability.

Contrasts use the altitude type `T`. Choose a wider or floating-point altitude
type in C++ when the contrast range requires it.

## Python usage

Python uses the same dense-buffer contracts:

```python
area = mmcfilters.Attribute.computeSingleTopologyAttribute(
    tree,
    mmcfilters.Attribute.AREA,
)
box_height = mmcfilters.Attribute.computeSingleTopologyAttribute(
    tree,
    mmcfilters.Attribute.BOX_HEIGHT,
)

filters = mmcfilters.AttributeFilters(tree)
keep = (area >= 4.0).tolist()
direct = filters.filteringDirectRule(keep)
pruned = filters.filteringByPruningMin(box_height, 2.0)

extinction = mmcfilters.ExtinctionValues(tree, area)
strongest = mmcfilters.ExtinctionSelectionPolicy.byTopK(8)
filtered = extinction.filtering(strongest)

uao = mmcfilters.UltimateAttributeOpening(tree, box_height)
uao.execute(maxCriterion=image.shape[0])
```

See [Python API](python-api.md) for array and dtype requirements.

## Edits and lifetime

Create helpers after the topology edits that should affect an operation. After a
topology mutation, recompute attributes and construct new `AttributeFilters`,
`ExtinctionValues`, and `UltimateAttributeOpening` objects. Stale objects reject
public reads explicitly.

## Related guides

- [Morphological trees](trees.md): ownership, `NodeId`, altitude, and mutation.
- [Attributes](attributes.md): dense node-indexed input buffers.
- [Saliency maps](saliency.md): edge-indexed hierarchy operators.
- [Editing API](editing-api.md): mutation and derived-state lifetime.
