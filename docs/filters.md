# Attribute Filters, Extinction Values, And UAO

This guide covers the public filtering layer in `mmcfilters/filters`. It
connects `AttributeFilters<T>`, `ExtinctionValues<T, Real>`,
`UltimateAttributeOpening<T, Real>`, and the MSER helper used by adaptive
filtering.

For the tree ownership, `NodeId`, proper-part, altitude, and mutation-version
model used here, see [Morphological Trees](trees.md). For computing the
node-indexed attributes consumed by these operators, see [Attributes](attributes.md).

## Public Model

The filtering layer reconstructs image-domain outputs from dense node-indexed
tree data. The common contracts are:

- tree nodes are addressed in internal dense `NodeId` space;
- attribute, criterion, score, and selection buffers have one element per
  internal node slot;
- image outputs use the original image-domain shape;
- object-style helpers snapshot the topology mutation version at construction
  time and reject use after topology mutation;
- plain attribute buffers are not versioned and must be recomputed after edits.

Use `WeightedMorphologicalTree<T>` when the operation needs owner state, such
as Python bindings, MSER-assisted filtering, `executeWithMSER`, image
reconstruction, or topology edits. Use `WeightedTreeView<T>` for read-only C++
filtering when the caller owns an external altitude buffer and does not need
MSER.

## Input Buffers

All node buffers are indexed by internal `NodeId`:

```text
buffer[node_id]
```

The expected sizes are:

- attributes: `tree.getNumInternalNodeSlots()` floating-point values;
- criteria: `tree.getNumInternalNodeSlots()` boolean values;
- scores: `tree.getNumInternalNodeSlots()` float values;
- UAO selection masks: `tree.getNumInternalNodeSlots()` byte/boolean values.

Dead internal slots may exist after edits. Buffers still keep the full slot
count; filtering traverses live nodes through the current tree topology.

Python arrays passed as attributes must be one-dimensional, C-contiguous
`np.float32` or `np.float64` arrays with length
`tree.numInternalNodeSlots`. Python boolean criteria are passed as lists or
vectors of the same length.

## AttributeFilters

`AttributeFilters<T>` groups the ordinary attribute-filter reconstruction
rules. It can be constructed from either a weighted owner or a non-owning
weighted view:

```cpp
#include <mmcfilters/attributes/Attributes.hpp>
#include <mmcfilters/filters/AttributeFilters.hpp>

using namespace mmcfilters;

auto boxHeightResult =
    AttributeComputation::computeSingleAttribute(weightedTree, BOX_HEIGHT);
const std::vector<float>& boxHeight = boxHeightResult.values();
auto boxHeight64 =
    AttributeComputation::computeSingleAttribute<double>(weightedTree, BOX_HEIGHT);

AttributeFilters<std::uint8_t> filters(weightedTree);

auto prunedMin = filters.filteringByPruningMin(boxHeight.data(), 2.0f);
auto prunedMax = filters.filteringByPruningMax(boxHeight.data(), 2.0f);
auto prunedMin64 = filters.filteringByPruningMin(boxHeight64.values().data(), 2.0);
```

When extracting one attribute column from a flat multi-attribute C++ result,
copy it to a dense per-node vector before calling filters. Multi-attribute
results are node-major, so one column is strided by `AttributeNames::NUM_ATTRIBUTES`.
A single-attribute result is already directly usable as a dense node buffer.

The object API allocates and returns output images. Static overloads write into
caller-owned output images and are useful in tight loops:

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

## Reconstruction Rules

Choose the reconstruction rule according to the intended attribute-filter
semantics:

- direct rule: accepted nodes use their own altitude and rejected nodes inherit
  the filtered level propagated from the parent;
- subtractive rule: accepted nodes add their altitude residue to the propagated
  parent level and rejected nodes inherit the parent level;
- subtractive score rule: each node residue is weighted by a dense float score
  before accumulation and the output is a float image;
- pruning-min rule: accepted branches remain traversable and rejected subtrees
  are reconstructed according to the pruning-min convention;
- pruning-max rule: rejected subtrees are detected bottom-up and reconstructed
  according to the pruning-max convention.

The criterion-based overloads receive the keep/reject decision directly. The
attribute-threshold overloads build the decision from a floating-point node
attribute and a threshold.

## Adaptive Criterion, MSER, And Depth Stability

`getAdaptiveCriterion(...)` adjusts a criterion with an MSER-style stability
analysis. It requires `WeightedMorphologicalTree<T>` ownership because MSER uses
the tree-owned altitude buffer to build delta neighbourhoods. Classical MSER
requires the topology to declare a globally monotone altitude order. The
standard max-tree and min-tree producers provide that capability:

```cpp
AttributeFilters<std::uint8_t> filters(weightedTree);

std::vector<bool> keep = /* one value per internal node slot */;
std::vector<bool> pruned = filters.getAdaptiveCriterion(
    keep,
    AltitudeDiff<std::uint8_t>{2});
```

`MSERComputer<T, Real>` is the advanced helper behind this path. For a positive
altitude delta, it pairs each node with an altitude-delta ascendant and
descendant, computes a variation score from an increasing attribute, and marks
strict local minima that pass the configured variation and attribute bounds. If
no attribute buffer is provided, it lazily computes `AREA`. Missing windows
produce `NaN`; nodes with `NaN` variation are not selected as MSERs.
Result getters such as `getVariation`, `getVariations`, and `getNumNodes` throw
`std::logic_error` until the first successful `computeMSER` call.

The variation score is:

```text
variation(x) = (attr(asc_delta(x)) - attr(desc_delta(x))) / attr(x)
```

Lower finite values are more stable. Adaptive pruning therefore moves a rejected
node to the minimum-variation representative among the node, its selected
ascendant, and its selected descendant.

MSER is not a general `WeightedTreeView<T>` operation. View-based filter objects
can run ordinary direct, subtractive, and pruning rules, but reject MSER-assisted
methods that need owner state.

The standard tree-of-shapes and self-dual residual-tree producers declare
`AltitudeOrder::UNCONSTRAINED`, so they do not provide the capability needed by
altitude MSER. Acceptance is based on this capability rather than the
descriptive tree kind. For unconstrained hierarchies, use depth stability:

```cpp
std::vector<bool> pruned = filters.getAdaptiveCriterionByDepth(
    keep,
    2);
```

`DepthStableRegionComputer<Real>` implements that rule. Here `depthDelta = 2`
means exactly two tree edges: climb two parent links for the ascendant and select
a descendant exactly two child links below the center. If several descendants
exist at that depth, the largest-area one is chosen; ties use the smallest
`NodeId`. This operator does not read altitude and should be described as
topological depth stability, not classical MSER. Its numeric score is still a
variation value; use `getVariation(...)`/`getVariations()` when inspecting it
directly. These result getters, the selected-node count, and the window-neighbour
getters throw `std::logic_error` until `computeByDepth(...)` succeeds.

## Extinction Values

For the scientific distinction among Cousty persistence, the former monotone
LCA projection, raster `contourMap`, and Xu shape-space saliency, see
[Hierarchy Saliency Maps](saliency.md). Its
[primary-reference mapping](saliency.md#primary-references-and-implementation-correspondence)
contains the complete citations and states which behavior is defined by each
paper versus added by this library. This section focuses on extinction selection,
filtering, and the filter-facing API.

`ExtinctionValues<T, Real>` ranks regional extrema by a scalar node attribute
on hierarchies that declare a globally monotone altitude order. Standard
max-tree and min-tree producers provide that capability. In this setting, the
regional extrema processed by the implementation are the tree leaves. The
attribute is a dense floating-point buffer indexed by internal `NodeId`; larger
attribute values are interpreted as stronger extrema. `Real` defaults to
`float`; select `double` when consuming double-precision attribute buffers.
Results are stored as `RegionalExtremaNode<Real>` records sorted by decreasing
extinction:

```cpp
#include <mmcfilters/filters/ExtinctionValues.hpp>

auto areaResult = AttributeComputation::computeSingleAttribute(weightedTree, AREA);
const std::vector<float>& area = areaResult.values();
auto area64Result = AttributeComputation::computeSingleAttribute<double>(weightedTree, AREA);

ExtinctionValues<std::uint8_t> extinction(weightedTree, area);
ExtinctionValues<std::uint8_t, double> extinction64(weightedTree, area64Result.values());
auto strongest = ExtinctionSelectionPolicy<float>::byTopK(8);
auto highExtinction = ExtinctionSelectionPolicy<float>::byThreshold(10.0f);

auto filtered = extinction.filtering(strongest);
auto filteredByThreshold = extinction.filtering(highExtinction);
auto contours = extinction.contourMap(strongest, ExtinctionContourScorePolicy::RankScore);

for (const RegionalExtremaNode<float>& item : extinction.getRegionalExtrema()) {
    NodeId leaf = item.leaf;
    NodeId cutoff = item.cutoffNode;
    float value = item.extinction;
}
```

`ExtinctionSelectionPolicy::byTopK(extremaToKeep)` selects the strongest extrema
by decreasing extinction ranking.
`ExtinctionSelectionPolicy::byThreshold(threshold)` selects every
extremum whose extinction value is greater than or equal to `threshold`.
`filtering(selection)` reconstructs an image from the selected extrema.
`contourMap(selection, scorePolicy)` writes an image-domain contour visualization
on the retained cutoff nodes. Use `ExtinctionContourScorePolicy::RankScore` for
dense rank scores, or `ExtinctionContourScorePolicy::ExtinctionValue` for raw
extinction values.
Rank counts must be non-negative and extinction thresholds must be finite. The
dominant extremum has no stronger merge point; its extinction is represented by
the explicit finite sentinel `std::numeric_limits<Real>::max()`.

`contourMap` is a visualization. For the edge-indexed saliency map of a hierarchy
in the quasi-flat-zone sense, use the formal extinction path:

```cpp
const std::vector<float>& valuation = extinction.getExtinctionValueAttribute();
auto edgeMap = extinction.computeFormalSaliencyEdgeMap();

std::vector<int> rankedValuation = extinction.computeRankedExtinctionValueAttribute();
auto rankedEdgeMap = extinction.computeRankedFormalSaliencyEdgeMap();
```

`getExtinctionValueAttribute()` returns a cached dense node
attribute: each regional-extremum leaf receives its extinction value, and every
non-leaf node receives the maximum extinction among the extrema contained in its
subtree. This max-propagated valuation is an intermediate quantity in the
hierarchical-watershed construction.

`computeFormalSaliencyEdgeMap()` follows Section 8.1 of Cousty et al. It builds
an altitude-ordered MST/BPTAO, assigns each binary merge the persistence

```text
persistence(node) = min(maxExtinction(left), maxExtinction(right))
```

and converts the persistence-weighted MST to its full-graph QFZ saliency map.
The dominant-extremum sentinel therefore does not propagate to ordinary merger
edges. The ranked variant ranks values already present on the final edge map.

The former `max descendants -> LCA` experiment remains available explicitly as
`computeMonotoneExtinctionProjection()`. It is a valid monotone hierarchy
projection, but it is not the Cousty hierarchical-watershed persistence method.

Standard trees of shapes and self-dual residual trees declare unconstrained
altitude order and are therefore rejected by `ExtinctionValues`. Their leaves
can be regional extrema, but the complete regional-extrema set is not generally
equivalent to `tree.getLeaves()`. The gate is the altitude-order capability,
not the descriptive kind. Support for these hierarchies is provided by the
separate `ShapeSpaceSaliency` API when the intended method is Xu shaping: it
treats all original tree nodes as a parent-child graph,
computes extrema and finite extinctions in that second shape space, and combines
the representative scores by maximum on the original region contours. See
[Morphological Trees](trees.md#shape-space-extinction-saliency).

## Ultimate Attribute Opening

`UltimateAttributeOpening<T, Real>` consumes a dense increasing-attribute buffer
and computes two image-domain outputs. `Real` defaults to `float`; use
`double` for double-precision attribute buffers:

- maximum contrast image: the largest selected altitude contrast for each
  pixel;
- associated image: the attribute index associated with that maximum contrast.

The basic workflow is:

```cpp
#include <mmcfilters/filters/UltimateAttributeOpening.hpp>

auto boxHeightResult = AttributeComputation::computeSingleAttribute(
    weightedTree,
    BOX_HEIGHT);
const std::vector<float>& boxHeight = boxHeightResult.values();
auto boxHeight64 = AttributeComputation::computeSingleAttribute<double>(
    weightedTree,
    BOX_HEIGHT);

UltimateAttributeOpening<std::uint8_t> uao(weightedTree, boxHeight);
UltimateAttributeOpening<std::uint8_t, double> uao64(weightedTree, boxHeight64.values());
uao.execute(maxCriterion);

auto contrast = uao.getMaxContrastImage();
auto associated = uao.getAssociatedImage();
auto colors = uao.getAssociatedColorImage();
```

`execute(maxCriterion)` treats all internal nodes as selectable candidates.
`execute(maxCriterion, selectedForFiltering)` accepts an explicit dense
selection mask. `executeWithMSER(maxCriterion, deltaMSER)` builds that mask from
classical altitude MSER and therefore requires a `WeightedMorphologicalTree<T>`
owner, not just a view. `executeWithDepthStability(maxCriterion, depthDelta)`
builds the mask from topological depth stability and does not use altitude for
the stability selection.

The increasing attribute is expected to encode the primitive scale used by UAO,
for example area, height, or another monotone criterion. The implementation
stores contrasts in the same type `T` as the altitude. Integral altitude types
therefore keep the historical API shape; use a wider or floating altitude type
in C++ when the contrast range matters.

## Python Usage

Python exposes the owner-oriented `uint8` workflow:

```python
import numpy as np
import mmcfilters

image = np.ascontiguousarray(image, dtype=np.uint8)
tree = mmcfilters.MorphologicalTreeFactory.createMaxTree(image, radius=1.5)

area = mmcfilters.Attribute.computeSingleAttribute(
    tree,
    mmcfilters.Attribute.AREA,
)
box_height = mmcfilters.Attribute.computeSingleAttribute(
    tree,
    mmcfilters.Attribute.BOX_HEIGHT,
)

filters = mmcfilters.AttributeFilters(tree)

keep_large = (area >= 4.0).tolist()
direct = filters.filteringDirectRule(keep_large)
subtractive = filters.filteringSubtractiveRule(keep_large)
pruned_min = filters.filteringByPruningMin(box_height, 2.0)
pruned_max = filters.filteringByPruningMax(box_height, 2.0)
score = filters.filteringSubtractiveScoreRule(area.tolist())
adaptive = filters.getAdaptiveCriterion(keep_large, delta=2)
depth_adaptive = filters.getAdaptiveCriterionByDepth(keep_large, depthDelta=2)

depth_stability = mmcfilters.DepthStableRegionComputer(tree)
depth_mask = depth_stability.computeByDepth(depthDelta=2)
depth_variation = depth_stability.getVariations()
```

Extinction values are available for hierarchies with a globally monotone
altitude-order capability, including the standard max-tree and min-tree
producers, either through `ExtinctionValues` directly or through
`AttributeFilters`:

```python
extinction = mmcfilters.ExtinctionValues(max_tree, area)
strongest = mmcfilters.ExtinctionSelectionPolicy.byTopK(8)
high_extinction = mmcfilters.ExtinctionSelectionPolicy.byThreshold(10.0)

filtered = extinction.filtering(strongest)
filtered_by_threshold = extinction.filtering(high_extinction)
contours = extinction.contourMap(
    strongest,
    mmcfilters.ExtinctionContourScorePolicy.RankScore,
)
records = extinction.getRegionalExtrema()
valuation = extinction.computeRankedExtinctionValueAttribute()
edge_map = extinction.computeFormalSaliencyEdgeMap(ranked=True)

filtered2 = filters.filteringByExtinction(area, strongest)
contours2 = filters.contourMapByExtinction(
    area,
    strongest,
    mmcfilters.ExtinctionContourScorePolicy.RankScore,
)
```

UAO follows the same dense attribute-buffer convention:

```python
uao = mmcfilters.UltimateAttributeOpening(tree, box_height)
uao.execute(maxCriterion=image.shape[0])

max_contrast = uao.getMaxContrastImage()
associated = uao.getAssociatedImage()
associated_color = uao.getAssociatedColoredImage()

uao.executeWithMSER(maxCriterion=image.shape[0], deltaMSER=2)
uao.executeWithDepthStability(maxCriterion=image.shape[0], depthDelta=2)
```

## Edits And Lifetime

Create filter helper objects after the topology edits that should affect the
operation:

```python
tree.mergeNodeIntoParent(node_id)

area = mmcfilters.Attribute.computeSingleAttribute(tree, mmcfilters.Attribute.AREA)
filters = mmcfilters.AttributeFilters(tree)
result = filters.filteringByPruningMin(area, 4.0)
```

Do not reuse `AttributeFilters`, `ExtinctionValues`, or
`UltimateAttributeOpening` objects after mutating the tree topology. The objects
record the mutation version and fail explicitly once the underlying topology
changes. Attribute arrays should also be recomputed after edits because their
values and live-node interpretation may have changed.

## Related Guides

- [Morphological Trees](trees.md): owner/view boundary, `NodeId`, proper parts,
  altitude, and topology mutation.
- [Attributes](attributes.md): computing the dense node-indexed buffers consumed
  by filters.
- [Python API Guide](python-api.md): Python construction, attributes, filters,
  contours, and Higra interop.
- [Incremental Contours](contours.md): contour extraction used by extinction
  saliency maps.
