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
auto output = ImageUInt8::create(tree.getNumRowsOfImage(), tree.getNumColsOfImage());

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

## Adaptive Criterion And MSER

`getAdaptiveCriterion(...)` adjusts a criterion with an MSER-style stability
analysis. It requires `WeightedMorphologicalTree<T>` ownership because MSER uses
the tree-owned altitude buffer to build delta neighbourhoods:

```cpp
AttributeFilters<std::uint8_t> filters(weightedTree);

std::vector<bool> keep = /* one value per internal node slot */;
std::vector<bool> pruned = filters.getAdaptiveCriterion(
    keep,
    AltitudeDiff<std::uint8_t>{2});
```

`MSERComputer<T, Real>` is the advanced-public helper behind this path. It
pairs each node with altitude-delta ascendant and descendant nodes, computes a
stability score from an increasing attribute, and marks strict local stability
minima that pass the configured variation and attribute bounds. If no attribute
buffer is provided, it lazily computes `AREA`.

MSER is not a general `WeightedTreeView<T>` operation. View-based filter objects
can run ordinary direct, subtractive, and pruning rules, but reject MSER-assisted
methods that need owner state.

## Extinction Values

`ExtinctionValues<T, Real>` ranks regional extrema by a scalar node attribute.
The attribute is a dense floating-point buffer indexed by internal `NodeId`.
`Real` defaults to `float`; select `double` when consuming double-precision
attribute buffers. Results are stored as `RegionalExtremaNode<Real>` records
sorted by decreasing extinction:

```cpp
#include <mmcfilters/filters/ExtinctionValues.hpp>

auto areaResult = AttributeComputation::computeSingleAttribute(weightedTree, AREA);
const std::vector<float>& area = areaResult.values();
auto area64Result = AttributeComputation::computeSingleAttribute<double>(weightedTree, AREA);

ExtinctionValues<std::uint8_t> extinction(weightedTree, area);
ExtinctionValues<std::uint8_t, double> extinction64(weightedTree, area64Result.values());
auto filtered = extinction.filtering(8);
auto saliency = extinction.saliencyMap(8, true);

for (const RegionalExtremaNode<float>& item : extinction.getExtinctionValues()) {
    NodeId leaf = item.leaf;
    NodeId cutoff = item.cutoffNode;
    float value = item.extinction;
}
```

`filtering(extremaToKeep)` reconstructs an image by retaining the strongest
extrema. `saliencyMap(extremaToKeep, unweighted)` writes saliency on compact
contours of the retained cutoff nodes. With `unweighted=true`, contours receive
rank-like scores; with `unweighted=false`, they receive extinction values.

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
MSER and therefore requires a `WeightedMorphologicalTree<T>` owner, not just a
view.

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
```

Extinction values are available either through `ExtinctionValues` directly or
through convenience methods on `AttributeFilters`:

```python
extinction = mmcfilters.ExtinctionValues(tree, area)
filtered = extinction.filtering(leafToKeep=8)
saliency = extinction.saliencyMap(leafToKeep=8, unweighted=True)
records = extinction.getExtinctionValues()

filtered2 = filters.filteringByExtinction(area, leafToKeep=8)
saliency2 = filters.saliencyMapByExtinction(
    area,
    leafToKeep=8,
    unweighted=True,
)
```

Pass `unweighted` explicitly when the saliency score convention matters. The
direct `ExtinctionValues.saliencyMap(...)` API and older convenience wrappers
have historically exposed different defaults.

UAO follows the same dense attribute-buffer convention:

```python
uao = mmcfilters.UltimateAttributeOpening(tree, box_height)
uao.execute(maxCriterion=image.shape[0])

max_contrast = uao.getMaxContrastImage()
associated = uao.getAssociatedImage()
associated_color = uao.getAssociatedColoredImage()

uao.executeWithMSER(maxCriterion=image.shape[0], deltaMSER=2)
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
