# Self-dual residual trees: implementation and demonstration guide

This page is the software companion to the paper **“Self-Dual Residual Trees.”** It does not repeat the
definitions or proofs. Its purpose is to locate the implementations, explain
how the mathematical objects are represented by `mmcfilters`, and point to the
executable demonstration.

## Public construction API

The public factories first build a max-tree and a min-tree with the same
symmetric adjacency, then pass both trees to the corresponding synchronized
builder.

### Python

```python
import mmcfilters
import numpy as np

image = np.ascontiguousarray(image, dtype=np.uint8)

unrestricted = (
    mmcfilters.MorphologicalTreeFactory.createSelfDualResidualTree(
        image,
        radius=1.0,  # shared 4-adjacency
    )
)

saturated = (
    mmcfilters.MorphologicalTreeFactory.createSaturatedSelfDualResidualTree(
        image,
        infinityPixel=0,
        radius=1.0,
    )
)

assert np.array_equal(unrestricted.reconstructionImage(), image)
assert np.array_equal(saturated.reconstructionImage(), image)
```

Python accepts two-dimensional, C-contiguous `np.uint8` images. A radius of
`1.0` selects 4-adjacency and `1.5` selects 8-adjacency. An explicit symmetric
`RegularGridAdjacency2D` can be supplied instead. Pixel identifiers, including
`infinityPixel`, use row-major order.

### C++

```cpp
#include <mmcfilters/trees/MorphologicalTreeFactory.hpp>

using namespace mmcfilters;

RegularGridAdjacency2D adjacency(rows, cols, 1.0);

auto unrestricted =
    MorphologicalTreeFactory::createSelfDualResidualTree(image, adjacency);

auto saturated =
    MorphologicalTreeFactory::createSaturatedSelfDualResidualTree(
        image, adjacency, NodeId{0});
```

The returned object is a `WeightedMorphologicalTree<T>`. For a non-root node
`n`, the implementation exposes the paper notation as follows:

| Paper object | Public query |
| --- | --- |
| residual support `X_k` | `tree.reconstructNode(n)` or `tree.getConnectedComponent(n)` |
| valuation `eta(X_k)` | `tree.getAltitude(n)` |
| parent `parent(X_k)` | `tree.getNodeParent(n)` |
| signed residual `r_k = eta(X_k) - eta(parent(X_k))` | `tree.getNodeResidue(n)` |
| proper part `rho(X_k)` | `tree.getProperParts(n)` |
| exact image represented by the valued tree | `tree.reconstructionImage()` |

Residual-tree altitudes alternate with polarity, so the returned hierarchy has
`AltitudeOrder::UNCONSTRAINED`. Algorithms that require a globally monotone
altitude must instead use a structural or otherwise increasing node valuation,
such as support area.

The default `SdrtTiePolicy::ContrastInvariantSpatial` implements the
polarity-independent area/spatial ordering used for the self-dual constructions
in the paper. `SdrtTiePolicy::MaxBeforeMinThenSpatial` is also available for
experiments that intentionally impose a polarity order.

## From the paper to the source code

The production path is split into small components so that selection,
synchronized evolution, saturation certification, and tree assembly can be
tested separately.

| Paper operation | Implementation |
| --- | --- |
| Construct the initial synchronized `T_max(I^0)` and `T_min(I^0)` | [`MorphologicalTreeFactory.hpp`](https://github.com/wonderalexandre/MorphologicalAttributeFilters/blob/main/mmcfilters/trees/MorphologicalTreeFactory.hpp) |
| Select the construction mode | [`UnrestrictedResidualTreeBuilder.hpp`](https://github.com/wonderalexandre/MorphologicalAttributeFilters/blob/main/mmcfilters/trees/sdrt/UnrestrictedResidualTreeBuilder.hpp) and [`SaturatedResidualTreeBuilder.hpp`](https://github.com/wonderalexandre/MorphologicalAttributeFilters/blob/main/mmcfilters/trees/sdrt/SaturatedResidualTreeBuilder.hpp) |
| Maintain the current flat zones of `I^k` | [`detail/FlatZonePartition.hpp`](https://github.com/wonderalexandre/MorphologicalAttributeFilters/blob/main/mmcfilters/trees/sdrt/detail/FlatZonePartition.hpp) |
| Maintain `M(I^k)` and apply the key `K_prec` | [`detail/ResidualTreeCandidateAgenda.hpp`](https://github.com/wonderalexandre/MorphologicalAttributeFilters/blob/main/mmcfilters/trees/sdrt/detail/ResidualTreeCandidateAgenda.hpp) and [`SdrtTiePolicy.hpp`](https://github.com/wonderalexandre/MorphologicalAttributeFilters/blob/main/mmcfilters/trees/sdrt/SdrtTiePolicy.hpp) |
| Gather `X_k`, its boundary, and its primal/dual owners | [`detail/ResidualTreeCandidatePreparation.hpp`](https://github.com/wonderalexandre/MorphologicalAttributeFilters/blob/main/mmcfilters/trees/sdrt/detail/ResidualTreeCandidatePreparation.hpp) |
| Evaluate saturated eligibility `chi_sat(X_k)` | [`detail/SaturatedResidualEligibility.hpp`](https://github.com/wonderalexandre/MorphologicalAttributeFilters/blob/main/mmcfilters/trees/sdrt/detail/SaturatedResidualEligibility.hpp), assisted by [`detail/SaturatedDynamicLca.hpp`](https://github.com/wonderalexandre/MorphologicalAttributeFilters/blob/main/mmcfilters/trees/sdrt/detail/SaturatedDynamicLca.hpp) |
| Move the selected extremum to its first merging level and update both component trees | [`detail/MinMaxResidualTreeEngine.hpp`](https://github.com/wonderalexandre/MorphologicalAttributeFilters/blob/main/mmcfilters/trees/sdrt/detail/MinMaxResidualTreeEngine.hpp) and [`adjust/DualMinMaxTreeIncrementalFilterLeaf.hpp`](https://github.com/wonderalexandre/MorphologicalAttributeFilters/blob/main/mmcfilters/trees/adjust/DualMinMaxTreeIncrementalFilterLeaf.hpp) |
| Record the selected supports and build the inclusion hierarchy | [`detail/ResidualTreeEventAssembler.hpp`](https://github.com/wonderalexandre/MorphologicalAttributeFilters/blob/main/mmcfilters/trees/sdrt/detail/ResidualTreeEventAssembler.hpp) |
| Validate proper-part ownership and exact reconstruction, then materialize the valued tree | [`detail/ResidualTreeMaterialization.hpp`](https://github.com/wonderalexandre/MorphologicalAttributeFilters/blob/main/mmcfilters/trees/sdrt/detail/ResidualTreeMaterialization.hpp) |

The common control loop is `MinMaxResidualTreeEngine<T, mode>`. At each event
it obtains the selected leaf from the ordered agenda, reads the first merging
level from the leaf's parent in the primal component tree, updates the primal
and dual trees together, merges the corresponding current flat zones, and
emits one residual-tree node. The saturated specialization adds complement
connectivity certification relative to `infinityPixel`; the unrestricted
specialization carries no exterior or saturation state.

## Demonstration notebook

[`Self_Dual_Residual_Trees_Tutorial.ipynb`](https://github.com/wonderalexandre/MorphologicalAttributeFilters/blob/main/notebooks/Self_Dual_Residual_Trees_Tutorial.ipynb)
