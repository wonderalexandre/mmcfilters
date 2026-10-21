# The Extremum Residual Tree: guide to the DGMM manuscript

This guide accompanies **“The Extremum Residual Tree: A Self-Dual Morphological
Tree Based on Levelings”**, prepared for submission to DGMM.

See [Illustrating Attribute Filtering with the Extremum Residual Tree and the Tree of Shapes](../../notebooks/Attribute_Filtering_Illustration_ExtResTree_ToS.ipynb)
for the Figure 1 example and attribute filtering on the cell image.

## Constructing an ext-res-tree

### Python

Install the library from the repository root:

```bash
python -m pip install .
```

Python accepts two-dimensional, C-contiguous `np.uint8` images. This example uses
the nine-point signal from Figure 1 as a one-row image:

```python
import mmcfilters
import numpy as np

image = np.array([[3, 3, 1, 2, 5, 4, 4, 2, 2]], dtype=np.uint8)
tree = mmcfilters.MorphologicalTreeFactory.create_extremum_residual_tree(image, radius=1.0)

assert np.array_equal(tree.reconstruct_from_node_altitudes(), image)
```

A radius of `1.0` gives 4-adjacency on a two-dimensional image and chain adjacency
on this one-row example. A radius of `1.5` gives 8-adjacency. An explicit symmetric
`RegularGridAdjacency2D` can also be supplied. See the [Python API guide](../python-api.md)
for the overloads.

### C++

The corresponding C++ example uses camelCase identifiers:

```cpp
#include <mmcfilters/trees/MorphologicalTreeFactory.hpp>

#include <cassert>
#include <cstdint>

int main() {
    using namespace mmcfilters;

    auto image = Image<std::uint8_t>::fromRaw(
        new std::uint8_t[]{3, 3, 1, 2, 5, 4, 4, 2, 2}, 1, 9);

    RegularGridAdjacency2D adjacency(1, 9, 1.0);
    auto tree = MorphologicalTreeFactory::createExtremumResidualTreeByCTs(
        image, adjacency);
    auto reconstruction = tree.reconstructFromNodeAltitudes();
    assert(reconstruction->isEqual(image));
}
```

The image owns the buffer passed to `fromRaw` and releases it with `delete[]`.
The result is a `ValuedMorphologicalTree<T>`. C++ supports the altitude types
accepted by the factory; the Python image contract above is more restrictive.
Both interfaces use increasing row-major pixel identifiers as the fixed spatial
order for equal-cardinality extrema.

## Algorithm 1 in the implementation

Algorithm 1 maintains a max-tree and a min-tree of the current image. A current
extremal flat zone is a leaf proper part in one of these trees, as in Eq. (29).
`ExtremumResidualTreeByCTsBuilder::build(image)` constructs the initial pair as
temporary state, then repeats the selection, recording, leveling, and
pending support updates until the image is constant. It returns an independent
`NativeHierarchyOwned<T>` with owned adjacency metadata. The factory
materializes this result as a `ValuedMorphologicalTree<T>`.

| Algorithm 1 | Implementation responsibility |
| --- | --- |
| Line 1: build the initial component trees | `ExtremumResidualTreeByCTsBuilder::build(image)`, using `ComponentTreeBuilder` and the shared native-hierarchy materialization |
| Lines 2–5: initialize and select the smallest current extremum | `ExtremumAgenda::initialize`, `ExtremumAgenda::next`, and `SelectedSupport::collectFromLeaf` |
| Lines 6–11: record the support and attach pending children | `ResidualTreeConstruction::addSelectedSupport` |
| Line 12: apply the elementary leveling | `ExtremumResidualTreeByCTsBuilder::applyElementaryLeveling`, using [DualMinMaxTreeIncrementalFilterLeaf.hpp](../../mmcfilters/trees/detail/DualMinMaxTreeIncrementalFilterLeaf.hpp) |
| Line 13: update pending supports in both component trees | `ResidualTreeConstruction::updatePendingSupports`, using transfers reported by the component-tree updater |
| Line 14: advance the scale index | The number of recorded nonterminal supports determines the next step |
| Lines 16–18: complete the terminal root and attach remaining supports | `ResidualTreeConstruction::completeTree` |

The algorithm is implemented by one public class,
[ExtremumResidualTreeByCTsBuilder](../../mmcfilters/trees/ExtremumResidualTreeByCTsBuilder.hpp).
