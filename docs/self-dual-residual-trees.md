# Self-dual residual trees: implementation and demonstration guide

This page is the software companion to the paper **“Self-Dual Residual Trees
by Synchronized Evolution of Component Trees.”** It does not repeat the
definitions or proofs. Its purpose is to locate the implementations, explain
how the mathematical objects are represented by `mmcfilters`, and point to the
executable demonstration.

The paper can link directly to this page at:

<https://github.com/wonderalexandre/MorphologicalAttributeFilters/blob/main/docs/self-dual-residual-trees.md>

For a reproducible archival reference, replace `main` with the release tag or
commit used for the paper experiments.

## Start here

| Resource | Location | Purpose |
| --- | --- | --- |
| Executable companion | [`notebooks/Self_Dual_Residual_Trees_Tutorial.ipynb`](https://github.com/wonderalexandre/MorphologicalAttributeFilters/blob/main/notebooks/Self_Dual_Residual_Trees_Tutorial.ipynb) | Recreates Figure 1, checks the structural properties, and reproduces the Figure 2 hydrant experiment. |
| Figure 2 input | [`dat/hydrant.png`](https://github.com/wonderalexandre/MorphologicalAttributeFilters/blob/main/dat/hydrant.png) | Original 363x352-pixel grayscale image used by the reproducible notebook experiment. |
| Public C++ factories | [`mmcfilters/trees/MorphologicalTreeFactory.hpp`](https://github.com/wonderalexandre/MorphologicalAttributeFilters/blob/main/mmcfilters/trees/MorphologicalTreeFactory.hpp) | Recommended C++ entry points for unrestricted and saturated residual trees. |
| Python bindings | [`pybinds/TreeBindings.cpp`](https://github.com/wonderalexandre/MorphologicalAttributeFilters/blob/main/pybinds/TreeBindings.cpp) | Exposes the two factories, the tie policy, and the weighted-tree query API. |
| User-facing API notes | [`docs/trees.md`](trees.md) and [`docs/python-api.md`](python-api.md) | Documents input contracts, adjacency choices, tree semantics, attributes, and reconstruction. |
| Main regression test | [`unit-tests/trees/sdrt/test_min_max_residual_tree_factory.cpp`](https://github.com/wonderalexandre/MorphologicalAttributeFilters/blob/main/unit-tests/trees/sdrt/test_min_max_residual_tree_factory.cpp) | Tests reconstruction, self-duality, construction policies, altitude types, edge cases, and custom adjacencies. |
| Construction benchmark | [`benchmarks/min_max_residual_tree_benchmark.cpp`](https://github.com/wonderalexandre/MorphologicalAttributeFilters/blob/main/benchmarks/min_max_residual_tree_benchmark.cpp) | Measures unrestricted and saturated construction on an input image. |

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

The two residual-tree factories implement the **shared-adjacency** constructions
of the paper. The complementary-adjacency result that coincides with the tree
of shapes is a structural theorem under an admissible topographic convention;
it is not exposed as a third mode of these residual-tree factories. The library
provides `MorphologicalTreeFactory::createTreeOfShapes` separately, and the
notebook uses it only for comparison.

## Demonstration notebook

[`Self_Dual_Residual_Trees_Tutorial.ipynb`](https://github.com/wonderalexandre/MorphologicalAttributeFilters/blob/main/notebooks/Self_Dual_Residual_Trees_Tutorial.ipynb)
is an executable companion rather than a second presentation of the theory. It:

1. recreates the exact 17x17 six-flat-zone image from Figure 1;
2. builds the unrestricted and saturated shared-4-adjacency trees;
3. reads the supports and signed residues from the public API;
4. recovers `I^0, ..., I^q` using the recurrence in Equation (5);
5. displays the Hasse diagrams and the crossing saturated sets of Figure 1(c);
6. checks the telescoping decomposition, nondecreasing support sizes,
   laminarity, non-empty proper parts, saturated eligibility, and contrast
   self-duality on the example;
7. uses the Figure 1 image for compact area-pruning and adjacency diagnostics;
   and
8. reproduces the complete Figure 2 hydrant experiment with the original
   363x352-pixel image, the three area thresholds `639`, `15333`, and `19166`,
   and all five hierarchies: Min4/Max8 tree of shapes, unrestricted and
   saturated residual trees with shared 8-adjacency, and unrestricted and
   saturated residual trees with centered shared 3x11 adjacency.

The notebook does not duplicate the proofs and does not implement a second
version of the constructor. It builds the final valued trees through the public
API and reconstructs their recorded event sequences for inspection.

## Running the notebook

From the repository root:

```bash
python -m pip install .
python -m pip install -r notebooks/requirements.txt
python -m jupyter lab notebooks/Self_Dual_Residual_Trees_Tutorial.ipynb
```

For a headless execution check:

```bash
mkdir -p build/notebook-runs
python -m jupyter nbconvert \
  --execute \
  --to notebook \
  --output-dir build/notebook-runs \
  notebooks/Self_Dual_Residual_Trees_Tutorial.ipynb
python scripts/validate_notebooks.py
```

## Tests and benchmark

Configure a development build with the Python bindings, tests, and benchmarks:

```bash
cmake -S . -B build \
  -DMMCFILTERS_BUILD_PYTHON=ON \
  -DMMCFILTERS_BUILD_TESTS=ON \
  -DMMCFILTERS_BUILD_BENCHMARKS=ON
cmake --build build --target \
  unit_min_max_residual_tree_factory \
  mmcfilters_min_max_residual_tree_benchmark
```

Run the focused C++ regression test:

```bash
ctest --test-dir build \
  -R unit_min_max_residual_tree_factory \
  --output-on-failure
```

Run the empirical construction benchmark on a grayscale image:

```bash
./build/benchmarks/mmcfilters_min_max_residual_tree_benchmark \
  path/to/image.png 5 both
```

The benchmark reports median construction time, node counts, and rejected
extrema. It is an empirical implementation benchmark; it does not add an
asymptotic complexity claim beyond the paper.
