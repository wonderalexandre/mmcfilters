# MorphologicalAttributeFilters

MorphologicalAttributeFilters is a C++/Python library for connected image filtering based on morphological trees, including component trees and trees of shapes.

The current codebase is organised around a dense `NodeId` API, direct proper-part ownership, incremental attribute computation, and Python bindings for interactive work.

## Main capabilities

- construction of max-trees, min-trees, and trees of shapes;
- dynamic `MorphologicalTree` operations on dense internal-node ids;
- incremental computation of geometric, statistical, topological, moment-based, bit-quads, and max-distance attributes;
- attribute filters, extinction values, contours, and Ultimate Attribute Opening;
- a C++ core that is effectively header-only, plus a pybind11 module for Python.

## Installation

From PyPI:

```bash
pip install mmcfilters
```

From source:

```bash
cmake -S . -B build -DMMCFILTERS_BUILD_PYTHON=ON
cmake --build build
```

To enable the regression suite or examples:

```bash
cmake -S . -B build \
  -DMMCFILTERS_BUILD_PYTHON=ON \
  -DMMCFILTERS_BUILD_TESTS=ON \
  -DMMCFILTERS_BUILD_EXAMPLES=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

## Quick Python example

```python
import numpy as np
import mmcfilters

img = np.random.randint(0, 255, size=(128, 128), dtype=np.uint8)

tree = mmcfilters.MorphologicalTree(img, True, 1.5)
weighted = mmcfilters.WeightedMorphologicalTree(img, True, 1.5)
root_id = tree.getRoot()
alive_ids = tree.getAliveNodeIds()
children_of_root = tree.getChildren(root_id)
proper_parts_of_root = tree.getProperParts(root_id)

component_id = tree.getSmallestComponent(0)
component_mask = tree.reconstructNode(component_id)
reconstructed_image = weighted.reconstructionImage()
area = mmcfilters.Attribute.computeSingleAttribute(tree, mmcfilters.Attribute.AREA)

visited = []
mmcfilters.Attribute.traversePostOrder(
    tree,
    lambda node_id: visited.append(("pre", node_id)),
    lambda parent_id, child_id: visited.append(("merge", parent_id, child_id)),
    lambda node_id: visited.append(("post", node_id)),
)
```

## Current `MorphologicalTree` boundary

`MorphologicalTree` currently provides:

- a dynamic hierarchy over dense internal `NodeId` values;
- direct ownership of proper parts;
- structural traversal and mutation operations;
- optional adjacency metadata plus explicit image-domain dimensions.

Weighted quantities are intentionally outside the topology-only tree:

- `WeightedMorphologicalTree` owns the dense node-altitude buffer;
- image reconstruction, node residues, and Higra `(parent, altitude)` export live on `WeightedMorphologicalTree`;
- attribute computation still runs first in the internal `MorphologicalTree` node-id space and is projected only at API boundaries when requested.

The main public C++/Python surface is centred on:

- `getRoot`
- `getAliveNodeIds`
- `getChildren`
- `getProperParts`
- `getNodeParent`
- `getPathBetweenNodes`
- `getSmallestComponent`
- `getIteratorBreadthFirstTraversal`
- `getNodeSubtree`

## Breaking Changes

Recent refactors changed the ownership model of incremental attribute results:

- `AttributeComputedIncrementally::{computeSingleAttribute, computeAttributes, computeSingleAttributeWithDelta}` now return move-only results whose buffers are owned by `std::vector<float>`.
- `DependencyMap` no longer owns cached buffers. It stores non-owning views, so caller-provided dependencies must be seeded with `computed.view()` and the owning result must stay alive during the computation.
- Python still returns mutable `numpy.ndarray` objects, but the NumPy capsule now owns the transferred C++ buffer instead of sharing a `std::shared_ptr<float[]>`.

Typical C++ migration pattern:

```cpp
auto computed = AttributeComputedIncrementally::computeSingleAttribute(tree, AREA);
DependencyMap deps;
deps[AREA] = computed.view();

auto combined = AttributeComputedIncrementally::computeAttributes(
    tree,
    {AREA, VOLUME},
    deps
);
```

## Repository guide

- [unit-tests/README.md](/Users/wonderalexandre/GitHub/MorphologicalAttributeFilters/unit-tests/README.md)
  Official regression suite for the active API.
- [notebooks/README.md](/Users/wonderalexandre/GitHub/MorphologicalAttributeFilters/notebooks/README.md)
  Notes for the maintained notebooks.
- [examples/README.md](/Users/wonderalexandre/GitHub/MorphologicalAttributeFilters/examples/README.md)
  Standalone examples, including the `EdtDIFT` PNG export example.
