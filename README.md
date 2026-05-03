# MorphologicalAttributeFilters

[![CI](https://github.com/wonderalexandre/MorphologicalAttributeFilters/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/wonderalexandre/MorphologicalAttributeFilters/actions/workflows/ci.yml)
[![Release](https://github.com/wonderalexandre/MorphologicalAttributeFilters/actions/workflows/release.yml/badge.svg)](https://github.com/wonderalexandre/MorphologicalAttributeFilters/actions/workflows/release.yml)

> Project status: transitional research library.
>
> This repository exists to support experiments around dynamic morphological
> trees, including component trees, trees of shapes, mutable tree topologies,
> converters, and incremental attribute computation. It is not intended to be a
> long-term general-purpose replacement for Higra, and this package may change
> substantially or stop existing in the future.

For general-purpose hierarchical image analysis, production workflows, and
stable APIs, prefer [Higra](https://github.com/higra/Higra). Use this library
when you specifically need dynamic tree operations, topology/proper-part
experiments, interoperability with Higra-style hierarchies, or the incremental
attribute framework provided here.

The current codebase is organised around a dense `NodeId` API, direct
proper-part ownership, incremental attribute computation, and Python bindings
for interactive work.

## Relationship with Higra

Higra should be preferred for stable, general-purpose component-tree and
hierarchical image-analysis workflows. This repository is mainly useful as an
experimental bridge for cases where the current Higra model is not the most
convenient fit, especially dynamic tree topology changes and project-specific
incremental attribute machinery.

## Current scope

This library currently concentrates code for:

- building max-trees, min-trees, and trees of shapes from images;
- representing trees with dense internal `NodeId` values;
- direct proper-part ownership and image-domain topology queries;
- dynamic `MorphologicalTree` operations and staged structural edits;
- incremental computation of geometric, statistical, topological, moment-based,
  bit-quads, and max-distance attributes;
- attribute filters, extinction values, contours, and Ultimate Attribute
  Opening;
- a C++ core that is effectively header-only, plus a pybind11 module for
  Python.

## Interoperability and converters

Although this project is not intended to replace Higra, it provides converters
that make it useful as an experimental interoperability layer.

The library can:

- build max-trees, min-trees, and trees of shapes directly from images;
- import static Higra-style hierarchies through `(parent, altitude)` arrays;
- preserve imported Higra node ids while the topology remains unedited;
- export the current live tree back to a compact Higra-style representation;
- move between topology-only `MorphologicalTree` use cases and
  `WeightedMorphologicalTree` workflows when explicit node altitudes are
  required;
- prototype dynamic tree edits and incremental attributes here, then export the
  resulting hierarchy for use elsewhere.

Imported Higra node ids are preserved only while the topology is not edited.
After edits, `exportHigraHierarchy()` creates a new compact hierarchy.

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

Installed C++ package:

```cmake
find_package(mmcfilters CONFIG REQUIRED)
target_link_libraries(my_target PRIVATE mmcfilters::core)
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

img = np.array(
    [
        [3, 3, 2, 2],
        [3, 4, 4, 2],
        [1, 4, 5, 2],
        [1, 1, 5, 0],
    ],
    dtype=np.uint8,
)

# radius=1.5 selects the 8-neighbourhood on a 2D square grid.
# Use radius=1.0 for 4-connectivity.
adjacency_radius = 1.5

# Case 1: build a topology-only max-tree for structural queries.
topology_tree = mmcfilters.MorphologicalTree.createMaxTree(
    img,
    radius=adjacency_radius,
)
root_node_id = topology_tree.getRoot()
root_children = topology_tree.getChildren(root_node_id)
root_direct_proper_parts = topology_tree.getProperParts(root_node_id)

# Case 2: inspect the component that owns one image pixel.
pixel_index = 10
pixel_component_id = topology_tree.getSmallestComponent(pixel_index)
pixel_component_pixels = list(topology_tree.connected_component_of(pixel_component_id))
pixel_component_mask = topology_tree.reconstructNode(pixel_component_id)

# Case 3: compute attributes that depend only on tree topology/support.
area_by_node = mmcfilters.Attribute.computeSingleAttribute(
    topology_tree,
    mmcfilters.Attribute.Type.AREA,
)

# Case 4: reuse the same topology as a weighted tree when node altitudes are required.
weighted_tree = mmcfilters.WeightedMorphologicalTree.createFromTopology(
    topology_tree,
    img,
)
max_dist_by_node = mmcfilters.Attribute.computeSingleAttribute(
    weighted_tree,
    mmcfilters.Attribute.Type.MAX_DIST,
)
reconstructed_image = weighted_tree.reconstructionImage()

# Case 5: export/import a Higra-style hierarchy for interoperability.
higra_parent, higra_altitude = weighted_tree.exportHigraHierarchy()
max_dist_by_higra = weighted_tree.project_node_values_to_exported_higra(
    max_dist_by_node,
    mmcfilters.Attribute.Type.MAX_DIST,
)
roundtrip_weighted_tree = mmcfilters.WeightedMorphologicalTree.createFromHigraParent(
    higra_parent,
    higra_altitude,
    img.shape[0],
    img.shape[1],
    mmcfilters.WeightedMorphologicalTree.MAX_TREE,
    radius=adjacency_radius,
)
```

## Current `MorphologicalTree` boundary

`MorphologicalTree` currently provides:

- a dynamic hierarchy over dense internal `NodeId` values;
- direct ownership of proper parts;
- structural traversal and safe local mutation operations;
- optional adjacency metadata plus explicit image-domain dimensions.

Weighted quantities are intentionally outside the topology-only tree:

- `WeightedMorphologicalTree` owns the dense node-altitude buffer;
- `WeightedMorphologicalTree` encapsulates its topology and exposes it only as `const MorphologicalTree&`;
- `WeightedMorphologicalTree.createFromTopology(...)` clones an existing topology and infers altitudes from an image, avoiding a second tree reconstruction;
- image reconstruction, node residues, and Higra `(parent, altitude)` export live on `WeightedMorphologicalTree`;
- attribute computation still runs first in the internal `MorphologicalTree` node-id space and is projected only at API boundaries when requested.

The editing API has three levels:

- read-only topology queries such as `getRoot`, `getAliveNodeIds`, `getChildren`, and `getProperParts`;
- safe public mutators, `pruneNode` and `mergeNodeIntoParent`, on both topology-only and weighted trees;
- staged structural edits through `MorphologicalTree::edit()` / `TreeEditor` and `WeightedMorphologicalTree::edit()` / `WeightedTreeEditor`.

Low-level topology rewiring is not public API. Full connected-tree validation
and weighted monotone-altitude validation run at editor `commit()` time. See
[docs/editing-api.md](docs/editing-api.md) for the detailed contract.

`NodeIdSpace::HIGRA`, `getNumHigraNodes()`, and `getHigraNodeId()` refer only to the original Higra node-id domain preserved by `createFromHigraParent`. Image-built trees and trees edited after import do not expose that domain. `exportHigraHierarchy()` always creates a new compact Higra representation of the current live rooted tree. Use `projectNodeValuesToExportedHigra()` / `project_node_values_to_exported_higra()` to project node-indexed attribute buffers to that exported compact layout without reimporting the tree; each `AttributeComputer` supplies the unit-component values for the exported leaves.

The main public C++/Python surface is centred on:

- `getRoot`
- `getAliveNodeIds`
- `getChildren`
- `getProperParts`
- `getNodeParent`
- `getPathBetweenNodes`
- `getSmallestComponent`
- `getConnectedComponent`
- `getIteratorBreadthFirstTraversal`
- `getNodeSubtree`

## Repository guide

Use this map to find the right entry point quickly:

- Core C++: [mmcfilters/](mmcfilters/)
- Python interface: [pybinds/mmcfilters.cpp](pybinds/mmcfilters.cpp) and
  [python/](python/)
- Tests: [unit-tests/](unit-tests/)
- Examples and notebooks: [examples/](examples/) and
  [notebooks/README.md](notebooks/README.md)
- Design notes: [docs/editing-api.md](docs/editing-api.md) and
  [docs/contours.md](docs/contours.md)
- Build and packaging: [CMakeLists.txt](CMakeLists.txt) and
  [pyproject.toml](pyproject.toml)

## Release process

Releases are automated by GitHub Actions. For a production release:

1. Make sure the `CI` and `Package` workflows are green on `main`.
2. Create and push a semantic version tag, for example `v1.0.1`.
3. The `Release` workflow validates that the tag matches the resolved package
   version, builds the source distribution and platform wheels, validates the
   package metadata, and attaches the artifacts to a GitHub Release.

The release wheel matrix targets Python 3.9 through 3.14 on:

- Linux manylinux x86_64;
- Windows x86_64;
- macOS arm64;
- macOS Intel x86_64.

PyPI publication is intentionally manual. Download the release artifacts from
the GitHub Release or from the `Release` workflow run, then upload them with:

```bash
python -m pip install --upgrade twine
python -m twine upload dist/*
```

Manual runs of the `Release` workflow also build downloadable artifacts without
creating a GitHub Release.
