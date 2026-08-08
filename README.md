# MorphologicalAttributeFilters

[![CI](https://github.com/wonderalexandre/MorphologicalAttributeFilters/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/wonderalexandre/MorphologicalAttributeFilters/actions/workflows/ci.yml)
[![Release](https://github.com/wonderalexandre/MorphologicalAttributeFilters/actions/workflows/release.yml/badge.svg)](https://github.com/wonderalexandre/MorphologicalAttributeFilters/actions/workflows/release.yml)

> Project status: research library for dynamic morphological-tree experiments.
>
> This package focuses on a proper-part tree model for finite partial
> partitions, optional image-domain geometry, tree editing, typed altitude contracts, incremental attributes,
> contours, and project-specific morphology research. It is not intended to be a
> general-purpose replacement for Higra.

## Scope

This repository is a research implementation for morphological-tree workflows
that need direct topology ownership, staged edits, typed altitude buffers, and
project-specific attribute machinery. Its central model is a tree over finite
partial partitions: proper parts have their own id domain, and internal
morphological nodes own them directly. Image trees optionally attach a regular
2D layout, but that layout is not required by the core topology. The implemented
models share the `MorphologicalTree` topology abstraction: rooted inclusion
topology, dense internal `NodeId` values, and explicit proper-part ownership.
Nodes and proper parts have independent id domains, so a node may have no
direct proper part when its support is supplied entirely by descendants.
Committed nodes must still have non-empty subtree support. Interpretation is
expressed through generic altitude-order and adjacency capabilities;
`MorphologicalTreeKind` is a descriptive diagnostic label.

Current functionality includes:

- max-tree, min-tree, and tree-of-shapes construction;
- import/export of Higra-style `(parent, altitude)` hierarchies;
- dynamic topology edits through safe mutators and staged editors;
- gray-level, shape, boundary, topology, and max-distance attributes;
- attribute filters, extinction values, hierarchy saliency maps, and Ultimate
  Attribute Opening;
- a C++20 header-oriented core plus a pybind11 Python package.

Higra remains the better fit for stable, general-purpose hierarchical
image-analysis workflows. Use this project when the experiment needs mutable
tree topology, direct owner-state access, or the local attribute/filter
machinery exposed here. See
[docs/attribute-catalog.md](docs/attribute-catalog.md) for the public descriptor
catalog and [docs/higra-interoperability.md](docs/higra-interoperability.md)
for import, export, and attribute-projection contracts.
The scientific contracts and operator distinctions for saliency are centralized
in [docs/saliency.md](docs/saliency.md). The unreleased API notes are recorded in
[CHANGELOG.md](CHANGELOG.md): `computeFormalSaliencyEdgeMap` keeps its name but
now uses the Cousty persistence MST/BPTAO construction. The former direct LCA
behavior is available only through `computeMonotoneExtinctionProjection`.
For an executable English introduction—from a hand-computable `3x3` example
through LCA, MST, extinction, contour visualization, and shape-space—see
[notebooks/Saliency_Maps_Tutorial.ipynb](notebooks/Saliency_Maps_Tutorial.ipynb).

Python currently follows the canonical 8-bit contract: factory inputs must be
C-contiguous `np.uint8` arrays and external altitude inputs must stay in
`[0, 255]`. C++ supports typed max/min construction through `Image<T>`,
typed Higra imports through `createFromHigraParent<T>`, and read-only
`WeightedTreeView<T>` altitude spans. Tree of Shapes construction is currently
`uint8_t`.

Self-dual residual trees are available through synchronized max-tree/min-tree
construction. The public factory exposes both the unrestricted hierarchy and
the version restricted to saturated regional extrema under one shared symmetric
adjacency. The underlying algorithm and its scientific reference backends are
developed in the `MorphoTreeDynamics` project. This branch integrates the
current synchronized production builder; independent reference constructions
are confined to differential tests. The obsolete historical implementation
previously removed from this repository has not been restored.

## Installation

From PyPI:

```bash
pip install mmcfilters
```

From source:

```bash
python -m pip install .
```

For a C++ source build and installation:

```bash
cmake -S . -B build -DMMCFILTERS_BUILD_PYTHON=OFF
cmake --build build
cmake --install build --prefix /path/to/prefix
```

Consume the installed C++ package with:

```cmake
find_package(mmcfilters CONFIG REQUIRED)
target_link_libraries(my_target PRIVATE mmcfilters::core)
```

In a repository checkout, enable the regression suite or C++ examples with:

```bash
cmake -S . -B build \
  -DMMCFILTERS_BUILD_PYTHON=ON \
  -DMMCFILTERS_BUILD_TESTS=ON \
  -DMMCFILTERS_BUILD_EXAMPLES=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

Internal benchmarks can be built with
`-DMMCFILTERS_BUILD_BENCHMARKS=ON`.

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

# Python factories require C-contiguous np.uint8 images. Use C++ for typed
# int32/float altitude trees.
img = np.ascontiguousarray(img, dtype=np.uint8)

# radius=1.5 selects the 8-neighbourhood on a 2D square grid.
# Use radius=1.0 for 4-connectivity.
adjacency_radius = 1.5

# Case 1: build a weighted max-tree. Topology queries are available on it.
weighted_tree = mmcfilters.MorphologicalTreeFactory.createMaxTree(
    img,
    radius=adjacency_radius,
)
root_node_id = weighted_tree.getRoot()
root_children = weighted_tree.getChildren(root_node_id)
root_direct_proper_parts = weighted_tree.getProperParts(root_node_id)

# Case 2: inspect the component that owns one image pixel.
pixel_index = 10
pixel_component_id = weighted_tree.getProperPartOwner(pixel_index)
pixel_component_pixels = list(weighted_tree.getConnectedComponent(pixel_component_id))
pixel_component_mask = weighted_tree.reconstructNode(pixel_component_id)

# Case 3: compute attributes that depend only on tree topology/support.
topology_names, topology_by_node = mmcfilters.Attribute.computeTopologyAttributes(
    weighted_tree,
    [mmcfilters.Attribute.AREA, mmcfilters.Attribute.BOX_HEIGHT],
)
area_by_node = topology_by_node[:, topology_names["AREA"]]
box_height_by_node = topology_by_node[:, topology_names["BOX_HEIGHT"]]

# Case 4: compute altitude-dependent attributes and reconstruct the image.
max_dist_by_node = mmcfilters.Attribute.computeSingleAttribute(
    weighted_tree,
    mmcfilters.Attribute.MAX_DIST,
)
reconstructed_image = weighted_tree.reconstructionImage()

# Case 5: run Ultimate Attribute Opening, the public UAO API.
uao = mmcfilters.UltimateAttributeOpening(weighted_tree, box_height_by_node)
uao.execute(img.shape[0])
max_contrast_image = uao.getMaxContrastImage()
associated_image = uao.getAssociatedImage()

# Case 6: export/import a Higra-style hierarchy for interoperability.
higra_parent, higra_altitude = weighted_tree.exportHigraHierarchy()
max_dist_by_higra = weighted_tree.project_node_values_to_exported_higra(
    max_dist_by_node,
    mmcfilters.Attribute.MAX_DIST,
)
roundtrip_weighted_tree = mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
    higra_parent,
    higra_altitude,
    img.shape[0],
    img.shape[1],
    mmcfilters.MorphologicalTreeKind.MAX_TREE,
    radius=adjacency_radius,
)
```

## Repository guide

Use this map to find the right entry point quickly:

API guides:

- Morphological tree model: [docs/trees.md](docs/trees.md)
- Hierarchy saliency maps and scientific contracts: [docs/saliency.md](docs/saliency.md)
- Attribute computation: [docs/attributes.md](docs/attributes.md)
- Attribute catalog: [docs/attribute-catalog.md](docs/attribute-catalog.md)
- Attribute filters, extinction values, and UAO: [docs/filters.md](docs/filters.md)
- Editing API and derived-state lifetime: [docs/editing-api.md](docs/editing-api.md)
- Higra interoperability: [docs/higra-interoperability.md](docs/higra-interoperability.md)
- Python interface: [docs/python-api.md](docs/python-api.md)
- Incremental pixel contours: [docs/contours.md](docs/contours.md)
- Geometric contour traces: [docs/contour-traces.md](docs/contour-traces.md)

Contributor design notes:

- Attribute computer architecture and extension:
  [docs/attribute-computer-architecture.md](docs/attribute-computer-architecture.md)
- Contour internals and benchmarks:
  [docs/contours.md](docs/contours.md)
- Recorded benchmark artifacts and provenance status:
  [docs/benchmarks/index.md](docs/benchmarks/index.md)

## Documentation

The `Documentation` workflow validates the public and internal Doxygen targets.
On pushes to `main`, it publishes only the public HTML output to GitHub Pages:
[wonderalexandre.github.io/MorphologicalAttributeFilters](https://wonderalexandre.github.io/MorphologicalAttributeFilters/).
The internal HTML output remains available as a workflow artifact.

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
