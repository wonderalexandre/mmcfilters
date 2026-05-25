# MorphologicalAttributeFilters

[![CI](https://github.com/wonderalexandre/MorphologicalAttributeFilters/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/wonderalexandre/MorphologicalAttributeFilters/actions/workflows/ci.yml)
[![Release](https://github.com/wonderalexandre/MorphologicalAttributeFilters/actions/workflows/release.yml/badge.svg)](https://github.com/wonderalexandre/MorphologicalAttributeFilters/actions/workflows/release.yml)

> Project status: research library for dynamic morphological-tree experiments.
>
> This package focuses on a proper-part tree model for image-domain partial
> partitions, tree editing, typed altitude contracts, incremental attributes,
> contours, and project-specific morphology research. It is not intended to be a
> general-purpose replacement for Higra.

## Scope

This repository is a research implementation for morphological-tree workflows
that need direct topology ownership, staged edits, typed altitude buffers, and
project-specific attribute machinery. Its central model is a tree over
image-domain partial partitions: image pixels are explicit proper parts, and
internal morphological nodes own those proper parts directly. The implemented
models share the `MorphologicalTree` topology abstraction: rooted inclusion
topology, dense internal `NodeId` values, and explicit proper-part ownership.

Current functionality includes:

- max-tree, min-tree, tree-of-shapes, and self-dual residual tree construction;
- import/export of Higra-style `(parent, altitude)` hierarchies;
- dynamic topology edits through safe mutators and staged editors;
- gray-level, shape, boundary, topology, and max-distance attributes;
- attribute filters, extinction values, and Ultimate Attribute Opening;
- a C++20 header-oriented core plus a pybind11 Python package.

Higra remains the better fit for stable, general-purpose hierarchical
image-analysis workflows. Use this project when the experiment needs mutable
tree topology, direct owner-state access, or the local attribute/filter
machinery exposed here. See
[docs/higra-interoperability.md](docs/higra-interoperability.md) for import,
export, and attribute-projection contracts.

Python currently follows the canonical 8-bit contract: factory inputs must be
C-contiguous `np.uint8` arrays and external altitude inputs must stay in
`[0, 255]`. C++ supports typed max/min and SDRT construction through `Image<T>`,
typed Higra imports through `createFromHigraParent<T>`, and read-only
`WeightedTreeView<T>` altitude spans. Tree of Shapes construction is currently
`uint8_t`.

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

For attribute computation in C++, prefer the aggregate public header:

```cpp
#include <mmcfilters/attributes/Attributes.hpp>
```

That header exposes the stable attribute-facing contracts: `Attribute`,
`AttributeGroup`, `AttributeNames`, `AttributeComputation`, and owning result
types. Use `AttributeComputation` with `WeightedMorphologicalTree<T>` or
`WeightedTreeView<T>` for altitude-aware requests, and use
`computeTopologyAttributes(...)` for topology/support-only requests. Concrete
attribute computers are advanced implementation components; ordinary consumers
should not include `mmcfilters/*/detail/` headers directly.

Canonical C++ attribute usage:

```cpp
#include <cstdint>
#include <array>
#include <vector>

#include <mmcfilters/attributes/Attributes.hpp>
#include <mmcfilters/trees/MorphologicalTreeFactory.hpp>
#include <mmcfilters/utils/Image.hpp>

using namespace mmcfilters;

std::array<std::uint8_t, 16> pixels{
    3, 3, 2, 2,
    3, 4, 4, 2,
    1, 4, 5, 2,
    1, 1, 5, 0};
auto image = ImageUInt8::fromExternal(pixels.data(), 4, 4);
auto weighted = MorphologicalTreeFactory::createMaxTree(image, 1.5);

auto [names, values] = AttributeComputation::computeAttributes(
    weighted,
    std::vector<AttributeOrGroup>{AREA, LEVEL, MAX_DIST});

auto [topologyNames, topologyValues] =
    AttributeComputation::computeTopologyAttributes(
        weighted.topology(),
        std::vector<AttributeOrGroup>{AREA, BOX_WIDTH, BALANCE_NODE});

auto [deltaNames, deltaValues] =
    AttributeComputation::computeSingleAttributeWithDelta(
        weighted,
        LEVEL,
        AltitudeDiff<std::uint8_t>{1},
        2);
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

For detailed API contracts, use:

- tree ownership and `NodeId` semantics: [docs/trees.md](docs/trees.md);
- attribute computation and descriptor semantics:
  [docs/attributes.md](docs/attributes.md);
- editing and validation rules: [docs/editing-api.md](docs/editing-api.md);
- Higra preserved vs exported node-id spaces:
  [docs/higra-interoperability.md](docs/higra-interoperability.md);
- Python facade examples: [docs/python-api.md](docs/python-api.md).

## Repository guide

Use this map to find the right entry point quickly:

- Core C++: [mmcfilters/](mmcfilters/)
- Morphological tree model: [docs/trees.md](docs/trees.md)
- Public attributes: [docs/attributes.md](docs/attributes.md). The C++ package
  is header-only; installed `detail` headers are implementation support, not
  compatibility-contract headers.
- Attribute filters, extinction values, and UAO:
  [docs/filters.md](docs/filters.md)
- Higra interoperability:
  [docs/higra-interoperability.md](docs/higra-interoperability.md)
- Python interface: [docs/python-api.md](docs/python-api.md),
  [pybinds/mmcfilters.cpp](pybinds/mmcfilters.cpp), and [python/](python/)
- Tests: [unit-tests/](unit-tests/)
- Examples and notebooks: [examples/README.md](examples/README.md) and
  [notebooks/README.md](notebooks/README.md)
- Design notes: [docs/editing-api.md](docs/editing-api.md) and
  [docs/contours.md](docs/contours.md)
- Build and packaging: [CMakeLists.txt](CMakeLists.txt) and
  [pyproject.toml](pyproject.toml)

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
