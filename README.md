# mmcfilters

[![PyPI](https://img.shields.io/pypi/v/mmcfilters.svg)](https://pypi.org/project/mmcfilters/)
[![Python](https://img.shields.io/pypi/pyversions/mmcfilters.svg)](https://pypi.org/project/mmcfilters/)
[![CI](https://github.com/wonderalexandre/mmcfilters/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/wonderalexandre/mmcfilters/actions/workflows/ci.yml)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue.svg)](https://wonderalexandre.github.io/mmcfilters/)
[![License: GPL-3.0](https://img.shields.io/badge/license-GPL--3.0--only-blue.svg)](https://github.com/wonderalexandre/mmcfilters/blob/main/LICENSE)

`mmcfilters` is an alpha research library for constructing and editing
morphological hierarchies, computing node attributes, and applying connected
filters and hierarchy operators. It provides a C++20 core and Python bindings.

The library represents hierarchies through a rooted connected-subset tree
model. A morphological tree of partial partitions is the stricter case in
which every node has a non-empty proper part. All supported image-domain trees
use the same topology model, with an explicit smallest-node map and optional
regular-grid geometry and adjacency.

Use `mmcfilters` for experiments that require mutable tree topology, direct
ownership of hierarchy state, incremental attributes, residual hierarchies, or
the filtering operators implemented by this project. For stable,
general-purpose hierarchical image processing, [Higra](https://github.com/higra/Higra)
will often be the better choice.

## Main capabilities

- construction of max-trees, min-trees, trees of shapes, and self-dual
  residual trees;
- generic rooted hierarchies with safe local and staged edits;
- a broad catalog of node attributes organized into gray-level, shape, moment,
  boundary, and tree-topology groups;
- a range of attribute-based filters and hierarchy operators, including pruning
  and non-pruning strategies;
- pixel contours and ordered geometric contour traces;
- import, export, and attribute projection for Higra-style hierarchies;
- typed altitude handling in C++ and a focused `np.uint8` Python API.

## Requirements

- Python 3.9–3.14 for the distributed Python package;
- NumPy 1.23 or newer;
- a C++20 compiler and CMake 3.20 or newer for source builds.

## Installation

Install the Python package from PyPI:

```bash
python -m pip install mmcfilters
```

Install from a source checkout:

```bash
python -m pip install .
```

For a C++-only build and installation:

```bash
cmake -S . -B build -DMMCFILTERS_BUILD_PYTHON=OFF
cmake --build build
cmake --install build --prefix /path/to/prefix
```

The installed C++ package exports the `mmcfilters::core` target:

```cmake
find_package(mmcfilters CONFIG REQUIRED)
target_link_libraries(my_target PRIVATE mmcfilters::core)
```

## Quick start

```python
import numpy as np
import mmcfilters

image = np.ascontiguousarray(
    [
        [3, 3, 2, 2],
        [3, 4, 4, 2],
        [1, 4, 5, 2],
        [1, 1, 5, 0],
    ],
    dtype=np.uint8,
)

# radius=1.5 selects 8-connectivity; use 1.0 for 4-connectivity.
tree = mmcfilters.MorphologicalTreeFactory.create_max_tree(image, radius=1.5)

root = tree.root
root_children = tree.children(root)

pixel = 10
smallest = tree.smallest_node(pixel)
support_pixels = list(tree.node_support(smallest))
component_mask = tree.reconstruct_node(smallest)

area = mmcfilters.Attribute.compute_single_topology_attribute(
    tree,
    mmcfilters.Attribute.AREA,
)
max_dist = mmcfilters.Attribute.compute_single_attribute(
    tree,
    mmcfilters.Attribute.MAX_DIST,
)
reconstructed_image = tree.reconstruct_from_node_altitudes()
```

## Documentation

### User guides

| Goal | Guide |
| --- | --- |
| Understand the tree model and construct hierarchies | [Morphological trees](https://wonderalexandre.github.io/mmcfilters/md_docs_2trees.html) |
| Use the Python interface | [Python API](https://wonderalexandre.github.io/mmcfilters/md_docs_2python-api.html) |
| Compute attributes | [Attributes](https://wonderalexandre.github.io/mmcfilters/md_docs_2attributes.html) |
| Apply filtering and hierarchy operators | [Filters](https://wonderalexandre.github.io/mmcfilters/md_docs_2filters.html) |
| Edit a tree safely | [Editing API](https://wonderalexandre.github.io/mmcfilters/md_docs_2editing-api.html) |
| Extract pixel contours or geometric traces | [Pixel contours](https://wonderalexandre.github.io/mmcfilters/md_docs_2contours.html) and [contour traces](https://wonderalexandre.github.io/mmcfilters/md_docs_2contour-traces.html) |
| Import from or export to Higra | [Higra interoperability](https://wonderalexandre.github.io/mmcfilters/md_docs_2higra-interoperability.html) |

### References

- [Attribute catalog](https://wonderalexandre.github.io/mmcfilters/md_docs_2attribute-catalog.html)
- [Generated C++ API](https://wonderalexandre.github.io/mmcfilters/annotated.html)

### Contributor guides

- [Attribute computer architecture](https://wonderalexandre.github.io/mmcfilters/md_docs_2attribute-computer-architecture.html)
- [Scientific benchmark builds](https://wonderalexandre.github.io/mmcfilters/md_docs_2scientific-benchmark-builds.html)

## License

mmcfilters is distributed under the
[GNU General Public License v3.0](https://github.com/wonderalexandre/mmcfilters/blob/main/LICENSE).
