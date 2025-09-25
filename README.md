# MorphologicalAttributeFilters

MorphologicalAttributeFilters is a C++/Python library for connected image filtering based on morphological trees (component
trees and trees of shapes). The code provides high-performance implementations together with Python bindings for interactive
experimentation.

## Key features

* Construction of morphological trees (component tree, tree of shapes) with different connectivities.
* Incremental computation of geometric, topological and radiometric attributes.
* Attribute-based filters (direct and subtractive rules, pruning, openings, ultimate attribute opening, etc.).
* Utilities for extinction values, primitive families, MSER and Bit-Quads.
* Pybind11 bindings that expose the high-level operations to Python.

## Installation

```bash
pip install mmcfilters
```

## Quick example

```python
import numpy as np
import mmcfilters
Type = mmcfilters.Attribute.Type

img = np.random.randint(0, 255, size=(128, 128), dtype=np.uint8)

maxtree = mmcfilters.MorphologicalTree(img, True, 1.5)
#mintree = mmcfilters.MorphologicalTree(img, False, 1.5)
#tos = mmcfilters.MorphologicalTree(img) 

filter = mmcfilters.AttributeFilters(maxtree)
area = mmcfilters.Attribute.computeSingleAttribute(maxtree, Type.AREA)

img_filtered = filter.filteringDirectRule(area > 50)
```
