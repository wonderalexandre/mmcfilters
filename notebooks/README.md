# Notebooks

This directory contains exploratory notebooks that are expected to follow the
current public Python API. Treat them as examples for experiments, not as a
stable product surface.

For the full Python contract and short runnable examples, see
[`docs/python-api.md`](../docs/python-api.md).

## Preferred tree model

For new notebook work, build trees through `mmcfilters.MorphologicalTreeFactory`.
Python currently accepts C-contiguous `np.uint8` images and returns the canonical
`WeightedMorphologicalTree` binding for image, Higra, and tree-of-shapes
construction. Use that object for topology queries, altitude access,
reconstruction, Higra export, and exported-Higra attribute projection.

Preferred builders:

- `mmcfilters.MorphologicalTreeFactory.createMaxTree(image, radius=1.5)`
- `mmcfilters.MorphologicalTreeFactory.createMinTree(image, radius=1.5)`
- `mmcfilters.MorphologicalTreeFactory.createTreeOfShapes(image)`

Preferred attribute calls:

- `mmcfilters.Attribute.computeTopologyAttributes(tree, [...])` for
  topology/support-only descriptors such as `AREA` and bounding boxes;
- `mmcfilters.Attribute.computeSingleAttribute(weighted, attr)` and
  `mmcfilters.Attribute.computeAttributes(weighted, [...])` for ordinary
  weighted descriptors.

Use `MorphologicalTreeKind` instead of numeric literals when importing
Higra-style hierarchies:

- `mmcfilters.MorphologicalTreeKind.MAX_TREE`
- `mmcfilters.MorphologicalTreeKind.MIN_TREE`
- `mmcfilters.MorphologicalTreeKind.TREE_OF_SHAPES`

Descriptor semantics and supported tree kinds are documented in
[`docs/attributes.md`](../docs/attributes.md).

## Current notebooks

The maintained notebooks in this directory are:

- `Attribute_Filters.ipynb`
- `Filter.ipynb`
- `Higra_MaxDist_Filtering.ipynb`
- `MaxDistExample.ipynb`
- `SimpleExamples.ipynb`
- `ToS_Contour_Example.ipynb`
- `UAO_Examples.ipynb`
