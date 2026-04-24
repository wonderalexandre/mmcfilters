# Notebooks

This directory contains notebooks that are expected to follow the current public API of the project.

## Preferred API

For new notebook work, prefer:

- dense internal-node ids (`NodeId`);
- `tree.getRoot()`
- `tree.getAliveNodeIds()`
- `tree.getChildren(nodeId)`
- `tree.getProperParts(nodeId)`
- `tree.getSmallestComponent(pixelId)`
- `tree.getAltitude(nodeId)`
- `tree.reconstructionImage()`
- `mmcfilters.Attribute.computeSingleAttribute(tree, mmcfilters.Attribute.AREA)`
- `mmcfilters.Attribute.computeAttributes(tree, [...])`
- `mmcfilters.Attribute.traversePostOrder(...)`
- `mmcfilters.MorphologicalTree(img, interpolation=mmcfilters.ToSInterpolation.SelfDual)` for Tree of Shapes construction

## Current notebooks

The maintained notebooks in this directory are:

- `Filter.ipynb`
- `Attribute_Filters.ipynb`
- `UAO_Exemplos.ipynb`

New notebooks should follow the same API and should not depend on removed legacy interfaces such as `NodeMT`, `ResidualTree`, or old parent-array helper layers.
