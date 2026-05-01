# Notebooks

This directory contains exploratory notebooks that are expected to follow the
current public Python API. Treat them as examples for experiments, not as a
stable product surface.

## Preferred tree model

For new notebook work, keep the topology-only and weighted use cases separate:

- use `mmcfilters.MorphologicalTree` for topology, traversal, proper-part
  ownership, and structural edits;
- use `mmcfilters.WeightedMorphologicalTree` when the workflow needs node
  altitudes, image reconstruction, node residues, Higra `(parent, altitude)`
  export, or attributes that require intensity values;
- use `WeightedMorphologicalTree.createFromTopology(topology_tree, image)` when
  a notebook already has a topology tree and needs the weighted view without
  rebuilding the hierarchy.

Preferred topology-only calls:

- `tree.getRoot()`
- `tree.getAliveNodeIds()`
- `tree.getChildren(node_id)`
- `tree.getProperParts(node_id)`
- `tree.getConnectedComponent(node_id)`
- `tree.getSmallestComponent(pixel_id)`
- `tree.reconstructNode(node_id)`

Preferred weighted calls:

- `weighted.getAltitude(node_id)`
- `weighted.reconstructionImage()`
- `weighted.exportHigraHierarchy()`
- `weighted.project_node_values_to_exported_higra(node_values, attribute)`

Preferred builders:

- `mmcfilters.MorphologicalTree.createMaxTree(image, radius=1.5)`
- `mmcfilters.MorphologicalTree.createMinTree(image, radius=1.5)`
- `mmcfilters.MorphologicalTree.createTreeOfShapes(image)`
- `mmcfilters.WeightedMorphologicalTree.createMaxTree(image, radius=1.5)`
- `mmcfilters.WeightedMorphologicalTree.createMinTree(image, radius=1.5)`
- `mmcfilters.WeightedMorphologicalTree.createTreeOfShapes(image)`

Preferred attribute calls:

- `mmcfilters.Attribute.computeSingleAttribute(tree, mmcfilters.Attribute.AREA)`
- `mmcfilters.Attribute.computeAttributes(tree, [mmcfilters.Attribute.AREA, ...])`
- `mmcfilters.Attribute.computeSingleAttribute(weighted, mmcfilters.Attribute.MAX_DIST)`

Use tree-type constants instead of numeric literals when importing Higra-style
hierarchies:

- `mmcfilters.MorphologicalTree.MAX_TREE`
- `mmcfilters.MorphologicalTree.MIN_TREE`
- `mmcfilters.MorphologicalTree.TREE_OF_SHAPES`
- `mmcfilters.WeightedMorphologicalTree.MAX_TREE`
- `mmcfilters.WeightedMorphologicalTree.MIN_TREE`
- `mmcfilters.WeightedMorphologicalTree.TREE_OF_SHAPES`

## Current notebooks

The maintained notebooks in this directory are:

- `Filter.ipynb`
- `Attribute_Filters.ipynb`
- `UAO_Exemplos.ipynb`

Other notebooks may remain useful as exploratory references, but new notebooks
should follow the API above and should not depend on removed legacy interfaces
such as `NodeMT`, `ResidualTree`, or old parent-array helper layers.
