#!/usr/bin/env python3

import importlib.util
import pathlib
import sys

import numpy as np


def load_native_module(build_dir: pathlib.Path):
    package_dir = build_dir / "python" / "mmcfilters"
    package_init = package_dir / "__init__.py"
    if not package_init.is_file():
        raise RuntimeError(f"package init not found: {package_init}")

    for name in list(sys.modules):
        if name == "mmcfilters" or name.startswith("mmcfilters."):
            sys.modules.pop(name, None)

    def handles_mmcfilters(finder):
        known_modules = {}
        known_modules.update(getattr(finder, "known_source_files", {}))
        known_modules.update(getattr(finder, "known_wheel_files", {}))
        return any(name == "mmcfilters" or name.startswith("mmcfilters.") for name in known_modules)

    sys.meta_path = [finder for finder in sys.meta_path if not handles_mmcfilters(finder)]

    spec = importlib.util.spec_from_file_location(
        "mmcfilters",
        package_init,
        submodule_search_locations=[str(package_dir)],
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules["mmcfilters"] = module
    spec.loader.exec_module(module)
    return module


def require(condition: bool, message: str):
    if not condition:
        raise RuntimeError(message)


def require_raises(fn, message: str):
    try:
        fn()
    except Exception:
        return
    raise RuntimeError(message)


def build_higra_hierarchy(tree_like):
    higra_parent = [
        tree_like.getSmallestComponent(pixel_id) + tree_like.numTotalProperParts
        for pixel_id in range(tree_like.numTotalProperParts)
    ]
    higra_parent.extend(
        tree_like.getNodeParent(node_id) + tree_like.numTotalProperParts
        for node_id in tree_like.getAliveNodeIds()
    )

    higra_altitude = [0] * len(higra_parent)
    for pixel_id in range(tree_like.numTotalProperParts):
        higra_altitude[pixel_id] = tree_like.getAltitude(tree_like.getSmallestComponent(pixel_id))
    for node_id in tree_like.getAliveNodeIds():
        higra_altitude[tree_like.numTotalProperParts + node_id] = tree_like.getAltitude(node_id)

    return higra_parent, higra_altitude


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: test_python_nodeid_api.py <build-dir>")

    build_dir = pathlib.Path(sys.argv[1]).resolve()
    mmcfilters = load_native_module(build_dir)
    require(hasattr(mmcfilters, "__version__"), "package import must expose __version__")
    require(hasattr(mmcfilters, "WeightedMorphologicalTree"), "package import must expose WeightedMorphologicalTree")
    require(mmcfilters.MorphologicalTree.MAX_TREE == 0, "MorphologicalTree MAX_TREE constant")
    require(mmcfilters.MorphologicalTree.MIN_TREE == 1, "MorphologicalTree MIN_TREE constant")
    require(mmcfilters.MorphologicalTree.TREE_OF_SHAPES == 2, "MorphologicalTree TREE_OF_SHAPES constant")
    require(
        mmcfilters.WeightedMorphologicalTree.MAX_TREE == mmcfilters.MorphologicalTree.MAX_TREE,
        "WeightedMorphologicalTree MAX_TREE constant",
    )
    require(
        mmcfilters.WeightedMorphologicalTree.MIN_TREE == mmcfilters.MorphologicalTree.MIN_TREE,
        "WeightedMorphologicalTree MIN_TREE constant",
    )
    require(
        mmcfilters.WeightedMorphologicalTree.TREE_OF_SHAPES == mmcfilters.MorphologicalTree.TREE_OF_SHAPES,
        "WeightedMorphologicalTree TREE_OF_SHAPES constant",
    )

    image = np.array(
        [
            [3, 3, 2, 2],
            [3, 4, 4, 2],
            [1, 4, 5, 2],
            [1, 1, 5, 0],
        ],
        dtype=np.uint8,
    )

    tree = mmcfilters.MorphologicalTree.createComponentTree(image, True)
    weighted = mmcfilters.WeightedMorphologicalTree.createComponentTree(image, True)
    max_tree = mmcfilters.MorphologicalTree.createMaxTree(image)
    min_tree = mmcfilters.MorphologicalTree.createMinTree(image)
    weighted_max_tree = mmcfilters.WeightedMorphologicalTree.createMaxTree(image)
    weighted_min_tree = mmcfilters.WeightedMorphologicalTree.createMinTree(image)

    require(max_tree.getAliveNodeIds() == tree.getAliveNodeIds(), "createMaxTree must match createComponentTree(..., True)")
    require(max_tree.treeType == mmcfilters.MorphologicalTree.MAX_TREE, "createMaxTree tree type")
    require(min_tree.treeType == mmcfilters.MorphologicalTree.MIN_TREE, "createMinTree tree type")
    require(weighted_max_tree.reconstructionImage().tolist() == weighted.reconstructionImage().tolist(), "weighted createMaxTree reconstruction")
    require(weighted_max_tree.treeType == mmcfilters.WeightedMorphologicalTree.MAX_TREE, "weighted createMaxTree tree type")
    require(weighted_min_tree.treeType == mmcfilters.WeightedMorphologicalTree.MIN_TREE, "weighted createMinTree tree type")

    weighted_from_topology = mmcfilters.WeightedMorphologicalTree.createFromTopology(tree, image)
    weighted_from_topology_snake = mmcfilters.WeightedMorphologicalTree.create_from_topology(tree, image)
    require(weighted_from_topology.getAliveNodeIds() == tree.getAliveNodeIds(), "createFromTopology must clone the topology structure")
    require(weighted_from_topology.treeType == tree.treeType, "createFromTopology tree type")
    require(weighted_from_topology.altitude == weighted.altitude, "createFromTopology altitude buffer")
    require(weighted_from_topology.reconstructionImage().tolist() == image.tolist(), "createFromTopology reconstruction")
    require(weighted_from_topology_snake.altitude == weighted.altitude, "create_from_topology alias altitude buffer")
    require_raises(
        lambda: mmcfilters.WeightedMorphologicalTree.createFromTopology(tree, np.array([[1]], dtype=np.uint8)),
        "createFromTopology must reject mismatched image domains",
    )
    weighted_from_topology.mergeNodeIntoParent(5)
    require(tree.getAliveNodeIds() == [0, 1, 2, 3, 4, 5], "createFromTopology must not mutate the source topology")

    require(tree.getRoot() == 0, "getRoot")
    require(tree.root == tree.getRoot(), "root property")
    require(tree.hasAdjacencyRelation is True, "max-tree should expose adjacency relation context")
    require(tree.getNodeParent(0) == 0, "root parent must point to itself")
    require(tree.getAliveNodeIds() == [0, 1, 2, 3, 4, 5], "alive NodeIds")
    require(tree.aliveNodeIds == tree.getAliveNodeIds(), "aliveNodeIds property")
    require(tree.alive_node_ids == tree.getAliveNodeIds(), "alive_node_ids property")
    contours = mmcfilters.ContourComputation.extraction(tree)
    require(contours.isMaterialized is False, "contours must start without global materialization in Python")
    leaf_id = tree.getLeafNodeIds()[0]
    leaf_contour_before_root_materialization = list(contours.getContour(leaf_id))
    require(contours.isContourMaterialized(leaf_id) is True, "getContour iteration must cache the requested leaf in Python")
    require(contours.isMaterialized is False, "getContour iteration must not materialize all contours in Python")
    list(contours.getContour(tree.getRoot()))
    require(contours.isMaterialized is True, "root getContour iteration must materialize all contours in Python")
    contours.materializeAll()
    require(contours.isMaterialized is True, "materializeAll must materialize all contours in Python")
    require(sorted(leaf_contour_before_root_materialization) == sorted(list(contours.getContour(leaf_id))), "Python incremental contour read must match materialized contour")
    require(tree.getLeafNodeIds() == [5], "leaf NodeIds")
    require(tree.leafNodeIds == tree.getLeafNodeIds(), "leafNodeIds property")
    require(tree.leaf_node_ids == tree.getLeafNodeIds(), "leaf_node_ids property")
    require(tree.getChildren(3) == [4], "children by NodeId")
    require(tree.childrenOf(3) == tree.getChildren(3), "childrenOf alias")
    require(tree.children_of(3) == tree.getChildren(3), "children_of alias")
    require(int(mmcfilters.Attribute.computeSingleAttribute(tree, mmcfilters.Attribute.AREA)[3]) == 8, "node area by attribute computer")
    require(mmcfilters.Attribute.describe(mmcfilters.Attribute.AREA).startswith("Area:"), "attribute description by pybind")
    require(tree.getNodeNumDescendants(2) == 3, "descendants count by NodeId")
    require(tree.getNodeNumSiblings(4) == 0, "siblings count by NodeId")
    require(tree.getNumProperParts(3) == 3, "direct proper-part count by NodeId")
    require(tree.getProperParts(3) == [0, 1, 4], "direct proper parts by NodeId")
    require(tree.properPartsOf(3) == tree.getProperParts(3), "properPartsOf alias")
    require(tree.proper_parts_of(3) == tree.getProperParts(3), "proper_parts_of alias")
    require(list(tree.getConnectedComponent(3)) == [0, 1, 4, 5, 6, 9, 10, 14], "connected component iterator by NodeId")
    require(list(tree.connectedComponentOf(3)) == list(tree.getConnectedComponent(3)), "connectedComponentOf alias")
    require(list(tree.connected_component_of(3)) == list(tree.getConnectedComponent(3)), "connected_component_of alias")
    require(list(weighted.connected_component_of(3)) == list(tree.getConnectedComponent(3)), "weighted connected_component_of alias")
    require(tree.parentOf(3) == tree.getNodeParent(3), "parentOf alias")
    require(tree.parent_of(3) == tree.getNodeParent(3), "parent_of alias")
    require(tree.smallestComponentOf(0) == tree.getSmallestComponent(0), "smallestComponentOf alias")
    require(tree.smallest_component_of(0) == tree.getSmallestComponent(0), "smallest_component_of alias")
    require(tree.nodeSubtreeOf(3) == tree.getNodeSubtree(3), "nodeSubtreeOf alias")
    require(tree.node_subtree_of(3) == tree.getNodeSubtree(3), "node_subtree_of alias")
    require(tree.descendantsOf(2) == tree.getDescendants(2), "descendantsOf alias")
    require(tree.descendants_of(2) == tree.getDescendants(2), "descendants_of alias")
    require(tree.reconstructNode(3).shape == (4, 4), "node reconstruction shape by NodeId")
    require(
        tree.reconstructNode(3).tolist()
        == [
            [255, 255, 0, 0],
            [255, 255, 255, 0],
            [0, 255, 255, 0],
            [0, 0, 255, 0],
        ],
        "node reconstruction values by NodeId",
    )
    require(int(tree.reconstructNode(3).sum()) == 255 * 8, "node reconstruction by NodeId")
    require(weighted.getRoot() == tree.getRoot(), "weighted getRoot")
    require(weighted.root == weighted.getRoot(), "weighted root property")
    require(weighted.getAliveNodeIds() == tree.getAliveNodeIds(), "weighted alive NodeIds")
    require(weighted.alive_node_ids == weighted.getAliveNodeIds(), "weighted alive_node_ids property")
    require(weighted.getAltitude(0) == 0, "weighted root altitude")
    require(weighted.getAltitude(5) == 5, "weighted getAltitude")
    require(weighted.altitudeOf(5) == weighted.getAltitude(5), "weighted altitudeOf alias")
    require(weighted.altitude_of(5) == weighted.getAltitude(5), "weighted altitude_of alias")
    require(weighted.getNodeResidue(5) == 1, "weighted node residue")
    require(weighted.residueOf(5) == weighted.getNodeResidue(5), "weighted residueOf alias")
    require(weighted.residue_of(5) == weighted.getNodeResidue(5), "weighted residue_of alias")
    require_raises(lambda: tree.getNodeParent(-1), "invalid getNodeParent must throw")
    require_raises(lambda: tree.getChildren(-1), "invalid getChildren must throw")
    require_raises(lambda: tree.children_of(-1), "invalid children_of must throw")
    require_raises(lambda: tree.getNumChildren(999), "invalid getNumChildren must throw")
    require_raises(lambda: tree.getNodeTimePreOrder(999), "invalid getNodeTimePreOrder must throw")
    require_raises(lambda: tree.getProperParts(999), "invalid getProperParts must throw")
    require_raises(lambda: list(tree.getConnectedComponent(999)), "invalid getConnectedComponent must throw")
    require_raises(lambda: weighted.getAltitude(-1), "invalid weighted getAltitude must throw")
    require_raises(lambda: weighted.altitude_of(-1), "invalid weighted altitude_of must throw")
    require_raises(lambda: weighted.getNodeResidue(999), "invalid weighted getNodeResidue must throw")
    require_raises(lambda: tree.mergeNodeIntoParent(-1), "mergeNodeIntoParent must reject invalid NodeId")
    require_raises(lambda: tree.mergeNodeIntoParent(tree.getRoot()), "mergeNodeIntoParent must reject root")
    require_raises(lambda: tree.pruneNode(-1), "pruneNode must reject invalid NodeId")
    require_raises(lambda: tree.pruneNode(tree.getRoot()), "pruneNode must reject root")
    require_raises(lambda: weighted.mergeNodeIntoParent(weighted.getRoot()), "weighted mergeNodeIntoParent must reject root")
    require_raises(lambda: weighted.pruneNode(weighted.getRoot()), "weighted pruneNode must reject root")
    require(weighted.altitude == [0, 1, 2, 3, 4, 5], "weighted altitude property")
    require(weighted.reconstructionImage().shape == (4, 4), "weighted reconstructionImage explicit shape")
    require(
        weighted.reconstructionImage().tolist()
        == [
            [3, 3, 2, 2],
            [3, 4, 4, 2],
            [1, 4, 5, 2],
            [1, 1, 5, 0],
        ],
        "weighted reconstructionImage explicit values",
    )
    weighted.setAltitude(5, 6)
    require(weighted.getAltitude(5) == 6, "weighted setAltitude")
    require(weighted.reconstructionImage().tolist()[2][2] == 6, "weighted reconstruction must reflect external altitude buffer")
    require(weighted.getRepresentativeProperPartsByFlood(3) == [0], "weighted representative proper parts by flood")
    weighted.setAltitudeBuffer([0, 1, 2, 3, 4, 5])
    require(
        weighted.reconstructionImage().tolist()
        == [
            [3, 3, 2, 2],
            [3, 4, 4, 2],
            [1, 4, 5, 2],
            [1, 1, 5, 0],
        ],
        "weighted setAltitudeBuffer round-trip",
    )
    require_raises(lambda: weighted.setAltitudeBuffer([0]), "weighted setAltitudeBuffer must reject wrong size")
    weighted.validateAltitudeBufferShape()
    weighted.validateMonotoneAltitude()
    weighted.altitude = [0, 1, 2, 6, 4, 5]
    require_raises(weighted.validateMonotoneAltitude, "weighted validateMonotoneAltitude must reject broken max-tree order")
    weighted.altitude = [0, 1, 2, 3, 4, 5]

    area_mapping = mmcfilters.Attribute.computeAttributeMapping(tree, mmcfilters.Attribute.AREA)
    require(area_mapping.shape == (4, 4), "attribute mapping shape")
    require(
        area_mapping.tolist()
        == [
            [8.0, 8.0, 12.0, 12.0],
            [8.0, 5.0, 5.0, 12.0],
            [15.0, 5.0, 2.0, 12.0],
            [15.0, 15.0, 2.0, 16.0],
        ],
        "area mapping by pybind",
    )
    weighted_area_mapping = mmcfilters.Attribute.computeAttributeMapping(weighted, mmcfilters.Attribute.AREA)
    require(np.array_equal(weighted_area_mapping, area_mapping), "weighted area mapping by pybind")
    area_attr = mmcfilters.Attribute.computeSingleAttribute(tree, mmcfilters.Attribute.AREA)
    require(area_attr.tolist() == [16.0, 15.0, 12.0, 8.0, 5.0, 2.0], "exact AREA attribute by NodeId")
    weighted_area_attr = mmcfilters.Attribute.computeSingleAttribute(weighted, mmcfilters.Attribute.AREA)
    require(weighted_area_attr.tolist() == area_attr.tolist(), "weighted exact AREA attribute by NodeId")
    default_delta_names, default_delta_attrs = mmcfilters.Attribute.computeSingleAttributeWithDelta(
        weighted, mmcfilters.Attribute.AREA, 1
    )
    require(int(default_delta_attrs[0, default_delta_names["AREA_ASC_1"]]) == 16, "default delta must use last-padding on missing asc")
    require(int(default_delta_attrs[5, default_delta_names["AREA_DESC_1"]]) == 2, "default delta must use last-padding on missing desc")
    delta_names, delta_attrs = mmcfilters.Attribute.computeSingleAttributeWithDelta(
        weighted, mmcfilters.Attribute.AREA, 1, "null-padding"
    )
    require(delta_names["AREA_ASC_1"] == 0, "delta names asc offset")
    require(delta_names["AREA"] == 1, "delta names center offset")
    require(delta_names["AREA_DESC_1"] == 2, "delta names desc offset")
    require(delta_attrs.shape == (tree.numNodes, 3), "delta attribute shape")
    require(np.isnan(delta_attrs[0, delta_names["AREA_ASC_1"]]), "null-padding missing asc must stay NaN")
    require(int(delta_attrs[0, delta_names["AREA"]]) == 16, "null-padding must preserve root center value")
    require(int(delta_attrs[3, delta_names["AREA_ASC_1"]]) == 12, "delta asc attribute by pybind")
    require(int(delta_attrs[3, delta_names["AREA"]]) == 8, "delta center attribute by pybind")
    require(int(delta_attrs[3, delta_names["AREA_DESC_1"]]) == 5, "delta desc attribute by pybind")
    require(int(delta_attrs[5, delta_names["AREA_ASC_1"]]) == 5, "null-padding must preserve available asc value")
    require(int(delta_attrs[5, delta_names["AREA"]]) == 2, "null-padding must preserve leaf center value")
    require(np.isnan(delta_attrs[5, delta_names["AREA_DESC_1"]]), "null-padding missing desc must stay NaN")
    weighted_delta_names, weighted_delta_attrs = mmcfilters.Attribute.computeSingleAttributeWithDelta(
        weighted, mmcfilters.Attribute.AREA, 1, "null-padding"
    )
    require(weighted_delta_names == delta_names, "weighted delta names")
    require(np.array_equal(weighted_delta_attrs, delta_attrs, equal_nan=True), "weighted delta attribute values")
    names, attrs = mmcfilters.Attribute.computeAttributes(
        weighted,
        [mmcfilters.Attribute.AREA, mmcfilters.Attribute.VOLUME, mmcfilters.Attribute.RELATIVE_VOLUME],
    )
    require(names["AREA"] == 0 and names["VOLUME"] == 1 and names["RELATIVE_VOLUME"] == 2, "attribute names map for combined attribute computation")
    require(attrs.shape == (tree.numNodes, 3), "combined attribute array shape")
    require(int(attrs[0, names["AREA"]]) == 16, "combined AREA attribute root")
    require(int(attrs[3, names["VOLUME"]]) == 31, "combined VOLUME attribute node 3")
    require(int(attrs[3, names["RELATIVE_VOLUME"]]) == 22, "combined RELATIVE_VOLUME attribute node 3")
    sparse = mmcfilters.MorphologicalTree.createComponentTree(image, True)
    sparse_weighted_for_delta = mmcfilters.WeightedMorphologicalTree.createComponentTree(image, True)
    sparse.mergeNodeIntoParent(4)
    sparse_weighted_for_delta.mergeNodeIntoParent(4)
    require(sparse.numNodes == 5, "sparse tree live node count")
    require(sparse.numInternalNodeSlots == 6, "sparse tree slot count")
    require_raises(lambda: sparse.getNodeParent(4), "dead-slot getNodeParent must throw")
    require_raises(lambda: sparse.getChildren(4), "dead-slot getChildren must throw")
    require_raises(lambda: sparse_weighted_for_delta.getAltitude(4), "dead-slot weighted getAltitude must throw")
    require_raises(lambda: sparse_weighted_for_delta.getNodeResidue(4), "dead-slot weighted getNodeResidue must throw")
    require_raises(lambda: sparse_weighted_for_delta.setAltitude(4, 7), "dead-slot weighted setAltitude must throw")
    require_raises(lambda: sparse.mergeNodeIntoParent(4), "dead-slot mergeNodeIntoParent must throw")
    require_raises(lambda: sparse.pruneNode(4), "dead-slot pruneNode must throw")
    sparse_area = mmcfilters.Attribute.computeSingleAttribute(sparse, mmcfilters.Attribute.AREA)
    require(sparse_area.shape == (sparse.numInternalNodeSlots,), "single attribute shape must follow internal slots")
    require(int(sparse_area[5]) == 2, "single attribute must preserve sparse slot values")
    sparse_delta_names, sparse_delta = mmcfilters.Attribute.computeSingleAttributeWithDelta(
        sparse_weighted_for_delta, mmcfilters.Attribute.AREA, 1, "null-padding"
    )
    require(sparse_delta.shape == (sparse.numInternalNodeSlots, 3), "delta attribute shape must follow internal slots")
    require(int(sparse_delta[3, sparse_delta_names["AREA"]]) == 8, "delta attribute must preserve live sparse slot values")
    require(np.isnan(sparse_delta[4]).all(), "delta attribute must keep dead sparse slot addressable")
    require(int(sparse_delta[5, sparse_delta_names["AREA_ASC_1"]]) == 8, "delta attribute must preserve sparse leaf asc value")
    require(int(sparse_delta[5, sparse_delta_names["AREA"]]) == 2, "delta attribute must preserve sparse leaf center value")
    require(np.isnan(sparse_delta[5, sparse_delta_names["AREA_DESC_1"]]), "delta attribute must preserve sparse leaf missing desc as NaN")
    sparse_weighted = mmcfilters.WeightedMorphologicalTree.createComponentTree(image, True)
    sparse_weighted.mergeNodeIntoParent(4)
    sparse_names, sparse_attrs = mmcfilters.Attribute.computeAttributes(
        sparse_weighted, [mmcfilters.Attribute.AREA, mmcfilters.Attribute.VOLUME]
    )
    require(sparse_attrs.shape == (sparse.numInternalNodeSlots, 2), "combined attribute shape must follow internal slots")
    require(int(sparse_attrs[5, sparse_names["AREA"]]) == 2, "combined attribute must preserve sparse slot values")

    exported_parent, exported_altitude = weighted.exportHigraHierarchy()
    exported_roundtrip = mmcfilters.WeightedMorphologicalTree.createFromHigraParent(
        exported_parent,
        exported_altitude,
        weighted.numRows,
        weighted.numCols,
        mmcfilters.WeightedMorphologicalTree.MAX_TREE,
        1.5,
    )
    exported_area = weighted.project_node_values_to_exported_higra(weighted_area_attr, mmcfilters.Attribute.AREA)
    roundtrip_exported_area = mmcfilters.Attribute.computeSingleAttribute(
        exported_roundtrip,
        mmcfilters.Attribute.AREA,
        mmcfilters.NodeIdSpace.HIGRA,
    )
    roundtrip_exported_area[: weighted.numTotalProperParts] = 1.0
    require(exported_area.shape == (len(exported_parent),), "exported-Higra projection shape")
    require(np.array_equal(exported_area, roundtrip_exported_area, equal_nan=True), "exported-Higra projection must match import projection")
    weighted_max_dist_attr = mmcfilters.Attribute.computeSingleAttribute(weighted, mmcfilters.Attribute.MAX_DIST)
    exported_area_and_max_dist = weighted.projectNodeValuesToExportedHigra(
        np.stack([weighted_area_attr, weighted_max_dist_attr], axis=1),
        [mmcfilters.Attribute.AREA, mmcfilters.Attribute.MAX_DIST],
    )
    require(exported_area_and_max_dist.shape == (len(exported_parent), 2), "2D exported-Higra projection shape")
    require(np.array_equal(exported_area_and_max_dist[:, 0], exported_area, equal_nan=True), "2D exported-Higra first column")
    require(float(exported_area_and_max_dist[0, 1]) == 0.0, "2D exported-Higra unit MAX_DIST value")
    require_raises(
        lambda: weighted.project_node_values_to_exported_higra(np.array([1.0], dtype=np.float32), mmcfilters.Attribute.AREA),
        "exported-Higra projection must reject wrong node-value size",
    )
    require_raises(
        lambda: weighted.project_node_values_to_exported_higra(weighted_area_attr, [mmcfilters.Attribute.AREA, mmcfilters.Attribute.MAX_DIST]),
        "exported-Higra projection must reject wrong attribute count",
    )

    higra_parent, higra_altitude = build_higra_hierarchy(weighted)

    require_raises(
        lambda: mmcfilters.WeightedMorphologicalTree.createFromHigraParent(
            higra_parent,
            higra_altitude,
            weighted.numRows,
            weighted.numCols,
            weighted.treeType,
        ),
        "Higra max/min import without explicit adjacency should be rejected",
    )

    from_higra = mmcfilters.WeightedMorphologicalTree.createFromHigraParent(
        higra_parent,
        higra_altitude,
        weighted.numRows,
        weighted.numCols,
        weighted.treeType,
        1.5,
    )
    require(from_higra.hasAdjacencyRelation is True, "Higra import with explicit adjacency must preserve it")
    require(from_higra.numHigraNodes == len(higra_parent), "Higra import must expose total Higra node count")
    require(from_higra.getHigraNodeId(3) == weighted.numTotalProperParts + 3, "slot->Higra mapping")
    require(not hasattr(from_higra, "hasHigraNodeIdMapping"), "Higra mapping predicate must not be public")
    require(not hasattr(from_higra, "getNodeIdFromHigra"), "Higra reverse mapping must not be public")
    higra_area = mmcfilters.Attribute.computeSingleAttribute(
        from_higra,
        mmcfilters.Attribute.AREA,
        mmcfilters.NodeIdSpace.HIGRA,
    )
    require(higra_area.shape == (from_higra.numHigraNodes,), "Higra-projected single attribute shape")
    require(int(higra_area[weighted.numTotalProperParts + 3]) == 8, "Higra-projected area by Higra node id")
    require(np.isnan(higra_area[5]), "Higra-projected leaf ids must remain NaN")
    higra_delta_names, higra_delta = mmcfilters.Attribute.computeSingleAttributeWithDelta(
        from_higra,
        mmcfilters.Attribute.AREA,
        1,
        "null-padding",
        mmcfilters.NodeIdSpace.HIGRA,
    )
    require(higra_delta.shape == (from_higra.numHigraNodes, 3), "Higra-projected delta attribute shape")
    require(int(higra_delta[weighted.numTotalProperParts + 3, higra_delta_names["AREA"]]) == 8, "Higra-projected delta center value")
    require(np.isnan(higra_delta[5]).all(), "Higra-projected delta leaf ids must remain addressable")
    higra_names, higra_attrs = mmcfilters.Attribute.computeAttributes(
        from_higra,
        [mmcfilters.Attribute.AREA, mmcfilters.Attribute.VOLUME],
        mmcfilters.NodeIdSpace.HIGRA,
    )
    require(higra_attrs.shape == (from_higra.numHigraNodes, 2), "Higra-projected combined attribute shape")
    require(int(higra_attrs[weighted.numTotalProperParts + 3, higra_names["AREA"]]) == 8, "Higra-projected combined area")

    from_higra_with_adj = mmcfilters.WeightedMorphologicalTree.createFromHigraParent(
        higra_parent,
        higra_altitude,
        weighted.numRows,
        weighted.numCols,
        weighted.treeType,
        1.5,
    )
    require(from_higra_with_adj.hasAdjacencyRelation is True, "explicit Higra adjacency must be preserved")

    exported_higra_parent, exported_higra_altitude = weighted.exportHigraHierarchy()
    require(len(exported_higra_parent) == weighted.numTotalProperParts + weighted.numNodes, "exported Higra hierarchy size")
    exported_roundtrip = mmcfilters.WeightedMorphologicalTree.createFromHigraParent(
        exported_higra_parent,
        exported_higra_altitude,
        weighted.numRows,
        weighted.numCols,
        weighted.treeType,
        1.5,
    )
    require(exported_roundtrip.numNodes == weighted.numNodes, "Higra export round-trip node count")
    require(exported_roundtrip.reconstructionImage().tolist() == weighted.reconstructionImage().tolist(), "Higra export round-trip reconstruction")
    reexported_higra_parent, reexported_higra_altitude = exported_roundtrip.exportHigraHierarchy()
    require(reexported_higra_parent == exported_higra_parent, "Higra export/import parent round-trip")
    require(reexported_higra_altitude == exported_higra_altitude, "Higra export/import altitude round-trip")

    sparse_higra_parent, sparse_higra_altitude = sparse_weighted.exportHigraHierarchy()
    require(len(sparse_higra_parent) == sparse_weighted.numTotalProperParts + sparse_weighted.numNodes, "sparse Higra export must compact dead slots")
    sparse_roundtrip = mmcfilters.WeightedMorphologicalTree.createFromHigraParent(
        sparse_higra_parent,
        sparse_higra_altitude,
        sparse_weighted.numRows,
        sparse_weighted.numCols,
        sparse_weighted.treeType,
        1.5,
    )
    require(sparse_roundtrip.numInternalNodeSlots == sparse_weighted.numNodes, "sparse Higra round-trip slot count")
    sparse_reexported_parent, sparse_reexported_altitude = sparse_roundtrip.exportHigraHierarchy()
    require(sparse_reexported_parent == sparse_higra_parent, "sparse Higra parent round-trip")
    require(sparse_reexported_altitude == sparse_higra_altitude, "sparse Higra altitude round-trip")

    rebuilt = mmcfilters.MorphologicalTree.createFromHigraParent(
        higra_parent,
        tree.numRows,
        tree.numCols,
        tree.treeType,
        1.5,
    )
    require(rebuilt.getRoot() == 0, "Higra topology import root alias")
    require(rebuilt.getAliveNodeIds() == [0, 1, 2, 3, 4, 5], "Higra topology import alive NodeIds")
    require(rebuilt.getChildren(3) == [4], "Higra topology import children")
    require(rebuilt.getSmallestComponent(10) == 5, "Higra topology import smallest component")
    require(rebuilt.getNodeParent(0) == 0, "Higra topology import root parent must point to itself")
    require(rebuilt.getPathToRootNodes(5) == [5, 4, 3, 2, 1, 0], "Higra topology import path to root")
    require(rebuilt.getPathBetweenNodes(5, 2) == [5, 4, 3, 2], "Higra topology import path between nodes")
    require(int(mmcfilters.Attribute.computeSingleAttribute(rebuilt, mmcfilters.Attribute.AREA)[3]) == 8, "Higra topology import area")

    visits = []
    mmcfilters.Attribute.traversePostOrder(
        tree,
        lambda node_id: visits.append(("pre", node_id)),
        lambda parent_id, child_id: visits.append(("merge", parent_id, child_id)),
        lambda node_id: visits.append(("post", node_id)),
    )
    require(visits[0] == ("pre", 0), "attribute traversal starts at root NodeId")
    require(any(step == ("merge", 3, 4) for step in visits), "attribute traversal merge by NodeId")

    attr = np.arange(tree.numNodes, dtype=np.float32)
    ext_values = mmcfilters.ExtinctionValues(tree, attr).getExtinctionValues()
    require(ext_values[0][0] == 5 and ext_values[0][1] == 0, "extinction values by NodeId")

    require(hasattr(mmcfilters.Attribute, "traversePostOrder"), "post-order traversal callback API should exist")
    require(not hasattr(mmcfilters, "NodeMT"), "NodeMT should be removed from Python API")
    require(hasattr(mmcfilters, "ContourRange"), "ContourRange should be exported")
    require(not hasattr(mmcfilters, "ContourProxy"), "legacy ContourProxy alias should be removed")
    require(hasattr(contours, "contoursByNode"), "contoursByNode should be the contour iteration API")
    require(not hasattr(contours, "contours"), "legacy contours() alias should be removed")
    require(not hasattr(contours, "isFullyMaterialized"), "isFullyMaterialized alias should be removed")
    require(tree.root == tree.getRoot(), "root property should expose a NodeId, not a legacy node handle")
    require(not hasattr(tree, "listNodes"), "legacy listNodes handle API should be removed")
    require(not hasattr(tree, "leaves"), "legacy leaves handle API should be removed")
    require(not hasattr(tree, "getAltitude"), "altitude access should move to WeightedMorphologicalTree")
    require(not hasattr(tree, "getNodeResidue"), "node residue access should move to WeightedMorphologicalTree")
    require(not hasattr(tree, "getRepresentativeProperPartsByFlood"), "representative-by-flood should move to WeightedMorphologicalTree")
    require(not hasattr(tree, "reconstructionImage"), "reconstructionImage should be hidden on topology-only trees")
    require(hasattr(weighted, "reconstructionImage"), "reconstructionImage should live on WeightedMorphologicalTree")
    require(not hasattr(mmcfilters.MorphologicalTree, "createFromHigra"), "topology-only MorphologicalTree should not expose weighted Higra import")
    require(not hasattr(mmcfilters.WeightedMorphologicalTree, "createFromHigra"), "legacy weighted Higra import alias should be removed")
    require(not hasattr(tree, "reconstructAltitude"), "reconstructAltitude should be removed from the tree API")
    require(not hasattr(mmcfilters, "reconstructionImage"), "module-level reconstructionImage should be removed")
    require(not hasattr(tree, "getNodeLevel"), "legacy altitude alias should be removed")
    require(not hasattr(tree, "getNodeArea"), "area should be computed via attributes, not stored in the tree")
    require(not hasattr(tree, "getNodeRepresentativeCNP"), "legacy representative proper-part alias should be removed")
    require(not hasattr(tree, "getNodeNumCNPs"), "legacy proper-part count alias should be removed")
    require(not hasattr(tree, "getCNPs"), "legacy proper-part traversal alias should be removed")
    require(not hasattr(tree, "getRepresentativeCNPs"), "legacy representative proper-part alias should be removed")
    require(not hasattr(tree, "getNodeRepresentativeProperPart"), "derived representative proper-part helper should be removed")
    require(not hasattr(tree, "getRepresentativeProperParts"), "derived representative proper-part helper should be removed")
    require(not hasattr(tree, "getPixelsOfCC"), "derived connected-component helper should be removed")
    require(not hasattr(tree, "getGlobalNodeIdFromLegacy"), "legacy slot->global conversion should be hidden from Python API")
    require(not hasattr(tree, "getLegacyNodeIdFromGlobal"), "legacy global->slot conversion should be hidden from Python API")
    require(not hasattr(tree, "rootNodeId"), "legacy rootNodeId property should be removed")
    require(not hasattr(tree, "nodeIdSpaceSize"), "legacy nodeIdSpaceSize property should be removed")
    require(not hasattr(tree, "allocateNode"), "low-level allocateNode mutator should be hidden from Python API")
    require(not hasattr(tree, "releaseNode"), "low-level releaseNode mutator should be hidden from Python API")
    require(not hasattr(tree, "attachNode"), "low-level attachNode mutator should be hidden from Python API")
    require(not hasattr(tree, "detachNode"), "low-level detachNode mutator should be hidden from Python API")
    require(not hasattr(tree, "moveNode"), "low-level moveNode mutator should be hidden from Python API")
    require(not hasattr(tree, "moveChildren"), "low-level moveChildren mutator should be hidden from Python API")
    require(not hasattr(tree, "moveProperPart"), "low-level moveProperPart mutator should be hidden from Python API")
    require(not hasattr(tree, "moveProperParts"), "low-level moveProperParts mutator should be hidden from Python API")
    require(not hasattr(tree, "setRootNode"), "low-level setRoot mutator should be hidden from Python API")
    require(not hasattr(weighted, "tree"), "weighted tree topology must not be exposed as a mutable Python handle")
    require(not hasattr(weighted, "topology"), "weighted topology accessor must stay C++-only and const")
    require(not hasattr(weighted, "edit"), "staged structural edit sessions must not be exposed without Python editor bindings")
    require(not hasattr(weighted, "attachNode"), "weighted low-level attachNode mutator should be hidden from Python API")
    require(not hasattr(weighted, "detachNode"), "weighted low-level detachNode mutator should be hidden from Python API")
    require(not hasattr(weighted, "moveNode"), "weighted low-level moveNode mutator should be hidden from Python API")
    require(not hasattr(weighted, "moveChildren"), "weighted low-level moveChildren mutator should be hidden from Python API")
    require(hasattr(weighted, "pruneNode"), "weighted safe prune mutator should stay public in Python API")
    require(hasattr(weighted, "mergeNodeIntoParent"), "weighted safe merge mutator should stay public in Python API")
    require(hasattr(mmcfilters.ExtinctionValues(tree, attr), "getExtinctionValues"), "extinction tuple API should be exposed under the canonical name")

    tos = mmcfilters.MorphologicalTree.createTreeOfShapes(
        np.array([[1, 2, 1], [2, 3, 2], [1, 2, 1]], dtype=np.uint8),
        mmcfilters.ToSInterpolation.Min4cMax8c,
    )
    weighted_tos = mmcfilters.WeightedMorphologicalTree.createTreeOfShapes(
        np.array([[1, 2, 1], [2, 3, 2], [1, 2, 1]], dtype=np.uint8),
        mmcfilters.ToSInterpolation.Min4cMax8c,
    )
    require(tos.hasAdjacencyRelation is False, "ToS should expose adjacency relation as optional/absent")
    require(tos.hasTreeOfShapesAdjacencyPolicy is True, "ToS should expose auxiliary min/max adjacency policy")
    require(tos.getTreeOfShapesMinTreeAdjacencyRadius() == 1.0, "Min4cMax8c min-tree radius")
    require(tos.getTreeOfShapesMaxTreeAdjacencyRadius() == 1.5, "Min4cMax8c max-tree radius")
    require(tos.getTreeOfShapesMinTreeAdjacencyRelation().radius == 1.0, "Min4cMax8c min-tree adjacency relation radius")
    require(tos.getTreeOfShapesMaxTreeAdjacencyRelation().radius == 1.5, "Min4cMax8c max-tree adjacency relation radius")
    require(tos.getTreeOfShapesMinTreeAdjacencyRelation().size == 5, "Min4cMax8c min-tree adjacency relation size")
    require(tos.getTreeOfShapesMaxTreeAdjacencyRelation().size == 9, "Min4cMax8c max-tree adjacency relation size")
    require(tos.numRows == 3, "ToS must expose image rows")
    require(tos.numCols == 3, "ToS must expose image cols")
    require(weighted_tos.reconstructionImage().shape == (3, 3), "weighted ToS reconstructionImage shape")
    tos_min8_max4 = mmcfilters.MorphologicalTree.createTreeOfShapes(
        np.array([[1, 2], [3, 0]], dtype=np.uint8),
        mmcfilters.ToSInterpolation.Min8cMax4c,
    )
    require(tos_min8_max4.getTreeOfShapesMinTreeAdjacencyRadius() == 1.5, "Min8cMax4c min-tree radius")
    require(tos_min8_max4.getTreeOfShapesMaxTreeAdjacencyRadius() == 1.0, "Min8cMax4c max-tree radius")
    require(tos_min8_max4.getTreeOfShapesMinTreeAdjacencyRelation().size == 9, "Min8cMax4c min-tree adjacency relation size")
    require(tos_min8_max4.getTreeOfShapesMaxTreeAdjacencyRelation().size == 5, "Min8cMax4c max-tree adjacency relation size")
    custom_seed_tos = mmcfilters.MorphologicalTree.createTreeOfShapes(
        np.array([[0, 1]], dtype=np.uint8),
        mmcfilters.ToSInterpolation.SelfDual,
        infinitySeedRow=2,
        infinitySeedCol=4,
    )
    require(custom_seed_tos.getRoot() != -1, "SelfDual ToS must accept a custom infinity seed on the outer boundary")
    internal_seed_tos = mmcfilters.MorphologicalTree.createTreeOfShapes(
        np.array([[0, 1]], dtype=np.uint8),
        mmcfilters.ToSInterpolation.SelfDual,
        infinitySeedRow=1,
        infinitySeedCol=1,
    )
    require(internal_seed_tos.getRoot() != -1, "SelfDual ToS must accept an internal custom infinity seed")
    require_raises(
        lambda: mmcfilters.MorphologicalTree.createTreeOfShapes(
            np.array([[0, 1]], dtype=np.uint8),
            mmcfilters.ToSInterpolation.SelfDual,
            infinitySeedRow=3,
            infinitySeedCol=5,
        ),
        "SelfDual ToS must reject an infinity seed outside the interpolated domain",
    )
    single_tos = mmcfilters.MorphologicalTree.createTreeOfShapes(np.array([[5]], dtype=np.uint8))
    single_weighted_tos = mmcfilters.WeightedMorphologicalTree.createTreeOfShapes(np.array([[5]], dtype=np.uint8))
    require(single_tos.numRows == 1 and single_tos.numCols == 1, "single-pixel default ToS dimensions")
    require(single_tos.getRoot() != -1, "single-pixel default ToS root")
    require(single_weighted_tos.reconstructionImage().tolist() == [[5]], "single-pixel default weighted ToS reconstruction")
    empty = np.empty((0, 0), dtype=np.uint8)
    require_raises(lambda: mmcfilters.MorphologicalTree.createComponentTree(empty, True), "empty component tree must throw")
    require_raises(lambda: mmcfilters.MorphologicalTree.createTreeOfShapes(empty), "empty tree of shapes must throw")
    require_raises(lambda: mmcfilters.WeightedMorphologicalTree.createComponentTree(empty, True), "empty weighted component tree must throw")
    require_raises(lambda: mmcfilters.WeightedMorphologicalTree.createTreeOfShapes(empty), "empty weighted tree of shapes must throw")

    print("python NodeId API ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
