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
        tree_like.getProperPartOwner(pixel_id) + tree_like.numTotalProperParts
        for pixel_id in range(tree_like.numTotalProperParts)
    ]
    higra_parent.extend(
        tree_like.getNodeParent(node_id) + tree_like.numTotalProperParts
        for node_id in tree_like.getAliveNodeIds()
    )

    higra_altitude = [0] * len(higra_parent)
    for pixel_id in range(tree_like.numTotalProperParts):
        higra_altitude[pixel_id] = tree_like.getAltitude(tree_like.getProperPartOwner(pixel_id))
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
    require(hasattr(mmcfilters, "MorphologicalTreeKind"), "package import must expose MorphologicalTreeKind")
    require(hasattr(mmcfilters, "SdrtTiePolicy"), "package import must expose SdrtTiePolicy")
    require(hasattr(mmcfilters, "HierarchySemantics"), "package import must expose HierarchySemantics")
    require(hasattr(mmcfilters, "AltitudeOrder"), "package import must expose AltitudeOrder")
    require(hasattr(mmcfilters, "AdjacencyMode"), "package import must expose AdjacencyMode")
    require(hasattr(mmcfilters, "GridDomain2D"), "package import must expose GridDomain2D")
    max_dist_requirements = mmcfilters.Attribute.requirements(
        mmcfilters.Attribute.MAX_DIST
    )
    require(
        max_dist_requirements
        == {
            "altitude": True,
            "gridDomain2D": True,
            "adjacency": "uniform",
            "monotoneAltitudeOrder": True,
            "altitudeForDirectionalAdjacency": False,
            "canonical4Or8Adjacency": False,
        },
        "Python attribute requirements must expose the C++ capability contract",
    )
    for api_module in (mmcfilters, mmcfilters._native):
        for removed_name in (
            "MorphologicalTree",
            "MorphologicalTreeBase",
            "AdjacencyRelation",
            "DualAdjacencyPolicy",
            "DirectionalAdjacency",
        ):
            require(
                not hasattr(api_module, removed_name),
                f"Python API must not expose removed compatibility name {removed_name}",
            )
        require(
            not hasattr(api_module, "WeightedMorphologicalTreeInt32"),
            "Python API must not expose WeightedMorphologicalTreeInt32 while Python stays uint8-only",
        )
        require(
            not hasattr(api_module, "WeightedMorphologicalTreeFloat32"),
            "Python API must not expose WeightedMorphologicalTreeFloat32 while Python stays uint8-only",
        )
        require(
            not hasattr(api_module, "WeightedTreeView"),
            "Python API must not expose WeightedTreeView while Python stays uint8-only",
        )
    require(not hasattr(mmcfilters.WeightedMorphologicalTree, "MAX_TREE"), "WeightedMorphologicalTree must not expose legacy integer tree-type constants")
    require(not hasattr(mmcfilters.WeightedMorphologicalTree, "MIN_TREE"), "WeightedMorphologicalTree must not expose legacy integer tree-type constants")
    require(not hasattr(mmcfilters.WeightedMorphologicalTree, "TREE_OF_SHAPES"), "WeightedMorphologicalTree must not expose legacy integer tree-type constants")

    image = np.array(
        [
            [3, 3, 2, 2],
            [3, 4, 4, 2],
            [1, 4, 5, 2],
            [1, 1, 5, 0],
        ],
        dtype=np.uint8,
    )

    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.createMaxTree(image.astype(np.int32)),
        "createMaxTree must reject non-uint8 integer arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.createMaxTree(image.astype(np.int64)),
        "createMaxTree must reject int64 arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.createMinTree(image.astype(np.float32)),
        "createMinTree must reject float arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.createMinTree(image.astype(np.float64)),
        "createMinTree must reject float64 arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.createTreeOfShapes(image.astype(np.int64)),
        "createTreeOfShapes must reject int64 arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.createTreeOfShapes(image.astype(np.float64)),
        "createTreeOfShapes must reject float64 arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.createTreeOfShapes(image.astype(bool)),
        "createTreeOfShapes must reject bool arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.createMaxTree(image.astype(object)),
        "createMaxTree must reject object arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.createMaxTree(image[:, ::-1]),
        "createMaxTree must reject non-contiguous uint8 arrays",
    )

    weighted = mmcfilters.MorphologicalTreeFactory.createMaxTree(image)
    # Weighted trees expose the public NodeId topology query API in Python.
    tree = weighted
    weighted_max_tree = mmcfilters.MorphologicalTreeFactory.createMaxTree(image)
    weighted_min_tree = mmcfilters.MorphologicalTreeFactory.createMinTree(image)
    residual_tree = mmcfilters.MorphologicalTreeFactory.createSelfDualResidualTree(
        image,
        radius=1.0,
    )
    saturated_residual_tree = (
        mmcfilters.MorphologicalTreeFactory.createSaturatedSelfDualResidualTree(
            image,
            infinityPixel=0,
            radius=1.0,
        )
    )

    require(weighted_max_tree.getAliveNodeIds() == tree.getAliveNodeIds(), "factory max-tree must expose the expected topology")
    require(weighted_max_tree.descriptiveKind == mmcfilters.MorphologicalTreeKind.MAX_TREE, "factory max-tree tree type")
    require(weighted_min_tree.descriptiveKind == mmcfilters.MorphologicalTreeKind.MIN_TREE, "factory min-tree tree type")
    require(
        residual_tree.descriptiveKind
        == mmcfilters.MorphologicalTreeKind.SELF_DUAL_RESIDUAL_TREE,
        "factory unrestricted residual-tree type",
    )
    require(
        saturated_residual_tree.descriptiveKind
        == mmcfilters.MorphologicalTreeKind.SELF_DUAL_RESIDUAL_TREE,
        "factory saturated residual-tree type",
    )
    require(
        residual_tree.altitudeOrder == mmcfilters.AltitudeOrder.UNCONSTRAINED,
        "residual-tree altitude order capability",
    )
    require(
        residual_tree.reconstructionImage().tolist() == image.tolist(),
        "unrestricted residual-tree exact reconstruction",
    )
    require(
        saturated_residual_tree.reconstructionImage().tolist() == image.tolist(),
        "saturated residual-tree exact reconstruction",
    )
    require(
        weighted_max_tree.altitudeOrder
        == mmcfilters.AltitudeOrder.INCREASING_FROM_ROOT,
        "factory max-tree altitude order capability",
    )
    require(
        weighted_min_tree.altitudeOrder
        == mmcfilters.AltitudeOrder.DECREASING_FROM_ROOT,
        "factory min-tree altitude order capability",
    )
    require(
        weighted_max_tree.adjacencyMode == mmcfilters.AdjacencyMode.UNIFORM,
        "factory max-tree uniform adjacency capability",
    )
    require(
        weighted_max_tree.getUniformGridAdjacency2D().radius == 1.5,
        "factory max-tree uniform adjacency relation",
    )
    rectangular_adjacency = mmcfilters.RegularGridAdjacency2D.rectangular(
        image.shape[0],
        image.shape[1],
        1,
        2,
    )
    require(
        rectangular_adjacency.shape
        == mmcfilters.RegularGridAdjacencyShape.StructuringElement,
        "rectangular adjacency shape",
    )
    require(
        rectangular_adjacency.size == 15,
        "rectangular adjacency stencil size",
    )
    require(
        rectangular_adjacency.neighborIndices(1, 1),
        "rectangular adjacency neighbor traversal",
    )
    custom_adjacency_tree = (
        mmcfilters.MorphologicalTreeFactory.createMaxTree(
            image,
            rectangular_adjacency,
        )
    )
    require(
        custom_adjacency_tree.hasUniformGridAdjacency2D is True,
        "custom max-tree regular-grid adjacency capability",
    )
    require(
        custom_adjacency_tree.getUniformGridAdjacency2D().offsets
        == rectangular_adjacency.offsets,
        "custom max-tree preserves its immutable stencil",
    )
    require_raises(
        lambda: mmcfilters.RegularGridAdjacency2D.fromStructuringElement(
            3,
            3,
            [(0, 0), (0, 1)],
        ),
        "asymmetric structuring element must be rejected",
    )
    require(weighted_max_tree.reconstructionImage().tolist() == weighted.reconstructionImage().tolist(), "weighted createMaxTree reconstruction")
    require(weighted_max_tree.descriptiveKind == mmcfilters.MorphologicalTreeKind.MAX_TREE, "weighted createMaxTree tree type")
    require(weighted_min_tree.descriptiveKind == mmcfilters.MorphologicalTreeKind.MIN_TREE, "weighted createMinTree tree type")
    require(not hasattr(mmcfilters.MorphologicalTreeFactory, "createComponentTree"), "createComponentTree should not be public")
    require(
        not hasattr(mmcfilters.MorphologicalTreeFactory, "createMaxTreeInt32"),
        "createMaxTreeInt32 should not be public while Python stays uint8-only",
    )
    require(
        not hasattr(mmcfilters.MorphologicalTreeFactory, "createMaxTreeFloat32"),
        "createMaxTreeFloat32 should not be public while Python stays uint8-only",
    )
    require(
        not hasattr(mmcfilters.MorphologicalTreeFactory, "createMinTreeInt32"),
        "createMinTreeInt32 should not be public while Python stays uint8-only",
    )
    require(
        not hasattr(mmcfilters.MorphologicalTreeFactory, "createMinTreeFloat32"),
        "createMinTreeFloat32 should not be public while Python stays uint8-only",
    )

    require(not hasattr(mmcfilters.MorphologicalTreeFactory, "createFromTopology"), "createFromTopology should not be public")
    require(not hasattr(mmcfilters.MorphologicalTreeFactory, "create_from_topology"), "create_from_topology should not be public")

    require(tree.getRoot() == 0, "getRoot")
    require(tree.root == tree.getRoot(), "root property")
    require(tree.hasUniformGridAdjacency2D is True, "max-tree should expose adjacency relation context")
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
    require(int(mmcfilters.Attribute.computeSingleAttribute(tree, mmcfilters.Attribute.AREA)[3]) == 8, "node area by attribute computer")
    require(mmcfilters.Attribute.describe(mmcfilters.Attribute.AREA).startswith("Area:"), "attribute description by pybind")
    require(tree.getNodeNumDescendants(2) == 3, "descendants count by NodeId")
    require(tree.getNodeNumSiblings(4) == 0, "siblings count by NodeId")
    require(tree.getNumProperParts(3) == 3, "direct proper-part count by NodeId")
    require(tree.getProperParts(3) == [0, 1, 4], "direct proper parts by NodeId")
    require(list(tree.getConnectedComponent(3)) == [0, 1, 4, 5, 6, 9, 10, 14], "connected component iterator by NodeId")
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
    require(weighted.getNodeResidue(5) == 1, "weighted node residue")
    require_raises(lambda: tree.getNodeParent(-1), "invalid getNodeParent must throw")
    require_raises(lambda: tree.getChildren(-1), "invalid getChildren must throw")
    require_raises(lambda: tree.getNumChildren(999), "invalid getNumChildren must throw")
    require_raises(lambda: tree.getNodeTimePreOrder(999), "invalid getNodeTimePreOrder must throw")
    require_raises(lambda: tree.getProperParts(999), "invalid getProperParts must throw")
    require_raises(lambda: list(tree.getConnectedComponent(999)), "invalid getConnectedComponent must throw")
    require_raises(lambda: weighted.getAltitude(-1), "invalid weighted getAltitude must throw")
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
    require_raises(lambda: weighted.setAltitude(5, 256), "weighted setAltitude must reject altitude above uint8")
    require_raises(lambda: weighted.setAltitude(5, -1), "weighted setAltitude must reject negative altitude")
    require_raises(lambda: weighted.setAltitude(5, np.float32(6.0)), "weighted setAltitude must reject numpy float altitude")
    require_raises(lambda: weighted.setAltitude(5, np.float64(6.0)), "weighted setAltitude must reject numpy float64 altitude")
    require_raises(lambda: weighted.setAltitude(5, True), "weighted setAltitude must reject bool altitude")
    weighted.setAltitudeBuffer([0, 1, 2, 3, 4, 5])
    weighted.setAltitudeBuffer(np.array([0, 1, 2, 3, 4, 5], dtype=np.uint8))
    strict_sample = 5
    strict_parent = weighted.getNodeParent(strict_sample)
    require_raises(
        lambda: weighted.setAltitude(
            strict_sample,
            weighted.getAltitude(strict_parent),
        ),
        "weighted setAltitude must reject equality with the parent",
    )
    equal_altitude = weighted.altitude
    equal_altitude[strict_sample] = equal_altitude[strict_parent]
    require_raises(
        lambda: weighted.setAltitudeBuffer(equal_altitude),
        "weighted setAltitudeBuffer must reject equality with the parent",
    )
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
    require_raises(lambda: weighted.setAltitudeBuffer([0, 1, 2, 3, 4, 256]), "weighted setAltitudeBuffer must reject altitude above uint8")
    require_raises(lambda: weighted.setAltitudeBuffer([0, 1, 2, 3, 4, -1]), "weighted setAltitudeBuffer must reject negative altitude")
    require_raises(
        lambda: weighted.setAltitudeBuffer(np.array([0, 1, 2, 3, 4, 5], dtype=np.int32)),
        "weighted setAltitudeBuffer must reject int32 altitude array",
    )
    require_raises(
        lambda: weighted.setAltitudeBuffer(np.array([0, 1, 2, 3, 4, 5], dtype=np.int64)),
        "weighted setAltitudeBuffer must reject int64 altitude array",
    )
    require_raises(
        lambda: weighted.setAltitudeBuffer(np.array([0, 1, 2, 3, 4, 5], dtype=np.float32)),
        "weighted setAltitudeBuffer must reject float altitude array",
    )
    require_raises(
        lambda: weighted.setAltitudeBuffer(np.array([0, 1, 2, 3, 4, 5], dtype=np.float64)),
        "weighted setAltitudeBuffer must reject float64 altitude array",
    )
    require_raises(
        lambda: weighted.setAltitudeBuffer(np.array([False, True, False, True, False, True], dtype=bool)),
        "weighted setAltitudeBuffer must reject bool altitude array",
    )
    require_raises(
        lambda: weighted.setAltitudeBuffer(np.array([0, 1, 2, 3, 4, 5], dtype=object)),
        "weighted setAltitudeBuffer must reject object altitude array",
    )
    require_raises(
        lambda: weighted.setAltitudeBuffer(np.array([0, 1, 2, 3, 4, 5], dtype=np.uint8).reshape(2, 3)),
        "weighted setAltitudeBuffer must reject 2D uint8 altitude array",
    )
    require_raises(
        lambda: weighted.setAltitudeBuffer(np.array([0, 1, 2, 3, 4, 5], dtype=np.uint8)[::-1]),
        "weighted setAltitudeBuffer must reject non-contiguous uint8 altitude array",
    )
    weighted.validateAltitudeBufferShape()
    weighted.validateMonotoneAltitude()
    require_raises(
        lambda: setattr(weighted, "altitude", [0, 1, 2, 6, 4, 5]),
        "weighted altitude property must reject broken max-tree order",
    )
    require_raises(lambda: setattr(weighted, "altitude", [0, 1, 2, 3, 4, 256]), "weighted altitude property must reject altitude above uint8")
    require_raises(lambda: setattr(weighted, "altitude", [0, 1, 2, 3, 4, -1]), "weighted altitude property must reject negative altitude")
    require_raises(
        lambda: setattr(weighted, "altitude", np.array([0, 1, 2, 3, 4, 5], dtype=np.float32)),
        "weighted altitude property must reject float altitude array",
    )
    require_raises(
        lambda: setattr(weighted, "altitude", np.array([0, 1, 2, 3, 4, 5], dtype=np.int32)),
        "weighted altitude property must reject int32 altitude array",
    )
    require_raises(
        lambda: setattr(weighted, "altitude", np.array([0, 1, 2, 3, 4, 5], dtype=np.int64)),
        "weighted altitude property must reject int64 altitude array",
    )
    require_raises(
        lambda: setattr(weighted, "altitude", np.array([False, True, False, True, False, True], dtype=bool)),
        "weighted altitude property must reject bool altitude array",
    )
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
    require(hasattr(mmcfilters.Attribute, "CONTOUR_PIXELS"), "Python Attribute must expose CONTOUR_PIXELS")
    require(hasattr(mmcfilters.Attribute, "CONTOUR_PERIMETER"), "Python Attribute must expose CONTOUR_PERIMETER")
    contour_names, contour_attrs = mmcfilters.Attribute.computeAttributes(
        weighted,
        [
            mmcfilters.Attribute.CONTOUR_PIXELS,
            mmcfilters.Attribute.CONTOUR_PERIMETER,
            mmcfilters.Attribute.CONTOUR_SIDE_NORTH,
        ],
    )
    require(
        contour_names == {"CONTOUR_PIXELS": 0, "CONTOUR_PERIMETER": 1, "CONTOUR_SIDE_NORTH": 2},
        "contour attribute names map",
    )
    require(contour_attrs.shape == (tree.numNodes, 3), "contour combined attribute shape")
    require(int(contour_attrs[0, contour_names["CONTOUR_PIXELS"]]) == 12, "root CONTOUR_PIXELS")
    require(int(contour_attrs[0, contour_names["CONTOUR_PERIMETER"]]) == 16, "root CONTOUR_PERIMETER")
    require(int(contour_attrs[0, contour_names["CONTOUR_SIDE_NORTH"]]) == 4, "root CONTOUR_SIDE_NORTH")
    contour_perimeter = mmcfilters.Attribute.computeSingleAttribute(weighted, mmcfilters.Attribute.CONTOUR_PERIMETER)
    require(int(contour_perimeter[0]) == 16, "single CONTOUR_PERIMETER by Python binding")
    exported_contour_perimeter = weighted.project_node_values_to_exported_higra(
        contour_perimeter,
        mmcfilters.Attribute.CONTOUR_PERIMETER,
    )
    contour_export_parent, _ = weighted.exportHigraHierarchy()
    contour_export_root = next(index for index, parent_id in enumerate(contour_export_parent) if index == parent_id)
    require(int(exported_contour_perimeter[0]) == 4, "unit CONTOUR_PERIMETER in exported Higra layout")
    require(
        int(exported_contour_perimeter[contour_export_root]) == 16,
        "internal root CONTOUR_PERIMETER in exported Higra layout",
    )
    all_names, all_attrs = mmcfilters.Attribute.computeAttributes(weighted, [mmcfilters.Attribute.ALL])
    require("CONTOUR_PIXELS" in all_names, "Attribute.ALL must include CONTOUR_PIXELS")
    require("CONTOUR_SIDE_SOUTH" in all_names, "Attribute.ALL must include CONTOUR_SIDE_SOUTH")
    require(all_attrs.shape[1] == len(all_names), "Attribute.ALL column count must match names")
    require(int(all_attrs[0, all_names["CONTOUR_PERIMETER"]]) == 16, "Attribute.ALL root CONTOUR_PERIMETER")

    require(not hasattr(mmcfilters, "BitquadDeltas"), "Python package must not expose internal BitquadDeltas")
    boundary_names, boundary_attrs = mmcfilters.Attribute.computeTopologyAttributes(
        weighted,
        [mmcfilters.Attribute.Group.BOUNDARY],
    )
    require("BITQUADS_AREA" in boundary_names, "topology boundary group must include BITQUADS_AREA")
    require("BITQUADS_PERIMETER" in boundary_names, "topology boundary group must include BITQUADS_PERIMETER")
    require("CONTOUR_SIDE_SOUTH" in boundary_names, "topology boundary group must include CONTOUR_SIDE_SOUTH")
    require(boundary_attrs.shape == (tree.numInternalNodeSlots, len(boundary_names)), "topology boundary group shape")

    sparse = mmcfilters.MorphologicalTreeFactory.createMaxTree(image)
    sparse_weighted_for_delta = mmcfilters.MorphologicalTreeFactory.createMaxTree(image)
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
    sparse_weighted = mmcfilters.MorphologicalTreeFactory.createMaxTree(image)
    sparse_weighted.mergeNodeIntoParent(4)
    sparse_names, sparse_attrs = mmcfilters.Attribute.computeAttributes(
        sparse_weighted, [mmcfilters.Attribute.AREA, mmcfilters.Attribute.VOLUME]
    )
    require(sparse_attrs.shape == (sparse.numInternalNodeSlots, 2), "combined attribute shape must follow internal slots")
    require(int(sparse_attrs[5, sparse_names["AREA"]]) == 2, "combined attribute must preserve sparse slot values")

    exported_parent, exported_altitude = weighted.exportHigraHierarchy()
    exported_roundtrip = mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
        exported_parent,
        exported_altitude,
        weighted.numRows,
        weighted.numCols,
        mmcfilters.MorphologicalTreeKind.MAX_TREE,
        1.5,
    )
    exported_area = weighted.project_node_values_to_exported_higra(weighted_area_attr, mmcfilters.Attribute.AREA)
    roundtrip_exported_area = mmcfilters.Attribute.computeSingleAttribute(
        exported_roundtrip,
        mmcfilters.Attribute.AREA,
        mmcfilters.NodeIdSpace.HIGRA,
    )
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
    invalid_high_higra_altitude = list(higra_altitude)
    invalid_high_higra_altitude[0] = 256
    invalid_low_higra_altitude = list(higra_altitude)
    invalid_low_higra_altitude[0] = -1
    uint8_higra_altitude = np.array(higra_altitude, dtype=np.uint8)

    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
            higra_parent,
            higra_altitude,
            weighted.numRows,
            weighted.numCols,
            mmcfilters.MorphologicalTreeKind.MAX_TREE,
        ),
        "Higra max/min import without explicit adjacency should be rejected",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
            higra_parent,
            invalid_high_higra_altitude,
            weighted.numRows,
            weighted.numCols,
            mmcfilters.MorphologicalTreeKind.MAX_TREE,
            1.5,
        ),
        "Higra import must reject altitude above uint8",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
            higra_parent,
            invalid_low_higra_altitude,
            weighted.numRows,
            weighted.numCols,
            mmcfilters.MorphologicalTreeKind.MAX_TREE,
            1.5,
        ),
        "Higra import must reject negative altitude",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
            higra_parent,
            np.array(higra_altitude, dtype=np.int32),
            weighted.numRows,
            weighted.numCols,
            mmcfilters.MorphologicalTreeKind.MAX_TREE,
            1.5,
        ),
        "Higra import must reject int32 altitude arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
            higra_parent,
            np.array(higra_altitude, dtype=np.int64),
            weighted.numRows,
            weighted.numCols,
            mmcfilters.MorphologicalTreeKind.MAX_TREE,
            1.5,
        ),
        "Higra import must reject int64 altitude arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
            higra_parent,
            np.array(higra_altitude, dtype=np.float32),
            weighted.numRows,
            weighted.numCols,
            mmcfilters.MorphologicalTreeKind.MAX_TREE,
            1.5,
        ),
        "Higra import must reject float altitude arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
            higra_parent,
            np.array(higra_altitude, dtype=np.float64),
            weighted.numRows,
            weighted.numCols,
            mmcfilters.MorphologicalTreeKind.MAX_TREE,
            1.5,
        ),
        "Higra import must reject float64 altitude arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
            higra_parent,
            np.array(higra_altitude, dtype=bool),
            weighted.numRows,
            weighted.numCols,
            mmcfilters.MorphologicalTreeKind.MAX_TREE,
            1.5,
        ),
        "Higra import must reject bool altitude arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
            higra_parent,
            np.array(higra_altitude, dtype=object),
            weighted.numRows,
            weighted.numCols,
            mmcfilters.MorphologicalTreeKind.MAX_TREE,
            1.5,
        ),
        "Higra import must reject object altitude arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
            higra_parent,
            uint8_higra_altitude.reshape(2, -1),
            weighted.numRows,
            weighted.numCols,
            mmcfilters.MorphologicalTreeKind.MAX_TREE,
            1.5,
        ),
        "Higra import must reject 2D uint8 altitude arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
            higra_parent,
            uint8_higra_altitude[::-1],
            weighted.numRows,
            weighted.numCols,
            mmcfilters.MorphologicalTreeKind.MAX_TREE,
            1.5,
        ),
        "Higra import must reject non-contiguous uint8 altitude arrays",
    )

    from_higra = mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
        higra_parent,
        uint8_higra_altitude,
        weighted.numRows,
        weighted.numCols,
        mmcfilters.MorphologicalTreeKind.MAX_TREE,
        1.5,
    )
    require(from_higra.hasUniformGridAdjacency2D is True, "Higra import with explicit adjacency must preserve it")
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
    require(int(higra_area[5]) == 1, "Higra-projected leaf ids must receive unit AREA values")
    higra_delta_names, higra_delta = mmcfilters.Attribute.computeSingleAttributeWithDelta(
        from_higra,
        mmcfilters.Attribute.AREA,
        1,
        "null-padding",
        mmcfilters.NodeIdSpace.HIGRA,
    )
    require(higra_delta.shape == (from_higra.numHigraNodes, 3), "Higra-projected delta attribute shape")
    require(int(higra_delta[weighted.numTotalProperParts + 3, higra_delta_names["AREA"]]) == 8, "Higra-projected delta center value")
    require(np.array_equal(higra_delta[5], np.ones(3)), "Higra-projected delta leaf ids must receive unit AREA values")
    higra_names, higra_attrs = mmcfilters.Attribute.computeAttributes(
        from_higra,
        [mmcfilters.Attribute.AREA, mmcfilters.Attribute.VOLUME],
        mmcfilters.NodeIdSpace.HIGRA,
    )
    require(higra_attrs.shape == (from_higra.numHigraNodes, 2), "Higra-projected combined attribute shape")
    require(int(higra_attrs[weighted.numTotalProperParts + 3, higra_names["AREA"]]) == 8, "Higra-projected combined area")
    require(int(higra_attrs[5, higra_names["AREA"]]) == 1, "Higra-projected combined leaf AREA")
    require(int(higra_attrs[5, higra_names["VOLUME"]]) == int(uint8_higra_altitude[5]), "Higra-projected combined leaf VOLUME")

    from_higra_with_adj = mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
        higra_parent,
        higra_altitude,
        weighted.numRows,
        weighted.numCols,
        mmcfilters.MorphologicalTreeKind.MAX_TREE,
        1.5,
    )
    require(from_higra_with_adj.hasUniformGridAdjacency2D is True, "explicit Higra adjacency must be preserved")

    exported_higra_parent, exported_higra_altitude = weighted.exportHigraHierarchy()
    require(len(exported_higra_parent) == weighted.numTotalProperParts + weighted.numNodes, "exported Higra hierarchy size")
    exported_roundtrip = mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
        exported_higra_parent,
        exported_higra_altitude,
        weighted.numRows,
        weighted.numCols,
        mmcfilters.MorphologicalTreeKind.MAX_TREE,
        1.5,
    )
    require(exported_roundtrip.numNodes == weighted.numNodes, "Higra export round-trip node count")
    require(exported_roundtrip.reconstructionImage().tolist() == weighted.reconstructionImage().tolist(), "Higra export round-trip reconstruction")
    reexported_higra_parent, reexported_higra_altitude = exported_roundtrip.exportHigraHierarchy()
    require(reexported_higra_parent == exported_higra_parent, "Higra export/import parent round-trip")
    require(reexported_higra_altitude == exported_higra_altitude, "Higra export/import altitude round-trip")

    sparse_higra_parent, sparse_higra_altitude = sparse_weighted.exportHigraHierarchy()
    require(len(sparse_higra_parent) == sparse_weighted.numTotalProperParts + sparse_weighted.numNodes, "sparse Higra export must compact dead slots")
    sparse_roundtrip = mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
        sparse_higra_parent,
        sparse_higra_altitude,
        sparse_weighted.numRows,
        sparse_weighted.numCols,
        mmcfilters.MorphologicalTreeKind.MAX_TREE,
        1.5,
    )
    require(sparse_roundtrip.numInternalNodeSlots == sparse_weighted.numNodes, "sparse Higra round-trip slot count")
    sparse_reexported_parent, sparse_reexported_altitude = sparse_roundtrip.exportHigraHierarchy()
    require(sparse_reexported_parent == sparse_higra_parent, "sparse Higra parent round-trip")
    require(sparse_reexported_altitude == sparse_higra_altitude, "sparse Higra altitude round-trip")

    rebuilt = mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
        higra_parent,
        higra_altitude,
        tree.numRows,
        tree.numCols,
        mmcfilters.MorphologicalTreeKind.MAX_TREE,
        1.5,
    )
    require(rebuilt.getRoot() == 0, "Higra topology import root alias")
    require(rebuilt.getAliveNodeIds() == [0, 1, 2, 3, 4, 5], "Higra topology import alive NodeIds")
    require(rebuilt.getChildren(3) == [4], "Higra topology import children")
    require(rebuilt.getProperPartOwner(10) == 5, "Higra topology import proper-part owner")
    require(rebuilt.getNodeParent(0) == 0, "Higra topology import root parent must point to itself")
    require(rebuilt.getPathToRootNodes(5) == [5, 4, 3, 2, 1, 0], "Higra topology import path to root")
    require(rebuilt.getPathBetweenNodes(5, 2) == [5, 4, 3, 2], "Higra topology import path between nodes")
    require(int(mmcfilters.Attribute.computeSingleAttribute(rebuilt, mmcfilters.Attribute.AREA)[3]) == 8, "Higra topology import area")

    attr = np.arange(tree.numNodes, dtype=np.float32)
    ext_values = mmcfilters.ExtinctionValues(tree, attr).getRegionalExtrema()
    require(ext_values[0][0] == 5 and ext_values[0][1] == 0, "extinction values by NodeId")

    require(not hasattr(mmcfilters.Attribute, "traversePostOrder"), "post-order traversal callback API should not be public")
    require(not hasattr(mmcfilters, "NodeMT"), "NodeMT should be removed from Python API")
    require(hasattr(mmcfilters, "ContourRange"), "ContourRange should be exported")
    require(not hasattr(mmcfilters, "ContourProxy"), "legacy ContourProxy alias should be removed")
    require(hasattr(contours, "contoursByNode"), "contoursByNode should be the contour iteration API")
    require(not hasattr(contours, "contours"), "legacy contours() alias should be removed")
    require(not hasattr(contours, "isFullyMaterialized"), "isFullyMaterialized alias should be removed")
    require(tree.root == tree.getRoot(), "root property should expose a NodeId, not a legacy node handle")
    require(not hasattr(tree, "listNodes"), "legacy listNodes handle API should be removed")
    require(not hasattr(tree, "leaves"), "legacy leaves handle API should be removed")
    require(hasattr(weighted, "reconstructionImage"), "reconstructionImage should live on WeightedMorphologicalTree")
    require(not hasattr(mmcfilters.WeightedMorphologicalTree, "createFromHigra"), "legacy weighted Higra import alias should be removed")
    require(not hasattr(tree, "reconstructAltitude"), "reconstructAltitude should be removed from the tree API")
    require(not hasattr(mmcfilters, "reconstructionImage"), "module-level reconstructionImage should be removed")
    require(not hasattr(tree, "getNodeLevel"), "legacy altitude alias should be removed")
    require(not hasattr(tree, "getNodeArea"), "area should be computed via attributes, not stored in the tree")
    require(not hasattr(tree, "getNodeRepresentativeCNP"), "legacy representative proper-part alias should be removed")
    require(not hasattr(tree, "getNodeNumCNPs"), "legacy proper-part count alias should be removed")
    require(not hasattr(tree, "getCNPs"), "legacy proper-part traversal alias should be removed")
    require(not hasattr(tree, "getSmallestComponent"), "legacy smallest-component alias should be removed")
    require(not hasattr(tree, "smallestComponentOf"), "legacy smallest-component alias should be removed")
    require(not hasattr(tree, "smallest_component_of"), "legacy smallest-component alias should be removed")
    require(not hasattr(tree, "childrenOf"), "Python childrenOf alias should be removed")
    require(not hasattr(tree, "children_of"), "Python children_of alias should be removed")
    require(not hasattr(tree, "properPartsOf"), "Python properPartsOf alias should be removed")
    require(not hasattr(tree, "proper_parts_of"), "Python proper_parts_of alias should be removed")
    require(not hasattr(tree, "connectedComponentOf"), "Python connectedComponentOf alias should be removed")
    require(not hasattr(tree, "connected_component_of"), "Python connected_component_of alias should be removed")
    require(not hasattr(tree, "parentOf"), "Python parentOf alias should be removed")
    require(not hasattr(tree, "parent_of"), "Python parent_of alias should be removed")
    require(not hasattr(tree, "properPartOwnerOf"), "Python properPartOwnerOf alias should be removed")
    require(not hasattr(tree, "proper_part_owner_of"), "Python proper_part_owner_of alias should be removed")
    require(not hasattr(tree, "nodeSubtreeOf"), "Python nodeSubtreeOf alias should be removed")
    require(not hasattr(tree, "node_subtree_of"), "Python node_subtree_of alias should be removed")
    require(not hasattr(tree, "descendantsOf"), "Python descendantsOf alias should be removed")
    require(not hasattr(tree, "descendants_of"), "Python descendants_of alias should be removed")
    require(not hasattr(weighted, "altitudeOf"), "Python altitudeOf alias should be removed")
    require(not hasattr(weighted, "altitude_of"), "Python altitude_of alias should be removed")
    require(not hasattr(weighted, "residueOf"), "Python residueOf alias should be removed")
    require(not hasattr(weighted, "residue_of"), "Python residue_of alias should be removed")
    require(not hasattr(tree, "getRepresentativeCNPs"), "legacy representative proper-part alias should be removed")
    require(not hasattr(tree, "getNodeRepresentativeProperPart"), "derived representative proper-part helper should be removed")
    require(not hasattr(tree, "getRepresentativeProperParts"), "derived representative proper-part helper should be removed")
    require(not hasattr(weighted, "getRepresentativeProperPartsByFlood"), "flood representative proper-part helper should be removed")
    require(not hasattr(tree, "getPixelsOfCC"), "derived connected-component helper should be removed")
    require(not hasattr(tree, "getGlobalNodeIdFromLegacy"), "legacy slot->global conversion should be hidden from Python API")
    require(not hasattr(tree, "getLegacyNodeIdFromGlobal"), "legacy global->slot conversion should be hidden from Python API")
    require(not hasattr(tree, "rootNodeId"), "legacy rootNodeId property should be removed")
    require(not hasattr(tree, "nodeIdSpaceSize"), "legacy nodeIdSpaceSize property should be removed")
    require(not hasattr(mmcfilters, "ComponentTreeCasf"), "ComponentTreeCasf alias should be removed")
    require(not hasattr(mmcfilters, "ComponentTreeCasfAttribute"), "ComponentTreeCasfAttribute alias should be removed")
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
    require(not hasattr(weighted, "setAltitudeUnchecked"), "unchecked altitude setter must not be public")
    require(not hasattr(weighted, "setAltitudeBufferUnchecked"), "unchecked altitude-buffer setter must not be public")
    require(not hasattr(weighted, "attachNode"), "weighted low-level attachNode mutator should be hidden from Python API")
    require(not hasattr(weighted, "detachNode"), "weighted low-level detachNode mutator should be hidden from Python API")
    require(not hasattr(weighted, "moveNode"), "weighted low-level moveNode mutator should be hidden from Python API")
    require(not hasattr(weighted, "moveChildren"), "weighted low-level moveChildren mutator should be hidden from Python API")
    require(hasattr(weighted, "pruneNode"), "weighted safe prune mutator should stay public in Python API")
    require(hasattr(weighted, "mergeNodeIntoParent"), "weighted safe merge mutator should stay public in Python API")
    require(hasattr(mmcfilters.ExtinctionValues(tree, attr), "getRegionalExtrema"), "regional-extrema tuple API should be exposed under the canonical name")

    tos = mmcfilters.MorphologicalTreeFactory.createTreeOfShapes(
        np.array([[1, 2, 1], [2, 3, 2], [1, 2, 1]], dtype=np.uint8),
        mmcfilters.ToSInterpolation.Min4cMax8c,
    )
    weighted_tos = mmcfilters.MorphologicalTreeFactory.createTreeOfShapes(
        np.array([[1, 2, 1], [2, 3, 2], [1, 2, 1]], dtype=np.uint8),
        mmcfilters.ToSInterpolation.Min4cMax8c,
    )
    require(tos.hasUniformGridAdjacency2D is False, "ToS should expose adjacency relation as optional/absent")
    require(tos.hasDirectionalGridAdjacency2D is True, "ToS should expose directional adjacency capability")
    require(tos.adjacencyMode == mmcfilters.AdjacencyMode.DIRECTIONAL, "ToS directional adjacency mode")
    require(tos.altitudeOrder == mmcfilters.AltitudeOrder.UNCONSTRAINED, "ToS unconstrained altitude order")
    require(tos.getDecreasingGridAdjacency2D().radius == 1.0, "Min4cMax8c decreasing adjacency radius")
    require(tos.getIncreasingGridAdjacency2D().radius == 1.5, "Min4cMax8c increasing adjacency radius")
    require(tos.getDecreasingGridAdjacency2D().size == 5, "Min4cMax8c decreasing adjacency size")
    require(tos.getIncreasingGridAdjacency2D().size == 9, "Min4cMax8c increasing adjacency size")
    require(tos.numRows == 3, "ToS must expose image rows")
    require(tos.numCols == 3, "ToS must expose image cols")
    require(weighted_tos.reconstructionImage().shape == (3, 3), "weighted ToS reconstructionImage shape")
    tos_boundary_names, tos_boundary_attrs = mmcfilters.Attribute.computeAttributes(
        weighted_tos, [mmcfilters.Attribute.Group.BOUNDARY]
    )
    tos_topology_boundary_names, tos_topology_boundary_attrs = mmcfilters.Attribute.computeTopologyAttributes(
        weighted_tos, [mmcfilters.Attribute.Group.BOUNDARY]
    )
    require("BITQUADS_AREA" in tos_boundary_names, "weighted ToS BOUNDARY group must expose BitQuads through computeAttributes")
    require("CONTOUR_SIDE_SOUTH" in tos_boundary_names, "weighted ToS BOUNDARY group must expose contour sides")
    require(tos_boundary_names == tos_topology_boundary_names, "weighted ToS topology BOUNDARY names must match the full pipeline")
    require(
        tos_boundary_attrs.shape == (weighted_tos.numInternalNodeSlots, len(tos_boundary_names)),
        "weighted ToS BOUNDARY shape must use internal node-id space",
    )
    require(
        np.allclose(tos_boundary_attrs, tos_topology_boundary_attrs, equal_nan=True),
        "weighted ToS topology BOUNDARY route must match the full pipeline",
    )
    tos_min8_max4 = mmcfilters.MorphologicalTreeFactory.createTreeOfShapes(
        np.array([[1, 2], [3, 0]], dtype=np.uint8),
        mmcfilters.ToSInterpolation.Min8cMax4c,
    )
    require(tos_min8_max4.getDecreasingGridAdjacency2D().radius == 1.5, "Min8cMax4c decreasing radius")
    require(tos_min8_max4.getIncreasingGridAdjacency2D().radius == 1.0, "Min8cMax4c increasing radius")
    require(tos_min8_max4.getDecreasingGridAdjacency2D().size == 9, "Min8cMax4c decreasing adjacency size")
    require(tos_min8_max4.getIncreasingGridAdjacency2D().size == 5, "Min8cMax4c increasing adjacency size")
    unpadded_options = mmcfilters.TreeOfShapesProducerOptions(
        interpolation=mmcfilters.ToSInterpolation.Min4cMax8c,
        padding=mmcfilters.ToSPaddingPolicy.NoPadding,
        infinitySeedRow=0,
        infinitySeedCol=0,
    )
    require(
        unpadded_options.padding == mmcfilters.ToSPaddingPolicy.NoPadding,
        "ToS producer options must expose the padding policy",
    )
    unpadded_input = np.array([[0, 2, 1], [2, 1, 0]], dtype=np.uint8)
    unpadded_tos = mmcfilters.MorphologicalTreeFactory.createTreeOfShapes(
        unpadded_input,
        unpadded_options,
    )
    require(
        unpadded_tos.numRows == 2 and unpadded_tos.numCols == 3,
        "unpadded ToS must publish the original image domain",
    )
    require(
        unpadded_tos.reconstructionImage().tolist() == unpadded_input.tolist(),
        "unpadded ToS reconstruction",
    )
    virtual_root_tos = mmcfilters.MorphologicalTreeFactory.createTreeOfShapes(
        np.array([[1, 1], [0, 0]], dtype=np.uint8),
        mmcfilters.ToSInterpolation.SelfDual,
    )
    virtual_root = virtual_root_tos.getRoot()
    upper_shape = virtual_root_tos.getProperPartOwner(0)
    lower_shape = virtual_root_tos.getProperPartOwner(2)
    require(virtual_root_tos.numNodes == 3, "SelfDual ToS must preserve the virtual root")
    require(virtual_root_tos.getNumProperParts(virtual_root) == 0, "virtual root must have an empty direct proper part")
    require(virtual_root_tos.isStructuralNode(virtual_root), "virtual root must be derived as structural")
    require(virtual_root_tos.getNumChildren(virtual_root) == 2, "virtual root must retain both child shapes")
    require(
        virtual_root_tos.reconstructionImage().tolist() == [[1, 1], [0, 0]],
        "virtual-root ToS reconstruction",
    )
    native_partial_partition = (
        mmcfilters.MorphologicalTreeFactory.createFromNativeTopology(
            [0, 0, 0],
            [1, 1, 2, 2],
            np.array([0, 1, 0], dtype=np.uint8),
            0,
            2,
            2,
            semantics=mmcfilters.HierarchySemantics(
                descriptiveKind=mmcfilters.MorphologicalTreeKind.GENERIC,
                directionalAdjacency=mmcfilters.DirectionalGridAdjacency2D(
                    mmcfilters.RegularGridAdjacency2D(2, 2, 1.0),
                    mmcfilters.RegularGridAdjacency2D(2, 2, 1.5),
                ),
            ),
        )
    )
    require(
        native_partial_partition.getNumProperParts(0) == 0,
        "generic native partial-partition root may have no direct proper part",
    )
    require(
        native_partial_partition.getNumChildren(0) == 2,
        "generic native partial-partition root children",
    )
    require(
        native_partial_partition.hasDirectionalGridAdjacency2D is True,
        "generic native partial-partition dual adjacency",
    )
    require(
        native_partial_partition.getDecreasingGridAdjacency2D().size == 5,
        "generic native partial-partition decreasing adjacency",
    )
    require(
        native_partial_partition.getIncreasingGridAdjacency2D().size == 9,
        "generic native partial-partition increasing adjacency",
    )
    require(
        native_partial_partition.reconstructionImage().tolist()
        == [[1, 1], [0, 0]],
        "generic native partial-partition reconstruction",
    )
    generic_semantics = mmcfilters.HierarchySemantics(
        altitudeOrder=mmcfilters.AltitudeOrder.INCREASING_FROM_ROOT,
    )
    generic_tree = mmcfilters.MorphologicalTreeFactory.createFromNativeTopology(
        [0, 0, 1],
        [2],
        np.array([0, 1, 2], dtype=np.uint8),
        0,
        1,
        1,
        semantics=generic_semantics,
    )
    require(
        generic_tree.descriptiveKind == mmcfilters.MorphologicalTreeKind.GENERIC,
        "generic native tree descriptive kind",
    )
    require(
        generic_tree.altitudeOrder
        == mmcfilters.AltitudeOrder.INCREASING_FROM_ROOT,
        "generic native tree altitude order",
    )
    require(
        generic_tree.adjacencyMode == mmcfilters.AdjacencyMode.NONE,
        "generic native tree adjacency mode",
    )
    require(
        generic_tree.hasGridDomain2D
        and generic_tree.gridDomain2D.rows == 1
        and generic_tree.gridDomain2D.cols == 1,
        "generic native tree explicit 2D domain",
    )
    require(
        generic_tree.isStructuralNode(0)
        and generic_tree.isStructuralNode(1)
        and not generic_tree.isStructuralNode(2),
        "generic native tree structural-node derivation",
    )
    for node_id in generic_tree.getAliveNodeIds():
        require(
            len(list(generic_tree.getConnectedComponent(node_id))) > 0,
            "every committed generic node must have non-empty subtree support",
        )
    abstract_tree = mmcfilters.MorphologicalTreeFactory.createFromNativeTopology(
        [0, 0, 0],
        [1, 2],
        np.array([10, 3, 20], dtype=np.uint8),
        0,
        semantics=mmcfilters.HierarchySemantics(),
    )
    require(
        abstract_tree.hasGridDomain2D is False,
        "abstract proper-part domain must not invent grid metadata",
    )
    require(
        abstract_tree.gridDomain2D is None,
        "abstract proper-part domain optional grid",
    )
    require(
        abstract_tree.descriptiveKind
        == mmcfilters.MorphologicalTreeKind.GENERIC,
        "canonical descriptive-kind property",
    )
    abstract_area = mmcfilters.Attribute.computeSingleAttribute(
        abstract_tree,
        mmcfilters.Attribute.AREA,
    )
    require(
        int(abstract_area[0]) == 2,
        "support attributes must work on an abstract proper-part domain",
    )
    abstract_gray_height = mmcfilters.Attribute.computeSingleAttribute(
        abstract_tree,
        mmcfilters.Attribute.GRAY_HEIGHT,
    )
    require(
        abstract_gray_height.tolist() == [10.0, 0.0, 0.0],
        "unconstrained GRAY_HEIGHT must use both subtree extrema",
    )
    require_raises(
        lambda: abstract_tree.reconstructionImage(),
        "abstract proper-part domain must reject image reconstruction",
    )
    require_raises(
        lambda: mmcfilters.Attribute.computeSingleAttribute(
            abstract_tree,
            mmcfilters.Attribute.BOX_WIDTH,
        ),
        "abstract proper-part domain must reject geometric attributes",
    )
    require_raises(
        lambda: abstract_tree.numRows,
        "abstract proper-part domain must not expose invented rows",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.createFromNativeTopology(
            [0, 0],
            [0],
            np.array([0, 1], dtype=np.uint8),
            0,
            1,
            1,
            semantics=mmcfilters.HierarchySemantics(),
        ),
        "native factory must reject an attached leaf with empty subtree support",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.createFromNativeTopology(
            [0, 0, 1],
            [2],
            np.array([0, 2, 1], dtype=np.uint8),
            0,
            1,
            1,
            semantics=mmcfilters.HierarchySemantics(
                altitudeOrder=mmcfilters.AltitudeOrder.INCREASING_FROM_ROOT,
            ),
        ),
        "native factory altitude must satisfy the declared generic order",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.createFromNativeTopology(
            [0, 0, 1],
            [2],
            np.array([0, 1, 1], dtype=np.uint8),
            0,
            1,
            1,
            semantics=mmcfilters.HierarchySemantics(
                altitudeOrder=mmcfilters.AltitudeOrder.INCREASING_FROM_ROOT,
            ),
        ),
        "native factory altitude must reject equality in a strict order",
    )
    custom_seed_tos = mmcfilters.MorphologicalTreeFactory.createTreeOfShapes(
        np.array([[0, 1]], dtype=np.uint8),
        mmcfilters.ToSInterpolation.SelfDual,
        infinitySeedRow=2,
        infinitySeedCol=4,
    )
    require(custom_seed_tos.getRoot() != -1, "SelfDual ToS must accept a custom infinity seed on the outer boundary")
    internal_seed_tos = mmcfilters.MorphologicalTreeFactory.createTreeOfShapes(
        np.array([[0, 1]], dtype=np.uint8),
        mmcfilters.ToSInterpolation.SelfDual,
        infinitySeedRow=1,
        infinitySeedCol=1,
    )
    require(internal_seed_tos.getRoot() != -1, "SelfDual ToS must accept an internal custom infinity seed")
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.createTreeOfShapes(
            np.array([[0, 1]], dtype=np.uint8),
            mmcfilters.ToSInterpolation.SelfDual,
            infinitySeedRow=3,
            infinitySeedCol=5,
        ),
        "SelfDual ToS must reject an infinity seed outside the interpolated domain",
    )
    single_tos = mmcfilters.MorphologicalTreeFactory.createTreeOfShapes(np.array([[5]], dtype=np.uint8))
    single_weighted_tos = mmcfilters.MorphologicalTreeFactory.createTreeOfShapes(np.array([[5]], dtype=np.uint8))
    require(single_tos.numRows == 1 and single_tos.numCols == 1, "single-pixel default ToS dimensions")
    require(single_tos.getRoot() != -1, "single-pixel default ToS root")
    require(single_weighted_tos.reconstructionImage().tolist() == [[5]], "single-pixel default weighted ToS reconstruction")

    empty = np.empty((0, 0), dtype=np.uint8)
    require_raises(lambda: mmcfilters.MorphologicalTreeFactory.createMaxTree(empty), "empty max-tree must throw")
    require_raises(lambda: mmcfilters.MorphologicalTreeFactory.createTreeOfShapes(empty), "empty tree of shapes must throw")
    require_raises(lambda: mmcfilters.MorphologicalTreeFactory.createMinTree(empty), "empty min-tree must throw")
    require_raises(lambda: mmcfilters.MorphologicalTreeFactory.createTreeOfShapes(empty), "empty weighted tree of shapes must throw")

    print("python NodeId API ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
