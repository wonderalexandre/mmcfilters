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
        tree_like.smallest_node(pixel_id) + tree_like.num_pixels
        for pixel_id in range(tree_like.num_pixels)
    ]
    higra_parent.extend(
        tree_like.parent(node_id) + tree_like.num_pixels
        for node_id in tree_like.alive_node_ids
    )

    higra_altitude = [0] * len(higra_parent)
    for pixel_id in range(tree_like.num_pixels):
        higra_altitude[pixel_id] = tree_like.node_altitude(tree_like.smallest_node(pixel_id))
    for node_id in tree_like.alive_node_ids:
        higra_altitude[tree_like.num_pixels + node_id] = tree_like.node_altitude(node_id)

    return higra_parent, higra_altitude


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: test_python_nodeid_api.py <build-dir>")

    build_dir = pathlib.Path(sys.argv[1]).resolve()
    mmcfilters = load_native_module(build_dir)

    def complementary_convention(rows, columns, min_radius, max_radius, domain_extension=None, infinity_pixel=0):
        extension = domain_extension if domain_extension is not None else mmcfilters.TopographicDomainExtension.NONE
        return mmcfilters.TopographicConvention(
            mmcfilters.ComplementaryGridImmersion(
                mmcfilters.ComplementaryAdjacencies(
                    mmcfilters.RegularGridAdjacency2D(rows, columns, min_radius),
                    mmcfilters.RegularGridAdjacency2D(rows, columns, max_radius),
                )
            ),
            extension,
            infinity_pixel,
            mmcfilters.TopographicAltitudeEncoding.EXACT_DOUBLED,
        )

    def self_dual_convention(domain_extension=None, infinity_pixel=0):
        extension = domain_extension or mmcfilters.TopographicDomainExtension.EXTERIOR_RING
        return mmcfilters.TopographicConvention(
            mmcfilters.SelfDualSpanImmersion(),
            extension,
            infinity_pixel,
            mmcfilters.TopographicAltitudeEncoding.EXACT_DOUBLED,
        )
    require(hasattr(mmcfilters, "__version__"), "package import must expose __version__")
    require(hasattr(mmcfilters, "ValuedMorphologicalTree"), "package import must expose ValuedMorphologicalTree")
    require(hasattr(mmcfilters, "MorphologicalTreeKind"), "package import must expose MorphologicalTreeKind")
    require(hasattr(mmcfilters, "SpatialOrder"), "package import must expose SpatialOrder")
    require(hasattr(mmcfilters, "RowMajorSpatialOrder"), "package import must expose RowMajorSpatialOrder")
    require(hasattr(mmcfilters, "SelfDualResidualKey"), "package import must expose SelfDualResidualKey")
    require(hasattr(mmcfilters, "SelfDualResidualSchedule"), "package import must expose SelfDualResidualSchedule")
    schedule = mmcfilters.SelfDualResidualSchedule(mmcfilters.SpatialOrder([3, 2, 1, 0]))
    keys = [
        mmcfilters.SelfDualResidualKey(2, 1),
        mmcfilters.SelfDualResidualKey(1, 2),
        mmcfilters.SelfDualResidualKey(1, 3),
    ]
    require(schedule.select_residual_candidate(keys) == 2, "Python residual schedule must use the configured spatial order")
    require(hasattr(mmcfilters, "MorphologicalTreeSemantics"), "package import must expose MorphologicalTreeSemantics")
    require(hasattr(mmcfilters, "NodeAltitudeOrder"), "package import must expose NodeAltitudeOrder")
    require(hasattr(mmcfilters, "SharedAdjacencyContext"), "package import must expose SharedAdjacencyContext")
    require(hasattr(mmcfilters, "TopographicConvention"), "package import must expose TopographicConvention")
    require(hasattr(mmcfilters, "GridDomain2D"), "package import must expose GridDomain2D")
    max_dist_requirements = mmcfilters.Attribute.requirements(
        mmcfilters.Attribute.MAX_DIST
    )
    require(
        max_dist_requirements
        == {
            "altitude": False,
            "grid_domain_2d": True,
            "adjacency": "none",
            "monotone_altitude_order": False,
            "altitude_for_directional_adjacency": False,
            "canonical_4_or_8_adjacency": False,
        },
        "Python attribute requirements must expose the C++ capability contract",
    )
    max_dist_exact_requirements = mmcfilters.Attribute.requirements(
        mmcfilters.Attribute.MAX_DIST_EXACT
    )
    require(
        max_dist_exact_requirements == max_dist_requirements,
        "Python MAX_DIST and MAX_DIST_EXACT must expose the same topology capability contract",
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
            not hasattr(api_module, "ValuedMorphologicalTreeInt32"),
            "Python API must not expose ValuedMorphologicalTreeInt32 while Python stays uint8-only",
        )
        require(
            not hasattr(api_module, "ValuedMorphologicalTreeFloat32"),
            "Python API must not expose ValuedMorphologicalTreeFloat32 while Python stays uint8-only",
        )
        require(
            not hasattr(api_module, "ValuedMorphologicalTreeView"),
            "Python API must not expose ValuedMorphologicalTreeView while Python stays uint8-only",
        )
    require(not hasattr(mmcfilters.ValuedMorphologicalTree, "MAX_TREE"), "ValuedMorphologicalTree must not expose legacy integer tree-type constants")
    require(not hasattr(mmcfilters.ValuedMorphologicalTree, "MIN_TREE"), "ValuedMorphologicalTree must not expose legacy integer tree-type constants")
    require(not hasattr(mmcfilters.ValuedMorphologicalTree, "TREE_OF_SHAPES"), "ValuedMorphologicalTree must not expose legacy integer tree-type constants")

    image = np.array(
        [
            [3, 3, 2, 2],
            [3, 4, 4, 2],
            [1, 4, 5, 2],
            [1, 1, 5, 0],
        ],
        dtype=np.uint8,
    )

    for canonical_factory in (
        "create_max_tree",
        "create_min_tree",
        "create_tree_of_shapes",
        "create_from_native_topology",
        "create_from_higra_parent",
    ):
        require(hasattr(mmcfilters.MorphologicalTreeFactory, canonical_factory), f"missing canonical Python factory {canonical_factory}")
    for legacy_factory in (
        "createMaxTree",
        "createMinTree",
        "createTreeOfShapes",
        "createFromNativeTopology",
        "createFromHigraParent",
    ):
        require(not hasattr(mmcfilters.MorphologicalTreeFactory, legacy_factory), f"legacy Python factory {legacy_factory} must be absent")

    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_max_tree(image.astype(np.int32)),
        "createMaxTree must reject non-uint8 integer arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_max_tree(image.astype(np.int64)),
        "createMaxTree must reject int64 arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_min_tree(image.astype(np.float32)),
        "createMinTree must reject float arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_min_tree(image.astype(np.float64)),
        "createMinTree must reject float64 arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(image.astype(np.int64)),
        "createTreeOfShapes must reject int64 arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(image.astype(np.float64)),
        "createTreeOfShapes must reject float64 arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(image.astype(bool)),
        "createTreeOfShapes must reject bool arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_max_tree(image.astype(object)),
        "createMaxTree must reject object arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_max_tree(image[:, ::-1]),
        "createMaxTree must reject non-contiguous uint8 arrays",
    )

    valued_tree = mmcfilters.MorphologicalTreeFactory.create_max_tree(image)
    # Valued trees expose the public NodeId topology query API in Python.
    tree = valued_tree
    valued_max_tree = mmcfilters.MorphologicalTreeFactory.create_max_tree(image)
    valued_min_tree = mmcfilters.MorphologicalTreeFactory.create_min_tree(image)
    residual_tree = mmcfilters.MorphologicalTreeFactory.create_unrestricted_residual_tree(
        image,
        radius=1.0,
    )
    custom_order_residual_tree = mmcfilters.MorphologicalTreeFactory.create_unrestricted_residual_tree(
        image,
        radius=1.0,
        spatial_order=mmcfilters.SpatialOrder(list(reversed(range(image.size)))),
    )
    saturated_residual_tree = (
        mmcfilters.MorphologicalTreeFactory.create_saturated_residual_tree(
            image,
            infinity_pixel=0,
            radius=1.0,
        )
    )

    require(valued_max_tree.alive_node_ids == tree.alive_node_ids, "factory max-tree must expose the expected topology")
    require(valued_max_tree.kind == mmcfilters.MorphologicalTreeKind.MAX_TREE, "factory max-tree tree type")
    require(valued_min_tree.kind == mmcfilters.MorphologicalTreeKind.MIN_TREE, "factory min-tree tree type")
    require(
        residual_tree.kind
        == mmcfilters.MorphologicalTreeKind.UNRESTRICTED_RESIDUAL_TREE,
        "factory unrestricted residual-tree type",
    )
    require(
        np.array_equal(custom_order_residual_tree.reconstruct_from_node_altitudes(), image),
        "custom spatial order must preserve exact residual reconstruction",
    )
    require(
        saturated_residual_tree.kind
        == mmcfilters.MorphologicalTreeKind.SATURATED_RESIDUAL_TREE,
        "factory saturated residual-tree type",
    )
    require(
        residual_tree.node_altitude_order == mmcfilters.NodeAltitudeOrder.UNCONSTRAINED,
        "residual-tree altitude order capability",
    )
    require(
        residual_tree.reconstruct_from_node_altitudes().tolist() == image.tolist(),
        "unrestricted residual-tree exact reconstruction",
    )
    require(
        saturated_residual_tree.reconstruct_from_node_altitudes().tolist() == image.tolist(),
        "saturated residual-tree exact reconstruction",
    )
    require(
        valued_max_tree.node_altitude_order
        == mmcfilters.NodeAltitudeOrder.INCREASING,
        "factory max-tree altitude order capability",
    )
    require(
        valued_min_tree.node_altitude_order
        == mmcfilters.NodeAltitudeOrder.DECREASING,
        "factory min-tree altitude order capability",
    )
    require(
        isinstance(valued_max_tree.shared_adjacency_context, mmcfilters.SharedAdjacencyContext),
        "factory max-tree shared-adjacency context",
    )
    require(
        valued_max_tree.shared_adjacency_context.adjacency.radius == 1.5,
        "factory max-tree shared adjacency relation",
    )
    rectangular_adjacency = mmcfilters.RegularGridAdjacency2D.rectangular(
        image.shape[0],
        image.shape[1],
        1,
        2,
    )
    require(
        rectangular_adjacency.shape
        == mmcfilters.RegularGridAdjacencyShape.STRUCTURING_ELEMENT,
        "rectangular adjacency shape",
    )
    require(
        rectangular_adjacency.size == 15,
        "rectangular adjacency stencil size",
    )
    require(
        rectangular_adjacency.neighbor_indices(1, 1),
        "rectangular adjacency neighbor traversal",
    )
    custom_adjacency_tree = (
        mmcfilters.MorphologicalTreeFactory.create_max_tree(
            image,
            rectangular_adjacency,
        )
    )
    require(
        isinstance(custom_adjacency_tree.shared_adjacency_context, mmcfilters.SharedAdjacencyContext),
        "custom max-tree shared-adjacency context",
    )
    require(
        custom_adjacency_tree.shared_adjacency_context.adjacency.offsets
        == rectangular_adjacency.offsets,
        "custom max-tree preserves its immutable stencil",
    )
    require_raises(
        lambda: mmcfilters.RegularGridAdjacency2D.from_structuring_element(
            3,
            3,
            [(0, 0), (0, 1)],
        ),
        "asymmetric structuring element must be rejected",
    )
    require(valued_max_tree.reconstruct_from_node_altitudes().tolist() == valued_tree.reconstruct_from_node_altitudes().tolist(), "valued_tree createMaxTree reconstruction")
    require(valued_max_tree.kind == mmcfilters.MorphologicalTreeKind.MAX_TREE, "valued_tree createMaxTree tree type")
    require(valued_min_tree.kind == mmcfilters.MorphologicalTreeKind.MIN_TREE, "valued_tree createMinTree tree type")
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

    require(tree.root == 0, "root")
    require(not callable(tree.root), "root is a property")
    require(isinstance(tree.shared_adjacency_context, mmcfilters.SharedAdjacencyContext), "max-tree should expose shared-adjacency context")
    require(tree.parent(0) == 0, "root parent must point to itself")
    require(tree.alive_node_ids == [0, 1, 2, 3, 4, 5], "alive NodeIds")
    require(tree.num_nodes == len(tree.alive_node_ids), "alive node cardinality")
    contours = mmcfilters.ContourComputation.extraction(tree)
    require(contours.is_materialized is False, "contours must start without global materialization in Python")
    leaf_id = tree.leaves[0]
    leaf_contour_before_root_materialization = list(contours.get_contour(leaf_id))
    require(contours.is_contour_materialized(leaf_id) is True, "get_contour iteration must cache the requested leaf in Python")
    require(contours.is_materialized is False, "get_contour iteration must not materialize all contours in Python")
    list(contours.get_contour(tree.root))
    require(contours.is_materialized is True, "root get_contour iteration must materialize all contours in Python")
    contours.materialize_all()
    require(contours.is_materialized is True, "materialize_all must materialize all contours in Python")
    require(sorted(leaf_contour_before_root_materialization) == sorted(list(contours.get_contour(leaf_id))), "Python incremental contour read must match materialized contour")
    require(tree.leaves == [5], "leaf NodeIds")
    require(tree.num_leaf_nodes == len(tree.leaves), "leaf cardinality")
    require(tree.children(3) == [4], "children by NodeId")
    require(int(mmcfilters.Attribute.compute_single_attribute(tree, mmcfilters.Attribute.AREA)[3]) == 8, "node area by attribute computer")
    require(mmcfilters.Attribute.describe(mmcfilters.Attribute.AREA).startswith("Area:"), "attribute description by pybind")
    require(tree.num_descendants(2) == 3, "descendants count by NodeId")
    require(tree.num_siblings(4) == 0, "siblings count by NodeId")
    require(tree.ancestors(5) == [5, 4, 3, 2, 1, 0], "ancestors include self and end at root")
    require(tree.ancestors(tree.root) == [tree.root], "root ancestor chain contains only self")
    require(tree.descendants(2) == [3, 4, 5], "descendants are strict")
    require(tree.subtree_nodes(2) == [2, 3, 4, 5], "subtree nodes include the queried node")
    require(tree.post_order() == [5, 4, 3, 2, 1, 0], "post-order traversal schedule")
    require(tree.breadth_first_traversal() == [0, 1, 2, 3, 4, 5], "breadth-first traversal")
    require(tree.lowest_common_ancestor(5, 2) == 2, "lowest common ancestor")
    require(sorted(tree.alive_node_ids, key=tree.dfs_entry_index) == tree.subtree_nodes(tree.root), "DFS entry order equals pre-order")
    require(sorted(tree.alive_node_ids, key=tree.dfs_exit_index) == tree.post_order(), "DFS exit order equals post-order")
    dfs_event_indices = [index for node in tree.alive_node_ids for index in (tree.dfs_entry_index(node), tree.dfs_exit_index(node))]
    require(sorted(dfs_event_indices) == list(range(2 * tree.num_nodes)), "DFS entry and exit indices form one event sequence")
    require(tree.dfs_entry_index(tree.root) == 0, "root DFS entry index")
    require(tree.dfs_exit_index(tree.root) == 2 * tree.num_nodes - 1, "root DFS exit index")
    for node in tree.alive_node_ids:
        entry_index = tree.dfs_entry_index(node)
        exit_index = tree.dfs_exit_index(node)
        require(entry_index < exit_index, "DFS entry precedes exit")
        require(tree.num_descendants(node) == (exit_index - entry_index - 1) // 2, "descendant count from DFS interval")
    require(tree.proper_part_cardinality(3) == 3, "direct proper-part count by NodeId")
    require(tree.proper_part(3) == [0, 1, 4], "direct proper parts by NodeId")
    require(list(tree.node_support(3)) == [0, 1, 4, 5, 6, 9, 10, 14], "connected component iterator by NodeId")
    require(tree.is_tree_of_partial_partitions(), "component tree must be a tree of partial partitions")
    tree.validate_tree_of_partial_partitions()
    require(len(tree.smallest_node_map) == tree.num_pixels, "smallest-node map size")
    for pixel in range(tree.num_pixels):
        require(pixel in tree.proper_part(tree.smallest_node(pixel)), "pixel must belong to its smallest node's proper part")
    require(tree.reconstruct_node(3).shape == (4, 4), "node reconstruction shape by NodeId")
    require(
        tree.reconstruct_node(3).tolist()
        == [
            [255, 255, 0, 0],
            [255, 255, 255, 0],
            [0, 255, 255, 0],
            [0, 0, 255, 0],
        ],
        "node reconstruction values by NodeId",
    )
    require(int(tree.reconstruct_node(3).sum()) == 255 * 8, "node reconstruction by NodeId")
    require(valued_tree.root == tree.root, "valued_tree root")
    require(not callable(valued_tree.root), "valued_tree root is a property")
    require(valued_tree.alive_node_ids == tree.alive_node_ids, "valued_tree alive NodeIds")
    require(valued_tree.num_nodes == len(valued_tree.alive_node_ids), "valued_tree alive-node cardinality")
    require(valued_tree.node_altitude(0) == 0, "valued_tree root altitude")
    require(valued_tree.node_altitude(5) == 5, "valued_tree nodeAltitude")
    require(valued_tree.node_residue(5) == 1, "valued_tree node residue")
    require_raises(lambda: tree.parent(-1), "invalid parent must throw")
    require_raises(lambda: tree.children(-1), "invalid children must throw")
    require_raises(lambda: tree.num_children(999), "invalid numChildren must throw")
    require_raises(lambda: tree.dfs_entry_index(999), "invalid dfsEntryIndex must throw")
    require_raises(lambda: tree.proper_part(999), "invalid properPart must throw")
    require_raises(lambda: list(tree.node_support(999)), "invalid node_support must throw")
    require_raises(lambda: valued_tree.node_altitude(-1), "invalid valued_tree nodeAltitude must throw")
    require_raises(lambda: valued_tree.node_residue(999), "invalid valued_tree nodeResidue must throw")
    require_raises(lambda: tree.merge_node_into_parent(-1), "mergeNodeIntoParent must reject invalid NodeId")
    require_raises(lambda: tree.merge_node_into_parent(tree.root), "mergeNodeIntoParent must reject root")
    require_raises(lambda: tree.prune_node(-1), "pruneNode must reject invalid NodeId")
    require_raises(lambda: tree.prune_node(tree.root), "pruneNode must reject root")
    require_raises(lambda: valued_tree.merge_node_into_parent(valued_tree.root), "valued_tree mergeNodeIntoParent must reject root")
    require_raises(lambda: valued_tree.prune_node(valued_tree.root), "valued_tree pruneNode must reject root")
    require(np.array_equal(valued_tree.node_altitudes, np.array([0, 1, 2, 3, 4, 5], dtype=np.uint8)), "valued_tree altitude property")
    require(valued_tree.reconstruct_from_node_altitudes().shape == (4, 4), "valued_tree reconstructFromNodeAltitudes explicit shape")
    require(
        valued_tree.reconstruct_from_node_altitudes().tolist()
        == [
            [3, 3, 2, 2],
            [3, 4, 4, 2],
            [1, 4, 5, 2],
            [1, 1, 5, 0],
        ],
        "valued_tree reconstructFromNodeAltitudes explicit values",
    )
    valued_tree.set_node_altitude(5, 6)
    require(valued_tree.node_altitude(5) == 6, "valued_tree setNodeAltitude")
    require(valued_tree.reconstruct_from_node_altitudes().tolist()[2][2] == 6, "valued_tree reconstruction must reflect external altitude buffer")
    require_raises(lambda: valued_tree.set_node_altitude(5, 256), "valued_tree setNodeAltitude must reject altitude above uint8")
    require_raises(lambda: valued_tree.set_node_altitude(5, -1), "valued_tree setNodeAltitude must reject negative altitude")
    require_raises(lambda: valued_tree.set_node_altitude(5, np.float32(6.0)), "valued_tree setNodeAltitude must reject numpy float altitude")
    require_raises(lambda: valued_tree.set_node_altitude(5, np.float64(6.0)), "valued_tree setNodeAltitude must reject numpy float64 altitude")
    require_raises(lambda: valued_tree.set_node_altitude(5, True), "valued_tree setNodeAltitude must reject bool altitude")
    valued_tree.node_altitudes = [0, 1, 2, 3, 4, 5]
    valued_tree.node_altitudes = np.array([0, 1, 2, 3, 4, 5], dtype=np.uint8)
    strict_sample = 5
    strict_parent = valued_tree.parent(strict_sample)
    require_raises(
        lambda: valued_tree.set_node_altitude(
            strict_sample,
            valued_tree.node_altitude(strict_parent),
        ),
        "valued_tree setNodeAltitude must reject equality with the parent",
    )
    equal_altitude = valued_tree.node_altitudes
    equal_altitude[strict_sample] = equal_altitude[strict_parent]
    require_raises(
        lambda: setattr(valued_tree, "node_altitudes", equal_altitude),
        "valued_tree setNodeAltitudes must reject equality with the parent",
    )
    require(
        valued_tree.reconstruct_from_node_altitudes().tolist()
        == [
            [3, 3, 2, 2],
            [3, 4, 4, 2],
            [1, 4, 5, 2],
            [1, 1, 5, 0],
        ],
        "valued_tree setNodeAltitudes round-trip",
    )
    require_raises(lambda: setattr(valued_tree, "node_altitudes", [0]), "node_altitudes must reject wrong size")
    require_raises(lambda: setattr(valued_tree, "node_altitudes", [0, 1, 2, 3, 4, 256]), "node_altitudes must reject altitude above uint8")
    require_raises(lambda: setattr(valued_tree, "node_altitudes", [0, 1, 2, 3, 4, -1]), "node_altitudes must reject negative altitude")
    require_raises(
        lambda: setattr(valued_tree, "node_altitudes", np.array([0, 1, 2, 3, 4, 5], dtype=np.int32)),
        "node_altitudes must reject int32 altitude array",
    )
    require_raises(
        lambda: setattr(valued_tree, "node_altitudes", np.array([0, 1, 2, 3, 4, 5], dtype=np.int64)),
        "node_altitudes must reject int64 altitude array",
    )
    require_raises(
        lambda: setattr(valued_tree, "node_altitudes", np.array([0, 1, 2, 3, 4, 5], dtype=np.float32)),
        "node_altitudes must reject float altitude array",
    )
    require_raises(
        lambda: setattr(valued_tree, "node_altitudes", np.array([0, 1, 2, 3, 4, 5], dtype=np.float64)),
        "node_altitudes must reject float64 altitude array",
    )
    require_raises(
        lambda: setattr(valued_tree, "node_altitudes", np.array([False, True, False, True, False, True], dtype=bool)),
        "node_altitudes must reject bool altitude array",
    )
    require_raises(
        lambda: setattr(valued_tree, "node_altitudes", np.array([0, 1, 2, 3, 4, 5], dtype=object)),
        "node_altitudes must reject object altitude array",
    )
    require_raises(
        lambda: setattr(valued_tree, "node_altitudes", np.array([0, 1, 2, 3, 4, 5], dtype=np.uint8).reshape(2, 3)),
        "node_altitudes must reject 2D uint8 altitude array",
    )
    require_raises(
        lambda: setattr(valued_tree, "node_altitudes", np.array([0, 1, 2, 3, 4, 5], dtype=np.uint8)[::-1]),
        "node_altitudes must reject non-contiguous uint8 altitude array",
    )
    valued_tree.validate_node_altitude_buffer_shape()
    valued_tree.validate_monotone_node_altitudes()
    require_raises(
        lambda: setattr(valued_tree, "node_altitudes", [0, 1, 2, 6, 4, 5]),
        "valued_tree altitude property must reject broken max-tree order",
    )
    require_raises(lambda: setattr(valued_tree, "node_altitudes", [0, 1, 2, 3, 4, 256]), "node_altitudes must reject altitude above uint8")
    require_raises(lambda: setattr(valued_tree, "node_altitudes", [0, 1, 2, 3, 4, -1]), "node_altitudes must reject negative altitude")
    require_raises(
        lambda: setattr(valued_tree, "node_altitudes", np.array([0, 1, 2, 3, 4, 5], dtype=np.float32)),
        "valued_tree altitude property must reject float altitude array",
    )
    require_raises(
        lambda: setattr(valued_tree, "node_altitudes", np.array([0, 1, 2, 3, 4, 5], dtype=np.int32)),
        "valued_tree altitude property must reject int32 altitude array",
    )
    require_raises(
        lambda: setattr(valued_tree, "node_altitudes", np.array([0, 1, 2, 3, 4, 5], dtype=np.int64)),
        "valued_tree altitude property must reject int64 altitude array",
    )
    require_raises(
        lambda: setattr(valued_tree, "node_altitudes", np.array([False, True, False, True, False, True], dtype=bool)),
        "valued_tree altitude property must reject bool altitude array",
    )
    valued_tree.node_altitudes = [0, 1, 2, 3, 4, 5]

    area_mapping = mmcfilters.Attribute.compute_attribute_mapping(tree, mmcfilters.Attribute.AREA)
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
    valued_tree_area_mapping = mmcfilters.Attribute.compute_attribute_mapping(valued_tree, mmcfilters.Attribute.AREA)
    require(np.array_equal(valued_tree_area_mapping, area_mapping), "valued_tree area mapping by pybind")
    area_attr = mmcfilters.Attribute.compute_single_attribute(tree, mmcfilters.Attribute.AREA)
    require(area_attr.tolist() == [16.0, 15.0, 12.0, 8.0, 5.0, 2.0], "exact AREA attribute by NodeId")
    valued_tree_area_attribute = mmcfilters.Attribute.compute_single_attribute(valued_tree, mmcfilters.Attribute.AREA)
    require(valued_tree_area_attribute.tolist() == area_attr.tolist(), "valued_tree exact AREA attribute by NodeId")
    default_sample_layout, default_sample_values = mmcfilters.Attribute.compute_sampled_node_attribute(
        valued_tree, mmcfilters.Attribute.AREA, 1, 1
    )
    require(
        int(default_sample_values[0, default_sample_layout["AREA_ANCESTOR_1"]]) == 16,
        "default sampling must repeat the nearest value for a missing ancestor",
    )
    require(
        int(default_sample_values[5, default_sample_layout["AREA_DESCENDANT_1"]]) == 2,
        "default sampling must repeat the nearest value for a missing descendant",
    )
    sample_layout, sampled_values = mmcfilters.Attribute.compute_sampled_node_attribute(
        valued_tree,
        mmcfilters.Attribute.AREA,
        1,
        1,
        missing_sample_policy=mmcfilters.MissingNodeAttributeSamplePolicy.NOT_A_NUMBER,
    )
    require(sample_layout["AREA_ANCESTOR_1"] == 0, "ancestor sample column")
    require(sample_layout["AREA"] == 1, "current-node sample column")
    require(sample_layout["AREA_DESCENDANT_1"] == 2, "descendant sample column")
    require(sampled_values.shape == (tree.num_nodes, 3), "sampled node-attribute shape")
    require(np.isnan(sampled_values[0, sample_layout["AREA_ANCESTOR_1"]]), "missing ancestor must stay NaN")
    require(int(sampled_values[0, sample_layout["AREA"]]) == 16, "sampling must preserve root value")
    require(int(sampled_values[3, sample_layout["AREA_ANCESTOR_1"]]) == 12, "ancestor sample by pybind")
    require(int(sampled_values[3, sample_layout["AREA"]]) == 8, "current-node sample by pybind")
    require(int(sampled_values[3, sample_layout["AREA_DESCENDANT_1"]]) == 5, "descendant sample by pybind")
    require(int(sampled_values[5, sample_layout["AREA_ANCESTOR_1"]]) == 5, "available ancestor sample")
    require(int(sampled_values[5, sample_layout["AREA"]]) == 2, "leaf current-node sample")
    require(np.isnan(sampled_values[5, sample_layout["AREA_DESCENDANT_1"]]), "missing descendant must stay NaN")
    repeated_layout, repeated_values = mmcfilters.Attribute.compute_sampled_node_attribute(
        valued_tree,
        mmcfilters.Attribute.AREA,
        1,
        1,
        missing_sample_policy=mmcfilters.MissingNodeAttributeSamplePolicy.NOT_A_NUMBER,
    )
    require(repeated_layout == sample_layout, "repeated sampled layout")
    require(np.array_equal(repeated_values, sampled_values, equal_nan=True), "repeated sampled values")

    missing_policy_results = {}
    for policy in (
        mmcfilters.MissingNodeAttributeSamplePolicy.REPEAT_NEAREST,
        mmcfilters.MissingNodeAttributeSamplePolicy.NOT_A_NUMBER,
        mmcfilters.MissingNodeAttributeSamplePolicy.ZERO,
    ):
        policy_layout, policy_values = mmcfilters.Attribute.compute_sampled_node_attribute(
            valued_tree,
            mmcfilters.Attribute.AREA,
            100,
            1,
            sampling_policy=mmcfilters.NodeAttributeSamplingPolicy.LARGEST_SUPPORT_DESCENDANT,
            missing_sample_policy=policy,
        )
        missing_policy_results[policy] = (policy_layout, policy_values)

    repeat_layout, repeat_values = missing_policy_results[mmcfilters.MissingNodeAttributeSamplePolicy.REPEAT_NEAREST]
    require(int(repeat_values[3, repeat_layout["AREA_ANCESTOR_1"]]) == 8, "RepeatNearest ancestor result")
    require(int(repeat_values[3, repeat_layout["AREA_DESCENDANT_1"]]) == 8, "RepeatNearest descendant result")
    nan_layout, nan_values = missing_policy_results[mmcfilters.MissingNodeAttributeSamplePolicy.NOT_A_NUMBER]
    require(np.isnan(nan_values[3, nan_layout["AREA_ANCESTOR_1"]]), "NotANumber ancestor result")
    require(np.isnan(nan_values[3, nan_layout["AREA_DESCENDANT_1"]]), "NotANumber descendant result")
    zero_layout, zero_values = missing_policy_results[mmcfilters.MissingNodeAttributeSamplePolicy.ZERO]
    require(int(zero_values[3, zero_layout["AREA_ANCESTOR_1"]]) == 0, "Zero ancestor result")
    require(int(zero_values[3, zero_layout["AREA_DESCENDANT_1"]]) == 0, "Zero descendant result")
    require_raises(
        lambda: mmcfilters.Attribute.compute_sampled_node_attribute(valued_tree, mmcfilters.Attribute.AREA, 0, 1),
        "altitude_step must be positive",
    )
    names, attrs = mmcfilters.Attribute.compute_attributes(
        valued_tree,
        [mmcfilters.Attribute.AREA, mmcfilters.Attribute.VOLUME, mmcfilters.Attribute.RELATIVE_VOLUME],
    )
    require(names["AREA"] == 0 and names["VOLUME"] == 1 and names["RELATIVE_VOLUME"] == 2, "attribute names map for combined attribute computation")
    require(attrs.shape == (tree.num_nodes, 3), "combined attribute array shape")
    require(int(attrs[0, names["AREA"]]) == 16, "combined AREA attribute root")
    require(int(attrs[3, names["VOLUME"]]) == 31, "combined VOLUME attribute node 3")
    require(int(attrs[3, names["RELATIVE_VOLUME"]]) == 22, "combined RELATIVE_VOLUME attribute node 3")
    require(hasattr(mmcfilters.Attribute, "CONTOUR_PIXELS"), "Python Attribute must expose CONTOUR_PIXELS")
    require(hasattr(mmcfilters.Attribute, "CONTOUR_PERIMETER"), "Python Attribute must expose CONTOUR_PERIMETER")
    contour_names, contour_attrs = mmcfilters.Attribute.compute_attributes(
        valued_tree,
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
    require(contour_attrs.shape == (tree.num_nodes, 3), "contour combined attribute shape")
    require(int(contour_attrs[0, contour_names["CONTOUR_PIXELS"]]) == 12, "root CONTOUR_PIXELS")
    require(int(contour_attrs[0, contour_names["CONTOUR_PERIMETER"]]) == 16, "root CONTOUR_PERIMETER")
    require(int(contour_attrs[0, contour_names["CONTOUR_SIDE_NORTH"]]) == 4, "root CONTOUR_SIDE_NORTH")
    contour_perimeter = mmcfilters.Attribute.compute_single_attribute(valued_tree, mmcfilters.Attribute.CONTOUR_PERIMETER)
    require(int(contour_perimeter[0]) == 16, "single CONTOUR_PERIMETER by Python binding")
    exported_contour_perimeter = valued_tree.project_node_values_to_exported_higra(
        contour_perimeter,
        mmcfilters.Attribute.CONTOUR_PERIMETER,
    )
    contour_export_parent, _ = valued_tree.export_higra_hierarchy()
    contour_export_root = next(index for index, parent_id in enumerate(contour_export_parent) if index == parent_id)
    require(int(exported_contour_perimeter[0]) == 4, "unit CONTOUR_PERIMETER in exported Higra layout")
    require(
        int(exported_contour_perimeter[contour_export_root]) == 16,
        "internal root CONTOUR_PERIMETER in exported Higra layout",
    )
    all_names, all_attrs = mmcfilters.Attribute.compute_attributes(valued_tree, [mmcfilters.Attribute.ALL])
    require("CONTOUR_PIXELS" in all_names, "Attribute.ALL must include CONTOUR_PIXELS")
    require("CONTOUR_SIDE_SOUTH" in all_names, "Attribute.ALL must include CONTOUR_SIDE_SOUTH")
    require(all_attrs.shape[1] == len(all_names), "Attribute.ALL column count must match names")
    require(int(all_attrs[0, all_names["CONTOUR_PERIMETER"]]) == 16, "Attribute.ALL root CONTOUR_PERIMETER")

    require(not hasattr(mmcfilters, "BitquadDeltas"), "Python package must not expose internal BitquadDeltas")
    boundary_names, boundary_attrs = mmcfilters.Attribute.compute_topology_attributes(
        valued_tree,
        [mmcfilters.Attribute.Group.BOUNDARY],
    )
    require("BITQUAD_AREA" in boundary_names, "topology boundary group must include BITQUAD_AREA")
    require("BITQUAD_PERIMETER" in boundary_names, "topology boundary group must include BITQUAD_PERIMETER")
    require("CONTOUR_SIDE_SOUTH" in boundary_names, "topology boundary group must include CONTOUR_SIDE_SOUTH")
    require(boundary_attrs.shape == (tree.num_internal_node_slots, len(boundary_names)), "topology boundary group shape")

    sparse = mmcfilters.MorphologicalTreeFactory.create_max_tree(image)
    sparse_valued_tree_for_sampling = mmcfilters.MorphologicalTreeFactory.create_max_tree(image)
    sparse.merge_node_into_parent(4)
    sparse_valued_tree_for_sampling.merge_node_into_parent(4)
    require(sparse.num_nodes == 5, "sparse tree live node count")
    require(sparse.num_internal_node_slots == 6, "sparse tree slot count")
    require_raises(lambda: sparse.parent(4), "dead-slot parent must throw")
    require_raises(lambda: sparse.children(4), "dead-slot children must throw")
    require_raises(lambda: sparse_valued_tree_for_sampling.node_altitude(4), "dead-slot valued_tree nodeAltitude must throw")
    require_raises(lambda: sparse_valued_tree_for_sampling.node_residue(4), "dead-slot valued_tree nodeResidue must throw")
    require_raises(lambda: sparse_valued_tree_for_sampling.set_node_altitude(4, 7), "dead-slot valued_tree setNodeAltitude must throw")
    require_raises(lambda: sparse.merge_node_into_parent(4), "dead-slot mergeNodeIntoParent must throw")
    require_raises(lambda: sparse.prune_node(4), "dead-slot pruneNode must throw")
    sparse_area = mmcfilters.Attribute.compute_single_attribute(sparse, mmcfilters.Attribute.AREA)
    require(sparse_area.shape == (sparse.num_internal_node_slots,), "single attribute shape must follow internal slots")
    require(int(sparse_area[5]) == 2, "single attribute must preserve sparse slot values")
    sparse_sample_layout, sparse_samples = mmcfilters.Attribute.compute_sampled_node_attribute(
        sparse_valued_tree_for_sampling,
        mmcfilters.Attribute.AREA,
        1,
        1,
        missing_sample_policy=mmcfilters.MissingNodeAttributeSamplePolicy.NOT_A_NUMBER,
    )
    require(sparse_samples.shape == (sparse.num_internal_node_slots, 3), "sampled values must follow internal slots")
    require(int(sparse_samples[3, sparse_sample_layout["AREA"]]) == 8, "sampling must preserve live sparse slot values")
    require(np.isnan(sparse_samples[4]).all(), "sampling must keep a dead sparse slot addressable")
    require(int(sparse_samples[5, sparse_sample_layout["AREA_ANCESTOR_1"]]) == 8, "sparse leaf ancestor sample")
    require(int(sparse_samples[5, sparse_sample_layout["AREA"]]) == 2, "sparse leaf current-node sample")
    require(np.isnan(sparse_samples[5, sparse_sample_layout["AREA_DESCENDANT_1"]]), "sparse leaf missing descendant")
    sparse_valued_tree = mmcfilters.MorphologicalTreeFactory.create_max_tree(image)
    sparse_valued_tree.merge_node_into_parent(4)
    sparse_names, sparse_attrs = mmcfilters.Attribute.compute_attributes(
        sparse_valued_tree, [mmcfilters.Attribute.AREA, mmcfilters.Attribute.VOLUME]
    )
    require(sparse_attrs.shape == (sparse.num_internal_node_slots, 2), "combined attribute shape must follow internal slots")
    require(int(sparse_attrs[5, sparse_names["AREA"]]) == 2, "combined attribute must preserve sparse slot values")

    exported_parent, exported_altitude = valued_tree.export_higra_hierarchy()
    exported_roundtrip = mmcfilters.MorphologicalTreeFactory.create_from_higra_parent(
        parent=exported_parent,
        node_altitudes=exported_altitude,
        rows=valued_tree.num_rows,
        columns=valued_tree.num_columns,
        kind=mmcfilters.MorphologicalTreeKind.MAX_TREE,
        radius=1.5,
    )
    exported_area = valued_tree.project_node_values_to_exported_higra(valued_tree_area_attribute, mmcfilters.Attribute.AREA)
    roundtrip_exported_area = mmcfilters.Attribute.compute_single_attribute(
        exported_roundtrip,
        mmcfilters.Attribute.AREA,
        mmcfilters.NodeIdSpace.HIGRA,
    )
    require(exported_area.shape == (len(exported_parent),), "exported-Higra projection shape")
    require(np.array_equal(exported_area, roundtrip_exported_area, equal_nan=True), "exported-Higra projection must match import projection")
    valued_tree_max_dist_attribute = mmcfilters.Attribute.compute_single_attribute(valued_tree, mmcfilters.Attribute.MAX_DIST)
    topology_max_dist_attribute = mmcfilters.Attribute.compute_single_topology_attribute(valued_tree, mmcfilters.Attribute.MAX_DIST)
    require(
        np.array_equal(valued_tree_max_dist_attribute, topology_max_dist_attribute, equal_nan=True),
        "Python MAX_DIST must use the same topology-only computation through valued and topology entry points",
    )
    valued_tree_max_dist_exact = mmcfilters.Attribute.compute_single_attribute(
        valued_tree, mmcfilters.Attribute.MAX_DIST_EXACT
    )
    topology_max_dist_exact = mmcfilters.Attribute.compute_single_topology_attribute(
        valued_tree, mmcfilters.Attribute.MAX_DIST_EXACT
    )
    require(
        np.array_equal(
            valued_tree_max_dist_exact,
            topology_max_dist_exact,
            equal_nan=True,
        ),
        "Python MAX_DIST_EXACT must use the same topology-only computation through valued and topology entry points",
    )
    exported_area_and_max_dist = valued_tree.project_node_values_to_exported_higra(
        np.stack([valued_tree_area_attribute, valued_tree_max_dist_attribute], axis=1),
        [mmcfilters.Attribute.AREA, mmcfilters.Attribute.MAX_DIST],
    )
    require(exported_area_and_max_dist.shape == (len(exported_parent), 2), "2D exported-Higra projection shape")
    require(np.array_equal(exported_area_and_max_dist[:, 0], exported_area, equal_nan=True), "2D exported-Higra first column")
    require(float(exported_area_and_max_dist[0, 1]) == 0.0, "2D exported-Higra unit MAX_DIST value")
    require_raises(
        lambda: valued_tree.project_node_values_to_exported_higra(np.array([1.0], dtype=np.float32), mmcfilters.Attribute.AREA),
        "exported-Higra projection must reject wrong node-value size",
    )
    require_raises(
        lambda: valued_tree.project_node_values_to_exported_higra(valued_tree_area_attribute, [mmcfilters.Attribute.AREA, mmcfilters.Attribute.MAX_DIST]),
        "exported-Higra projection must reject wrong attribute count",
    )

    higra_parent, higra_altitude = build_higra_hierarchy(valued_tree)
    invalid_high_higra_altitude = list(higra_altitude)
    invalid_high_higra_altitude[0] = 256
    invalid_low_higra_altitude = list(higra_altitude)
    invalid_low_higra_altitude[0] = -1
    uint8_higra_altitude = np.array(higra_altitude, dtype=np.uint8)

    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_from_higra_parent(
            higra_parent,
            higra_altitude,
            valued_tree.num_rows,
            valued_tree.num_columns,
            mmcfilters.MorphologicalTreeKind.MAX_TREE,
        ),
        "Higra max/min import without explicit adjacency should be rejected",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_from_higra_parent(
            higra_parent,
            invalid_high_higra_altitude,
            valued_tree.num_rows,
            valued_tree.num_columns,
            mmcfilters.MorphologicalTreeKind.MAX_TREE,
            1.5,
        ),
        "Higra import must reject altitude above uint8",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_from_higra_parent(
            higra_parent,
            invalid_low_higra_altitude,
            valued_tree.num_rows,
            valued_tree.num_columns,
            mmcfilters.MorphologicalTreeKind.MAX_TREE,
            1.5,
        ),
        "Higra import must reject negative altitude",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_from_higra_parent(
            higra_parent,
            np.array(higra_altitude, dtype=np.int32),
            valued_tree.num_rows,
            valued_tree.num_columns,
            mmcfilters.MorphologicalTreeKind.MAX_TREE,
            1.5,
        ),
        "Higra import must reject int32 altitude arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_from_higra_parent(
            higra_parent,
            np.array(higra_altitude, dtype=np.int64),
            valued_tree.num_rows,
            valued_tree.num_columns,
            mmcfilters.MorphologicalTreeKind.MAX_TREE,
            1.5,
        ),
        "Higra import must reject int64 altitude arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_from_higra_parent(
            higra_parent,
            np.array(higra_altitude, dtype=np.float32),
            valued_tree.num_rows,
            valued_tree.num_columns,
            mmcfilters.MorphologicalTreeKind.MAX_TREE,
            1.5,
        ),
        "Higra import must reject float altitude arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_from_higra_parent(
            higra_parent,
            np.array(higra_altitude, dtype=np.float64),
            valued_tree.num_rows,
            valued_tree.num_columns,
            mmcfilters.MorphologicalTreeKind.MAX_TREE,
            1.5,
        ),
        "Higra import must reject float64 altitude arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_from_higra_parent(
            higra_parent,
            np.array(higra_altitude, dtype=bool),
            valued_tree.num_rows,
            valued_tree.num_columns,
            mmcfilters.MorphologicalTreeKind.MAX_TREE,
            1.5,
        ),
        "Higra import must reject bool altitude arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_from_higra_parent(
            higra_parent,
            np.array(higra_altitude, dtype=object),
            valued_tree.num_rows,
            valued_tree.num_columns,
            mmcfilters.MorphologicalTreeKind.MAX_TREE,
            1.5,
        ),
        "Higra import must reject object altitude arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_from_higra_parent(
            higra_parent,
            uint8_higra_altitude.reshape(2, -1),
            valued_tree.num_rows,
            valued_tree.num_columns,
            mmcfilters.MorphologicalTreeKind.MAX_TREE,
            1.5,
        ),
        "Higra import must reject 2D uint8 altitude arrays",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_from_higra_parent(
            higra_parent,
            uint8_higra_altitude[::-1],
            valued_tree.num_rows,
            valued_tree.num_columns,
            mmcfilters.MorphologicalTreeKind.MAX_TREE,
            1.5,
        ),
        "Higra import must reject non-contiguous uint8 altitude arrays",
    )

    from_higra = mmcfilters.MorphologicalTreeFactory.create_from_higra_parent(
        higra_parent,
        uint8_higra_altitude,
        valued_tree.num_rows,
        valued_tree.num_columns,
        mmcfilters.MorphologicalTreeKind.MAX_TREE,
        1.5,
    )
    require(isinstance(from_higra.shared_adjacency_context, mmcfilters.SharedAdjacencyContext), "Higra import must preserve shared adjacency")
    require(from_higra.num_higra_nodes == len(higra_parent), "Higra import must expose total Higra node count")
    require(from_higra.higra_node_id(3) == valued_tree.num_pixels + 3, "slot->Higra mapping")
    require(not hasattr(from_higra, "hasHigraNodeIdMapping"), "Higra mapping predicate must not be public")
    require(not hasattr(from_higra, "getNodeIdFromHigra"), "Higra reverse mapping must not be public")
    higra_area = mmcfilters.Attribute.compute_single_attribute(
        from_higra,
        mmcfilters.Attribute.AREA,
        mmcfilters.NodeIdSpace.HIGRA,
    )
    require(higra_area.shape == (from_higra.num_higra_nodes,), "Higra-projected single attribute shape")
    require(int(higra_area[valued_tree.num_pixels + 3]) == 8, "Higra-projected area by Higra node id")
    require(int(higra_area[5]) == 1, "Higra-projected leaf ids must receive unit AREA values")
    higra_sample_layout, higra_samples = mmcfilters.Attribute.compute_sampled_node_attribute(
        from_higra,
        mmcfilters.Attribute.AREA,
        1,
        1,
        missing_sample_policy=mmcfilters.MissingNodeAttributeSamplePolicy.NOT_A_NUMBER,
        output_space=mmcfilters.NodeIdSpace.HIGRA,
    )
    require(higra_samples.shape == (from_higra.num_higra_nodes, 3), "Higra-projected sampled attribute shape")
    require(int(higra_samples[valued_tree.num_pixels + 3, higra_sample_layout["AREA"]]) == 8, "Higra-projected current-node sample")
    require(np.array_equal(higra_samples[5], np.ones(3)), "Higra-projected leaf ids must receive unit AREA samples")
    higra_names, higra_attrs = mmcfilters.Attribute.compute_attributes(
        from_higra,
        [mmcfilters.Attribute.AREA, mmcfilters.Attribute.VOLUME],
        mmcfilters.NodeIdSpace.HIGRA,
    )
    require(higra_attrs.shape == (from_higra.num_higra_nodes, 2), "Higra-projected combined attribute shape")
    require(int(higra_attrs[valued_tree.num_pixels + 3, higra_names["AREA"]]) == 8, "Higra-projected combined area")
    require(int(higra_attrs[5, higra_names["AREA"]]) == 1, "Higra-projected combined leaf AREA")
    require(int(higra_attrs[5, higra_names["VOLUME"]]) == int(uint8_higra_altitude[5]), "Higra-projected combined leaf VOLUME")

    from_higra_with_adj = mmcfilters.MorphologicalTreeFactory.create_from_higra_parent(
        higra_parent,
        higra_altitude,
        valued_tree.num_rows,
        valued_tree.num_columns,
        mmcfilters.MorphologicalTreeKind.MAX_TREE,
        1.5,
    )
    require(isinstance(from_higra_with_adj.shared_adjacency_context, mmcfilters.SharedAdjacencyContext), "explicit Higra adjacency must be preserved")

    exported_higra_parent, exported_higra_altitude = valued_tree.export_higra_hierarchy()
    require(len(exported_higra_parent) == valued_tree.num_pixels + valued_tree.num_nodes, "exported Higra hierarchy size")
    exported_roundtrip = mmcfilters.MorphologicalTreeFactory.create_from_higra_parent(
        exported_higra_parent,
        exported_higra_altitude,
        valued_tree.num_rows,
        valued_tree.num_columns,
        mmcfilters.MorphologicalTreeKind.MAX_TREE,
        1.5,
    )
    require(exported_roundtrip.num_nodes == valued_tree.num_nodes, "Higra export round-trip node count")
    require(exported_roundtrip.reconstruct_from_node_altitudes().tolist() == valued_tree.reconstruct_from_node_altitudes().tolist(), "Higra export round-trip reconstruction")
    reexported_higra_parent, reexported_higra_altitude = exported_roundtrip.export_higra_hierarchy()
    require(reexported_higra_parent == exported_higra_parent, "Higra export/import parent round-trip")
    require(reexported_higra_altitude == exported_higra_altitude, "Higra export/import altitude round-trip")

    sparse_higra_parent, sparse_higra_altitude = sparse_valued_tree.export_higra_hierarchy()
    require(len(sparse_higra_parent) == sparse_valued_tree.num_pixels + sparse_valued_tree.num_nodes, "sparse Higra export must compact dead slots")
    sparse_roundtrip = mmcfilters.MorphologicalTreeFactory.create_from_higra_parent(
        sparse_higra_parent,
        sparse_higra_altitude,
        sparse_valued_tree.num_rows,
        sparse_valued_tree.num_columns,
        mmcfilters.MorphologicalTreeKind.MAX_TREE,
        1.5,
    )
    require(sparse_roundtrip.num_internal_node_slots == sparse_valued_tree.num_nodes, "sparse Higra round-trip slot count")
    sparse_reexported_parent, sparse_reexported_altitude = sparse_roundtrip.export_higra_hierarchy()
    require(sparse_reexported_parent == sparse_higra_parent, "sparse Higra parent round-trip")
    require(sparse_reexported_altitude == sparse_higra_altitude, "sparse Higra altitude round-trip")

    rebuilt = mmcfilters.MorphologicalTreeFactory.create_from_higra_parent(
        higra_parent,
        higra_altitude,
        tree.num_rows,
        tree.num_columns,
        mmcfilters.MorphologicalTreeKind.MAX_TREE,
        1.5,
    )
    require(rebuilt.root == 0, "Higra topology import root alias")
    require(rebuilt.alive_node_ids == [0, 1, 2, 3, 4, 5], "Higra topology import alive NodeIds")
    require(rebuilt.children(3) == [4], "Higra topology import children")
    require(rebuilt.smallest_node(10) == 5, "Higra topology import smallest node")
    require(rebuilt.parent(0) == 0, "Higra topology import root parent must point to itself")
    require(rebuilt.ancestors(5) == [5, 4, 3, 2, 1, 0], "Higra topology import path to root")
    require(rebuilt.path_between_nodes(5, 2) == [5, 4, 3, 2], "Higra topology import path between nodes")
    require(int(mmcfilters.Attribute.compute_single_attribute(rebuilt, mmcfilters.Attribute.AREA)[3]) == 8, "Higra topology import area")

    attr = np.arange(tree.num_nodes, dtype=np.float32)
    ext_values = mmcfilters.ExtinctionValues(tree, attr).get_regional_extrema()
    require(ext_values[0][0] == 5 and ext_values[0][1] == 0, "extinction values by NodeId")

    require(not hasattr(mmcfilters.Attribute, "traversePostOrder"), "post-order traversal callback API should not be public")
    require(not hasattr(mmcfilters, "NodeMT"), "NodeMT should be removed from Python API")
    require(hasattr(mmcfilters, "ContourRange"), "ContourRange should be exported")
    require(not hasattr(mmcfilters, "ContourProxy"), "legacy ContourProxy alias should be removed")
    require(hasattr(contours, "contours_by_node"), "contours_by_node should be the contour iteration API")
    require(not hasattr(contours, "contours"), "legacy contours() alias should be removed")
    require(not hasattr(contours, "isFullyMaterialized"), "isFullyMaterialized alias should be removed")
    require(tree.root == tree.root, "root property should expose a NodeId, not a legacy node handle")
    require(not hasattr(tree, "listNodes"), "legacy listNodes handle API should be removed")
    require(tree.leaves == [5], "leaves is the canonical structural query")
    require(hasattr(valued_tree, "reconstruct_from_node_altitudes"), "reconstruct_from_node_altitudes should live on ValuedMorphologicalTree")
    require(not hasattr(valued_tree, "reconstructFromNodeAltitudes"), "legacy reconstructFromNodeAltitudes should be removed")
    require(not hasattr(mmcfilters.ValuedMorphologicalTree, "createFromHigra"), "legacy valued_tree Higra import alias should be removed")
    require(not hasattr(tree, "reconstructAltitude"), "reconstructAltitude should be removed from the tree API")
    require(not hasattr(mmcfilters, "reconstructFromNodeAltitudes"), "module-level reconstructFromNodeAltitudes should be removed")
    zero_baseline_tree = mmcfilters.MorphologicalTreeFactory.create_max_tree(np.array([[37]], dtype=np.uint8), 1.5)
    require(zero_baseline_tree.node_residue(zero_baseline_tree.root) == 37,
            "root residue should equal the root altitude under the fixed zero reconstruction baseline")
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
    require(not hasattr(tree, "nodeSubtreeOf"), "Python nodeSubtreeOf alias should be removed")
    require(not hasattr(tree, "node_subtree_of"), "Python node_subtree_of alias should be removed")
    require(not hasattr(tree, "descendantsOf"), "Python descendantsOf alias should be removed")
    require(not hasattr(tree, "descendants_of"), "Python descendants_of alias should be removed")
    require(not hasattr(valued_tree, "altitudeOf"), "Python altitudeOf alias should be removed")
    require(not hasattr(valued_tree, "altitude_of"), "Python altitude_of alias should be removed")
    require(not hasattr(valued_tree, "residueOf"), "Python residueOf alias should be removed")
    require(not hasattr(valued_tree, "residue_of"), "Python residue_of alias should be removed")
    require(not hasattr(tree, "getRepresentativeCNPs"), "legacy representative proper-part alias should be removed")
    require(not hasattr(tree, "getNodeRepresentativeProperPart"), "derived representative proper-part helper should be removed")
    require(not hasattr(tree, "getRepresentativeProperParts"), "derived representative proper-part helper should be removed")
    require(not hasattr(valued_tree, "getRepresentativeProperPartsByFlood"), "flood representative proper-part helper should be removed")
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
    require(not hasattr(tree, "movePixelToProperPart"), "low-level movePixelToProperPart mutator should be hidden from Python API")
    require(not hasattr(tree, "mergeProperParts"), "low-level mergeProperParts mutator should be hidden from Python API")
    require(not hasattr(tree, "setRootNode"), "low-level setRoot mutator should be hidden from Python API")
    require(not hasattr(valued_tree, "tree"), "valued tree topology must not be exposed as a mutable Python handle")
    require(not hasattr(valued_tree, "topology"), "valued_tree topology accessor must stay C++-only and const")
    require(not hasattr(valued_tree, "edit"), "staged structural edit sessions must not be exposed without Python editor bindings")
    require(not hasattr(valued_tree, "setAltitudeUnchecked"), "unchecked altitude setter must not be public")
    require(not hasattr(valued_tree, "setAltitudeBufferUnchecked"), "unchecked altitude-buffer setter must not be public")
    require(not hasattr(valued_tree, "attachNode"), "valued_tree low-level attachNode mutator should be hidden from Python API")
    require(not hasattr(valued_tree, "detachNode"), "valued_tree low-level detachNode mutator should be hidden from Python API")
    require(not hasattr(valued_tree, "moveNode"), "valued_tree low-level moveNode mutator should be hidden from Python API")
    require(not hasattr(valued_tree, "moveChildren"), "valued_tree low-level moveChildren mutator should be hidden from Python API")
    require(hasattr(valued_tree, "prune_node"), "safe prune mutator should stay public in Python API")
    require(hasattr(valued_tree, "merge_node_into_parent"), "safe merge mutator should stay public in Python API")
    require(not hasattr(valued_tree, "pruneNode"), "legacy pruneNode should be removed")
    require(not hasattr(valued_tree, "mergeNodeIntoParent"), "legacy mergeNodeIntoParent should be removed")
    require(hasattr(mmcfilters.ExtinctionValues(tree, attr), "get_regional_extrema"), "regional-extrema tuple API should be exposed under the canonical name")

    tos = mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(
        np.array([[1, 2, 1], [2, 3, 2], [1, 2, 1]], dtype=np.uint8),
        complementary_convention(3, 3, 1.0, 1.5),
    )
    valued_tree_of_shapes = mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(
        np.array([[1, 2, 1], [2, 3, 2], [1, 2, 1]], dtype=np.uint8),
        complementary_convention(3, 3, 1.0, 1.5),
    )
    require(tos.shared_adjacency_context is None, "ToS should not expose a shared-adjacency context")
    require(isinstance(tos.topographic_convention, mmcfilters.TopographicConvention), "ToS should retain a topographic convention")
    require(tos.node_altitude_order == mmcfilters.NodeAltitudeOrder.UNCONSTRAINED, "ToS unconstrained altitude order")
    tos_adjacencies = tos.topographic_convention.immersion.complementary_adjacencies
    require(tos_adjacencies.min_adjacency.radius == 1.0, "Min4cMax8c minimum adjacency radius")
    require(tos_adjacencies.max_adjacency.radius == 1.5, "Min4cMax8c maximum adjacency radius")
    require(tos_adjacencies.min_adjacency.size == 5, "Min4cMax8c minimum adjacency size")
    require(tos_adjacencies.max_adjacency.size == 9, "Min4cMax8c maximum adjacency size")
    require(tos.num_rows == 3, "ToS must expose image rows")
    require(tos.num_columns == 3, "ToS must expose image columns")
    require(valued_tree_of_shapes.reconstruct_from_node_altitudes().shape == (3, 3), "valued_tree ToS reconstructFromNodeAltitudes shape")
    require(valued_tree_of_shapes.node_altitudes.dtype == np.uint16, "valued_tree ToS exact altitude dtype")
    require(valued_tree_of_shapes.reconstruct_from_node_altitudes().dtype == np.uint16, "valued_tree ToS reconstruction dtype")

    # The default convention publishes unchanged 8-bit source levels over the
    # canonical 4/8 complementary-grid immersion.
    tos_source = np.array([[1, 2, 1], [2, 3, 2], [1, 2, 1]], dtype=np.uint8)
    default_specification = mmcfilters.TopographicConvention()
    require(
        isinstance(default_specification.immersion, mmcfilters.CanonicalComplementaryGridImmersion),
        "default convention must select the canonical complementary grid",
    )
    require(
        default_specification.immersion.pairing == mmcfilters.ComplementaryPairing.MIN4_MAX8,
        "default convention must select minimum-4/maximum-8 connectivity",
    )
    require(
        default_specification.domain_extension == mmcfilters.TopographicDomainExtension.NONE,
        "default convention must not pad the source domain",
    )
    require(default_specification.infinity_pixel == 0, "default convention infinity pixel must be zero")
    require(
        default_specification.altitude_encoding == mmcfilters.TopographicAltitudeEncoding.UINT8,
        "default convention must publish uint8 altitudes",
    )
    default_tos = mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(tos_source)
    require(default_tos.node_altitudes.dtype == np.uint8, "default ToS must publish uint8 altitudes")
    require(
        default_tos.topographic_convention.altitude_encoding == mmcfilters.TopographicAltitudeEncoding.UINT8,
        "default ToS must declare the 8-bit altitude encoding",
    )
    require(
        default_tos.topographic_convention.domain_extension == mmcfilters.TopographicDomainExtension.NONE,
        "default ToS must not pad the source domain",
    )
    require(default_tos.topographic_convention.infinity_pixel == 0, "default ToS infinity pixel must be zero")
    require(
        default_tos.reconstruct_from_node_altitudes().tolist() == tos_source.tolist(),
        "default ToS reconstruction must reproduce the source image",
    )
    # The canonical immersion is resolved, so the retained convention exposes
    # explicit adjacencies bound to the source domain.
    default_adjacencies = default_tos.topographic_convention.immersion.complementary_adjacencies
    require(default_adjacencies.min_adjacency.size == 5, "resolved default minimum adjacency must be 4-connected")
    require(default_adjacencies.max_adjacency.size == 9, "resolved default maximum adjacency must be 8-connected")

    # The 8-bit encoding is exact: it is the doubled hierarchy halved.
    doubled_tos = mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(
        tos_source, complementary_convention(3, 3, 1.0, 1.5)
    )
    require(default_tos.num_nodes == doubled_tos.num_nodes, "8-bit ToS node count must match the doubled hierarchy")
    require(
        (2 * default_tos.node_altitudes.astype(np.uint16)).tolist() == doubled_tos.node_altitudes.tolist(),
        "8-bit ToS altitudes must be the doubled altitudes halved",
    )

    # A self-dual span immersion may place its exterior median on a half level.
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(
            tos_source,
            mmcfilters.TopographicConvention(
                mmcfilters.SelfDualSpanImmersion(),
                mmcfilters.TopographicDomainExtension.EXTERIOR_RING,
            ),
        ),
        "self-dual span immersion must reject the 8-bit altitude encoding",
    )

    inverse_pairing_tos = mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(
        tos_source,
        mmcfilters.TopographicConvention(
            mmcfilters.CanonicalComplementaryGridImmersion(mmcfilters.ComplementaryPairing.MIN8_MAX4)
        ),
    )
    inverse_pairing_adjacencies = inverse_pairing_tos.topographic_convention.immersion.complementary_adjacencies
    require(inverse_pairing_adjacencies.min_adjacency.size == 9, "resolved MIN8_MAX4 minimum adjacency must be 8-connected")
    require(inverse_pairing_adjacencies.max_adjacency.size == 5, "resolved MIN8_MAX4 maximum adjacency must be 4-connected")
    tos_boundary_names, tos_boundary_attrs = mmcfilters.Attribute.compute_attributes(
        valued_tree_of_shapes, [mmcfilters.Attribute.Group.BOUNDARY]
    )
    tos_topology_boundary_names, tos_topology_boundary_attrs = mmcfilters.Attribute.compute_topology_attributes(
        valued_tree_of_shapes, [mmcfilters.Attribute.Group.BOUNDARY]
    )
    require("BITQUAD_AREA" in tos_boundary_names, "valued_tree ToS BOUNDARY group must expose Bitquad through compute_attributes")
    require("CONTOUR_SIDE_SOUTH" in tos_boundary_names, "valued_tree ToS BOUNDARY group must expose contour sides")
    require(tos_boundary_names == tos_topology_boundary_names, "valued_tree ToS topology BOUNDARY names must match the full pipeline")
    require(
        tos_boundary_attrs.shape == (valued_tree_of_shapes.num_internal_node_slots, len(tos_boundary_names)),
        "valued_tree ToS BOUNDARY shape must use internal node-id space",
    )
    require(
        np.allclose(tos_boundary_attrs, tos_topology_boundary_attrs, equal_nan=True),
        "valued_tree ToS topology BOUNDARY route must match the full pipeline",
    )
    tos_min8_max4 = mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(
        np.array([[1, 2], [3, 0]], dtype=np.uint8),
        complementary_convention(2, 2, 1.5, 1.0),
    )
    inverse_adjacencies = tos_min8_max4.topographic_convention.immersion.complementary_adjacencies
    require(inverse_adjacencies.min_adjacency.radius == 1.5, "Min8cMax4c minimum radius")
    require(inverse_adjacencies.max_adjacency.radius == 1.0, "Min8cMax4c maximum radius")
    require(inverse_adjacencies.min_adjacency.size == 9, "Min8cMax4c minimum adjacency size")
    require(inverse_adjacencies.max_adjacency.size == 5, "Min8cMax4c maximum adjacency size")
    unpadded_convention = complementary_convention(
        2, 3, 1.0, 1.5, mmcfilters.TopographicDomainExtension.NONE
    )
    require(
        unpadded_convention.domain_extension == mmcfilters.TopographicDomainExtension.NONE,
        "topographic convention must expose the domain extension",
    )
    unpadded_input = np.array([[0, 2, 1], [2, 1, 0]], dtype=np.uint8)
    unpadded_tos = mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(
        unpadded_input,
        unpadded_convention,
    )
    require(
        unpadded_tos.num_rows == 2 and unpadded_tos.num_columns == 3,
        "unpadded ToS must publish the original image domain",
    )
    require(
        unpadded_tos.reconstruct_from_node_altitudes().tolist() == (2 * unpadded_input.astype(np.uint16)).tolist(),
        "unpadded ToS reconstruction",
    )
    virtual_root_tos = mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(
        np.array([[1, 1], [0, 0]], dtype=np.uint8),
        self_dual_convention(),
    )
    virtual_root = virtual_root_tos.root
    upper_shape = virtual_root_tos.smallest_node(0)
    lower_shape = virtual_root_tos.smallest_node(2)
    require(virtual_root_tos.num_nodes == 3, "SelfDual ToS must preserve the virtual root")
    require(virtual_root_tos.proper_part_cardinality(virtual_root) == 0, "virtual root must have an empty direct proper part")
    require(virtual_root_tos.has_empty_proper_part(virtual_root), "virtual root must be derived as structural")
    require(virtual_root_tos.num_children(virtual_root) == 2, "virtual root must retain both child shapes")
    require(
        virtual_root_tos.reconstruct_from_node_altitudes().tolist() == [[2, 2], [0, 0]],
        "virtual-root ToS reconstruction",
    )
    require(virtual_root_tos.node_altitude(virtual_root) == 1, "virtual-root ToS exact half-level altitude")
    require(
        virtual_root_tos.node_altitude(virtual_root) != virtual_root_tos.node_altitude(lower_shape),
        "virtual-root ToS must not contain an equal-altitude lower edge",
    )
    native_partial_partition = (
        mmcfilters.MorphologicalTreeFactory.create_from_native_topology(
            [0, 0, 0],
            [1, 1, 2, 2],
            np.array([0, 1, 0], dtype=np.uint8),
            0,
            2,
            2,
            semantics=mmcfilters.MorphologicalTreeSemantics(
                kind=mmcfilters.MorphologicalTreeKind.TREE_OF_SHAPES,
                construction_context=mmcfilters.TopographicConvention(
                    mmcfilters.ComplementaryGridImmersion(
                        mmcfilters.ComplementaryAdjacencies(
                            mmcfilters.RegularGridAdjacency2D(2, 2, 1.0),
                            mmcfilters.RegularGridAdjacency2D(2, 2, 1.5),
                        )
                    )
                ),
            ),
        )
    )
    require(
        native_partial_partition.proper_part_cardinality(0) == 0,
        "generic native partial-partition root may have no direct proper part",
    )
    require(
        native_partial_partition.num_children(0) == 2,
        "generic native partial-partition root children",
    )
    require(
        isinstance(native_partial_partition.topographic_convention, mmcfilters.TopographicConvention),
        "native partial-partition topographic convention",
    )
    require(
        native_partial_partition.topographic_convention.immersion.complementary_adjacencies.min_adjacency.size == 5,
        "native partial-partition minimum adjacency",
    )
    require(
        native_partial_partition.topographic_convention.immersion.complementary_adjacencies.max_adjacency.size == 9,
        "native partial-partition maximum adjacency",
    )
    require(
        native_partial_partition.reconstruct_from_node_altitudes().tolist()
        == [[1, 1], [0, 0]],
        "generic native partial-partition reconstruction",
    )
    generic_semantics = mmcfilters.MorphologicalTreeSemantics(
        node_altitude_order=mmcfilters.NodeAltitudeOrder.INCREASING,
    )
    generic_tree = mmcfilters.MorphologicalTreeFactory.create_from_native_topology(
        parent=[0, 0, 1],
        smallest_node_map=[2],
        node_altitudes=np.array([0, 1, 2], dtype=np.uint8),
        root=0,
        rows=1,
        columns=1,
        semantics=generic_semantics,
    )
    require(
        generic_tree.kind == mmcfilters.MorphologicalTreeKind.GENERIC,
        "generic native tree descriptive kind",
    )
    require(
        generic_tree.node_altitude_order
        == mmcfilters.NodeAltitudeOrder.INCREASING,
        "generic native tree altitude order",
    )
    require(
        isinstance(generic_tree.construction_context, mmcfilters.NoConstructionContext),
        "generic native tree records no construction context",
    )
    require(
        generic_tree.has_grid_domain_2d
        and generic_tree.grid_domain_2d.rows == 1
        and generic_tree.grid_domain_2d.columns == 1,
        "generic native tree explicit 2D domain",
    )
    require(
        generic_tree.has_empty_proper_part(0)
        and generic_tree.has_empty_proper_part(1)
        and not generic_tree.has_empty_proper_part(2),
        "generic native tree structural-node derivation",
    )
    require(
        not generic_tree.is_tree_of_partial_partitions(),
        "generic tree with empty proper parts is not a tree of partial partitions",
    )
    require_raises(
        generic_tree.validate_tree_of_partial_partitions,
        "tree-of-partial-partitions validation must reject empty proper parts",
    )
    for node_id in generic_tree.alive_node_ids:
        require(
            len(list(generic_tree.node_support(node_id))) > 0,
            "every committed generic node must have non-empty subtree support",
        )
    abstract_tree = mmcfilters.MorphologicalTreeFactory.create_from_native_topology(
        [0, 0, 0],
        [1, 2],
        np.array([10, 3, 20], dtype=np.uint8),
        0,
        semantics=mmcfilters.MorphologicalTreeSemantics(),
    )
    require(
        abstract_tree.has_grid_domain_2d is False,
        "abstract pixel domain must not invent grid metadata",
    )
    require(
        abstract_tree.grid_domain_2d is None,
        "abstract pixel domain optional grid",
    )
    require(
        abstract_tree.kind
        == mmcfilters.MorphologicalTreeKind.GENERIC,
        "canonical descriptive-kind property",
    )
    abstract_area = mmcfilters.Attribute.compute_single_attribute(
        abstract_tree,
        mmcfilters.Attribute.AREA,
    )
    require(
        int(abstract_area[0]) == 2,
        "support attributes must work on an abstract pixel domain",
    )
    abstract_gray_height = mmcfilters.Attribute.compute_single_attribute(
        abstract_tree,
        mmcfilters.Attribute.GRAY_LEVEL_HEIGHT,
    )
    require(
        abstract_gray_height.tolist() == [10.0, 0.0, 0.0],
        "unconstrained GRAY_LEVEL_HEIGHT must use both subtree extrema",
    )
    require_raises(
        lambda: abstract_tree.reconstruct_from_node_altitudes(),
        "abstract pixel domain must reject image reconstruction",
    )
    require_raises(
        lambda: mmcfilters.Attribute.compute_single_attribute(
            abstract_tree,
            mmcfilters.Attribute.BOX_WIDTH,
        ),
        "abstract pixel domain must reject geometric attributes",
    )
    require_raises(
        lambda: abstract_tree.num_rows,
        "abstract pixel domain must not expose invented rows",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_from_native_topology(
            [0, 0],
            [0],
            np.array([0, 1], dtype=np.uint8),
            0,
            1,
            1,
            semantics=mmcfilters.MorphologicalTreeSemantics(),
        ),
        "native factory must reject an attached leaf with empty subtree support",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_from_native_topology(
            [0, 0, 1],
            [2],
            np.array([0, 2, 1], dtype=np.uint8),
            0,
            1,
            1,
            semantics=mmcfilters.MorphologicalTreeSemantics(
                node_altitude_order=mmcfilters.NodeAltitudeOrder.INCREASING,
            ),
        ),
        "native factory altitude must satisfy the declared generic order",
    )
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_from_native_topology(
            [0, 0, 1],
            [2],
            np.array([0, 1, 1], dtype=np.uint8),
            0,
            1,
            1,
            semantics=mmcfilters.MorphologicalTreeSemantics(
                node_altitude_order=mmcfilters.NodeAltitudeOrder.INCREASING,
            ),
        ),
        "native factory altitude must reject equality in a strict order",
    )
    boundary_infinity_tos = mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(
        np.array([[0, 1]], dtype=np.uint8),
        self_dual_convention(infinity_pixel=14),
    )
    require(boundary_infinity_tos.root != -1, "SelfDual ToS must accept a custom infinity pixel on the outer boundary")
    internal_infinity_tos = mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(
        np.array([[0, 1]], dtype=np.uint8),
        self_dual_convention(infinity_pixel=6),
    )
    require(internal_infinity_tos.root != -1, "SelfDual ToS must accept an internal custom infinity pixel")
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(
            np.array([[0, 1]], dtype=np.uint8),
            self_dual_convention(infinity_pixel=20),
        ),
        "SelfDual ToS must reject an infinity pixel outside the interpolated domain",
    )
    single_tos = mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(np.array([[5]], dtype=np.uint8))
    single_valued_tree_of_shapes = mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(np.array([[5]], dtype=np.uint8))
    require(single_tos.num_rows == 1 and single_tos.num_columns == 1, "single-pixel default ToS dimensions")
    require(single_tos.root != -1, "single-pixel default ToS root")
    # The default convention publishes unchanged source levels, so the
    # reconstruction is the source pixel itself rather than its doubled unit.
    require(single_valued_tree_of_shapes.reconstruct_from_node_altitudes().tolist() == [[5]], "single-pixel default valued_tree ToS reconstruction")

    empty = np.empty((0, 0), dtype=np.uint8)
    require_raises(lambda: mmcfilters.MorphologicalTreeFactory.create_max_tree(empty), "empty max-tree must throw")
    require_raises(lambda: mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(empty), "empty tree of shapes must throw")
    require_raises(lambda: mmcfilters.MorphologicalTreeFactory.create_min_tree(empty), "empty min-tree must throw")
    require_raises(lambda: mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(empty), "empty valued tree of shapes must throw")

    # Concise call forms. Symbolic attribute names are accepted wherever an
    # Attribute value is, and resolve to exactly the same computation.
    concise_tree = mmcfilters.MorphologicalTreeFactory.create_max_tree(image, radius=1.5)
    require(
        np.array_equal(
            mmcfilters.Attribute.compute_single_attribute(concise_tree, "AREA"),
            mmcfilters.Attribute.compute_single_attribute(concise_tree, mmcfilters.Attribute.AREA),
        ),
        "symbolic attribute name must match the enum value",
    )
    require(
        np.array_equal(
            mmcfilters.Attribute.compute_single_topology_attribute(concise_tree, "AREA"),
            mmcfilters.Attribute.compute_single_topology_attribute(concise_tree, mmcfilters.Attribute.AREA),
        ),
        "symbolic name must match the enum value for topology attributes",
    )
    named_layout, named_values = mmcfilters.Attribute.compute_attributes(concise_tree, ["AREA", "VOLUME"])
    typed_layout, typed_values = mmcfilters.Attribute.compute_attributes(
        concise_tree, [mmcfilters.Attribute.AREA, mmcfilters.Attribute.VOLUME]
    )
    require(list(named_layout) == list(typed_layout), "symbolic attribute list must produce the same layout")
    require(np.array_equal(named_values, typed_values), "symbolic attribute list must produce the same values")

    # Group names and mixed sequences are accepted too.
    group_layout, _ = mmcfilters.Attribute.compute_attributes(concise_tree, ["BOUNDARY"])
    typed_group_layout, _ = mmcfilters.Attribute.compute_attributes(concise_tree, [mmcfilters.Attribute.Group.BOUNDARY])
    require(list(group_layout) == list(typed_group_layout), "symbolic group name must match the enum group")
    mixed_layout, _ = mmcfilters.Attribute.compute_attributes(concise_tree, [mmcfilters.Attribute.AREA, "VOLUME"])
    require(list(mixed_layout) == ["AREA", "VOLUME"], "mixed enum and symbolic sequences must be accepted")

    require(
        mmcfilters.Attribute.describe("AREA") == mmcfilters.Attribute.describe(mmcfilters.Attribute.AREA),
        "describe must accept a symbolic name",
    )
    require(
        mmcfilters.Attribute.requirements("AREA") == mmcfilters.Attribute.requirements(mmcfilters.Attribute.AREA),
        "requirements must accept a symbolic name",
    )
    require(
        np.array_equal(
            mmcfilters.Attribute.compute_attribute_mapping(concise_tree, "AREA"),
            mmcfilters.Attribute.compute_attribute_mapping(concise_tree, mmcfilters.Attribute.AREA),
        ),
        "compute_attribute_mapping must accept a symbolic name",
    )

    # Matching is exact, and a rejected name reports the near matches.
    require_raises(lambda: mmcfilters.Attribute.compute_single_attribute(concise_tree, "area"), "symbolic names must be case sensitive")
    require_raises(lambda: mmcfilters.Attribute.compute_single_attribute(concise_tree, "NOT_AN_ATTRIBUTE"), "unknown names must be rejected")
    require_raises(lambda: mmcfilters.Attribute.compute_attributes(concise_tree, "AREA"), "a bare str must not stand in for a sequence")

    # The self-dual span convention is available as a single call.
    require(
        mmcfilters.self_dual_span_convention().altitude_encoding == mmcfilters.TopographicAltitudeEncoding.EXACT_DOUBLED,
        "the self-dual span helper must default to doubled units",
    )
    unpadded_eight_bit = mmcfilters.self_dual_span_convention(
        domain_extension=mmcfilters.TopographicDomainExtension.NONE,
        altitude_encoding=mmcfilters.TopographicAltitudeEncoding.UINT8,
    )
    require(isinstance(unpadded_eight_bit.immersion, mmcfilters.SelfDualSpanImmersion), "the helper must select the span immersion")

    # Convention fields may be passed directly to the factory.
    field_tree = mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(
        tos_source,
        immersion=mmcfilters.SelfDualSpanImmersion(),
        domain_extension=mmcfilters.TopographicDomainExtension.NONE,
        altitude_encoding=mmcfilters.TopographicAltitudeEncoding.UINT8,
    )
    convention_tree = mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(tos_source, unpadded_eight_bit)
    require(field_tree.num_nodes == convention_tree.num_nodes, "convention fields must build the same tree as the convention")
    require(
        np.array_equal(field_tree.node_altitudes, convention_tree.node_altitudes),
        "convention fields must publish the same altitudes as the convention",
    )
    require(field_tree.node_altitudes.dtype == np.uint8, "unpadded self-dual span must publish uint8 altitudes")
    require_raises(
        lambda: mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(
            tos_source,
            mmcfilters.TopographicConvention(),
            altitude_encoding=mmcfilters.TopographicAltitudeEncoding.UINT8,
        ),
        "a complete convention and individual fields must not be combined",
    )

    print("python NodeId API ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
