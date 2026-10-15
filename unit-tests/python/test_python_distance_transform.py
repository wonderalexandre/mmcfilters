#!/usr/bin/env python3

import importlib.util
import os
import pathlib
import sys

import numpy as np


CONTRACT_CHECKED = os.environ.get("MMCFILTERS_CONTRACT_MODE", "CHECKED").upper() == "CHECKED"


def require(condition: bool, message: str):
    if not condition:
        raise RuntimeError(message)


def require_raises(fn, message: str):
    if not CONTRACT_CHECKED:
        return
    try:
        fn()
    except Exception:
        return
    raise RuntimeError(message)


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
        return any(
            name == "mmcfilters" or name.startswith("mmcfilters.")
            for name in known_modules
        )

    sys.meta_path = [finder for finder in sys.meta_path if not handles_mmcfilters(finder)]
    spec = importlib.util.spec_from_file_location(
        "mmcfilters", package_init, submodule_search_locations=[str(package_dir)]
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules["mmcfilters"] = module
    spec.loader.exec_module(module)
    return module


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: test_python_distance_transform.py <build-dir>")

    mmcfilters = load_native_module(pathlib.Path(sys.argv[1]).resolve())
    image = np.full((5, 5), 7, dtype=np.uint8)
    tree = mmcfilters.MorphologicalTreeFactory.create_max_tree(image)
    expected = np.array(
        [
            [0, 0, 0, 0, 0],
            [0, 1, 1, 1, 0],
            [0, 1, 4, 1, 0],
            [0, 1, 1, 1, 0],
            [0, 0, 0, 0, 0],
        ],
        dtype=np.int64,
    )

    exact_field = mmcfilters.DistanceTransform.compute_node_field(
        tree, tree.root, "exact"
    )
    require(exact_field.node_id == tree.root, "exact field node id")
    require(
        exact_field.accuracy == mmcfilters.DistanceTransformAccuracy.EXACT,
        "exact field accuracy",
    )
    require(len(exact_field.support_pixels) == 25, "exact compact support size")
    require(
        np.array_equal(exact_field.to_dense(), expected),
        "exact dense squared-distance field",
    )

    approximate_field = mmcfilters.DistanceTransform.compute_node_field(
        tree, tree.root, mmcfilters.DistanceTransformAccuracy.APPROXIMATE
    )
    require(
        approximate_field.accuracy
        == mmcfilters.DistanceTransformAccuracy.APPROXIMATE,
        "approximate field accuracy",
    )
    require(
        np.array_equal(approximate_field.to_dense(), expected),
        "approximate dense squared-distance field",
    )

    exact_profiles = mmcfilters.DistanceTransform.compute_profiles(tree)
    require(
        exact_profiles.accuracy == mmcfilters.DistanceTransformAccuracy.EXACT,
        "default profile accuracy",
    )
    require(
        exact_profiles.histogram(tree.root) == [(0, 16), (1, 8), (4, 1)],
        "exact root histogram",
    )
    require(
        exact_profiles.survival_profile(tree.root) == [(0, 25), (1, 9), (4, 1)],
        "exact root survival profile",
    )

    approximate_profiles = mmcfilters.DistanceTransform.compute_profiles(
        tree, "approximate"
    )
    require(
        approximate_profiles.histogram(tree.root)
        == exact_profiles.histogram(tree.root),
        "constant-square approximate profile",
    )

    general_image = np.array(
        [[3, 1, 4], [1, 5, 9], [2, 6, 5]], dtype=np.uint8
    )
    hierarchy_cases = {
        "max-tree": mmcfilters.MorphologicalTreeFactory.create_max_tree(general_image),
        "min-tree": mmcfilters.MorphologicalTreeFactory.create_min_tree(general_image),
        "extremum residual tree": (
            mmcfilters.MorphologicalTreeFactory.create_extremum_residual_tree(
                general_image, radius=1.0
            )
        ),
        "tree of shapes": (
            mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(general_image)
        ),
    }
    for hierarchy_name, hierarchy in hierarchy_cases.items():
        for accuracy in (
            mmcfilters.DistanceTransformAccuracy.EXACT,
            mmcfilters.DistanceTransformAccuracy.APPROXIMATE,
        ):
            field = mmcfilters.DistanceTransform.compute_node_field(
                hierarchy, hierarchy.root, accuracy
            )
            root_support_size = sum(
                1 for _ in hierarchy.node_support(hierarchy.root)
            )
            require(
                len(field.support_pixels) == root_support_size,
                f"{hierarchy_name} root field support",
            )
            require(
                field.to_dense().shape == general_image.shape,
                f"{hierarchy_name} dense field shape",
            )
            profiles = mmcfilters.DistanceTransform.compute_profiles(
                hierarchy, accuracy
            )
            for node_id in hierarchy.alive_node_ids:
                support_size = sum(1 for _ in hierarchy.node_support(node_id))
                require(
                    sum(count for _, count in profiles.histogram(node_id))
                    == support_size,
                    f"{hierarchy_name} profile support count for node {node_id}",
                )

    exact_names, exact_values = mmcfilters.Attribute.compute_topology_attributes(
        tree,
        [
            mmcfilters.Attribute.MAX_DIST_EXACT,
            mmcfilters.Attribute.MAX_SQUARED_DIST_EXACT,
        ],
        dtype=np.float64,
    )
    exact_root = dict(zip(exact_names, exact_values[tree.root]))
    require(exact_root["MAX_DIST_EXACT"] == 2.0, "exact maximum distance")
    require(
        exact_root["MAX_SQUARED_DIST_EXACT"] == 4.0,
        "exact maximum squared distance",
    )

    require_raises(
        lambda: mmcfilters.DistanceTransform.compute_node_field(
            tree, tree.root, "unknown"
        ),
        "unknown accuracy must fail",
    )
    require_raises(
        lambda: exact_profiles.histogram(-1),
        "invalid profile node id must fail",
    )
    require(
        not hasattr(mmcfilters, "EXACT"),
        "accuracy values must not pollute the top-level module",
    )
    require(
        not hasattr(mmcfilters, "APPROXIMATE"),
        "accuracy values must stay nested in their enum",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
