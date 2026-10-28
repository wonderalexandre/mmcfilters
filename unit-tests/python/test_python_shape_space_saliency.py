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
        return any(
            name == "mmcfilters" or name.startswith("mmcfilters.")
            for name in known_modules
        )

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


def require_raises(fn, message: str, expected_exception=ValueError):
    try:
        fn()
    except expected_exception:
        return
    raise RuntimeError(message)


def make_nested_chain(mmcfilters, radius=1.0):
    parent = [4, 5, 6, 6, 5, 6, 6]
    altitude = np.array([3, 2, 1, 1, 3, 2, 1], dtype=np.uint8)
    return mmcfilters.MorphologicalTreeFactory.create_from_higra_parent(
        parent,
        altitude,
        1,
        4,
        mmcfilters.MorphologicalTreeKind.MAX_TREE,
        radius,
    )


def check_composite_result_storage(mmcfilters):
    import gc

    retained = []
    for dtype in (np.float32, np.float64):
        for polarity in (mmcfilters.ShapeSpaceExtremaPolarity.MINIMA, mmcfilters.ShapeSpaceExtremaPolarity.MAXIMA):
            tree = make_nested_chain(mmcfilters)
            scores = np.array([1, 5, 2], dtype=dtype)
            extinction = mmcfilters.ShapeSpaceSaliency.compute_extinction_values(tree, scores, polarity)
            result = mmcfilters.ShapeSpaceSaliency.compute(tree, scores, polarity)
            require(isinstance(extinction, mmcfilters.ShapeSpaceExtinctionResult), "typed extinction result")
            require(isinstance(result, mmcfilters.ShapeSpaceSaliencyResult), "typed composite saliency result")
            np.testing.assert_array_equal(result.node_scores, extinction.node_scores)
            require(result.edge_map is result.edge_map, "composite edge-map access reuses the owned component")
            expected_pixels = result.edge_map.to_pixel_image()
            for owner in (extinction, result):
                values = owner.node_scores
                require(values.dtype == dtype and values.shape == (tree.num_internal_node_slots,),
                        "shape-space scores retain precision and dense node-slot layout")
                require(not values.flags.writeable and values.base is owner, "node scores retain the immutable result")
                require(np.shares_memory(values, owner.node_scores), "score access does not copy the buffer")
                require_raises(lambda: values.__setitem__(0, 99), "scores reject writes")
                require_raises(lambda: values.setflags(write=True), "score writes cannot be re-enabled")
                require(isinstance(owner.extrema, tuple) and owner.extrema is owner.extrema,
                        "extrema are cached immutable records")
                records = tuple((item.representative, item.birth_level, item.death_level, item.extinction)
                                for item in owner.extrema)
                for item in owner.extrema:
                    require(isinstance(item, mmcfilters.ShapeSpaceExtremum), "each extremum has a typed record")
                    require_raises(lambda: setattr(item, "extinction", -1), "extremum fields are read-only", AttributeError)
                require_raises(lambda: setattr(owner, "node_scores", scores), "score properties are read-only", AttributeError)
                retained.append((values, values.copy(), owner.extrema, records))

            edge_map = result.edge_map
            edge_values = edge_map.values
            np.testing.assert_array_equal(edge_map.values,
                                          mmcfilters.ShapeSpaceSaliency.project_contour_scores(tree, extinction.node_scores).values)
            require(edge_values.base is edge_map and not edge_values.flags.writeable, "nested map views retain their component")
            expected_edges = edge_values.copy()
            require_raises(lambda: setattr(result, "edge_map", edge_map), "the composite edge-map property is read-only", AttributeError)
            scores[:] = -99
            tree.merge_node_into_parent(tree.leaves[0])
            for values, expected, records, expected_records in retained:
                np.testing.assert_array_equal(values, expected)
                require(tuple((item.representative, item.birth_level, item.death_level, item.extinction)
                              for item in records) == expected_records, "records remain snapshots after tree edits")
            # Keep only the nested map and array references after deleting the composite.
            del tree, result, extinction, owner, scores
            gc.collect()
            np.testing.assert_array_equal(edge_map.values, expected_edges)
            np.testing.assert_array_equal(edge_map.to_pixel_image(), expected_pixels)
            del edge_map
            gc.collect()
            np.testing.assert_array_equal(edge_values, expected_edges)

    del values, item
    gc.collect()
    for values, expected, records, expected_records in retained:
        np.testing.assert_array_equal(values, expected)
        require(tuple((item.representative, item.birth_level, item.death_level, item.extinction)
                      for item in records) == expected_records, "record lifetime is independent of the result object")

    def detached_map():
        local_tree = make_nested_chain(mmcfilters)
        return mmcfilters.ShapeSpaceSaliency.compute(
            local_tree, np.array([1, 4, 10], dtype=np.float64),
            mmcfilters.ShapeSpaceExtremaPolarity.MINIMA,
        ).edge_map

    nested = detached_map()
    gc.collect()
    np.testing.assert_array_equal(nested.values, [9, 0, 0])
    detached_values = nested.values
    del nested
    gc.collect()
    np.testing.assert_array_equal(detached_values, [9, 0, 0])

    singleton = mmcfilters.MorphologicalTreeFactory.create_max_tree(np.array([[7]], dtype=np.uint8))
    single = mmcfilters.ShapeSpaceSaliency.compute(singleton, np.array([2], dtype=np.float32),
                                                mmcfilters.ShapeSpaceExtremaPolarity.MINIMA)
    np.testing.assert_array_equal(single.node_scores, [0])
    require(len(single.edge_map) == 0, "a one-pixel Xu result has no graph edges")
    require(len(single.extrema) == 1 and single.extrema[0].extinction == 0, "a constant attribute has a zero extinction interval")


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: test_python_shape_space_saliency.py <build-dir>")

    mmcfilters = load_native_module(pathlib.Path(sys.argv[1]).resolve())
    check_composite_result_storage(mmcfilters)

    require(hasattr(mmcfilters, "ShapeSpaceExtremaPolarity"), "polarity enum export")
    require(hasattr(mmcfilters, "ShapeSpaceSaliency"), "shape-space saliency export")
    require(
        not hasattr(mmcfilters.ShapeSpaceSaliency, "computeExtinctionValues"),
        "camelCase compatibility alias must not be exposed",
    )
    require(
        not hasattr(mmcfilters.ShapeSpaceSaliency, "projectContourScores"),
        "camelCase compatibility alias must not be exposed",
    )

    tree = make_nested_chain(mmcfilters)
    scores32 = np.array([1.0, 4.0, 10.0], dtype=np.float32)
    extinction = mmcfilters.ShapeSpaceSaliency.compute_extinction_values(
        tree=tree,
        attribute=scores32,
        polarity=mmcfilters.ShapeSpaceExtremaPolarity.MINIMA,
    )
    require(isinstance(extinction, mmcfilters.ShapeSpaceExtinctionResult), "typed extinction result")
    require(extinction.node_scores.dtype == np.float32, "float32 extinction dtype")
    require(
        np.array_equal(extinction.node_scores, np.array([9.0, 0.0, 0.0], dtype=np.float32)),
        "float32 extinction values",
    )
    require(len(extinction.extrema) == 1, "nested-chain minimum count")
    extremum = extinction.extrema[0]
    require(isinstance(extremum, mmcfilters.ShapeSpaceExtremum), "typed extremum record")
    require(extremum.representative == 0, "minimum representative")
    require(abs(extremum.birth_level - 1.0) < 1e-6, "minimum birth level")
    require(abs(extremum.death_level - 10.0) < 1e-6, "minimum death level")
    require(abs(extremum.extinction - 9.0) < 1e-6, "minimum extinction")

    projection = mmcfilters.ShapeSpaceSaliency.project_contour_scores(
        tree,
        extinction.node_scores,
    )
    require(projection.values.dtype == np.float32, "float32 projection dtype")
    require(
        np.array_equal(projection.values, np.array([9.0, 0.0, 0.0], dtype=np.float32)),
        "maximum-on-contours projection",
    )

    one_shot = mmcfilters.ShapeSpaceSaliency.compute(
        tree=tree,
        attribute=scores32,
        polarity=mmcfilters.ShapeSpaceExtremaPolarity.MINIMA,
    )
    require(isinstance(one_shot, mmcfilters.ShapeSpaceSaliencyResult), "typed one-shot saliency result")
    require(np.array_equal(one_shot.node_scores, extinction.node_scores), "one-shot node scores")
    require(np.array_equal(one_shot.edge_map.values, projection.values), "one-shot edge map")

    scores64 = np.array([1.0, 5.0, 2.0], dtype=np.float64)
    maxima = mmcfilters.ShapeSpaceSaliency.compute(
        tree,
        scores64,
        mmcfilters.ShapeSpaceExtremaPolarity.MAXIMA,
        radius=1.0,
    )
    require(maxima.node_scores.dtype == np.float64, "float64 extinction dtype")
    require(
        np.array_equal(maxima.node_scores, np.array([0.0, 4.0, 0.0], dtype=np.float64)),
        "maximum-oriented extinction values",
    )
    require(
        np.array_equal(maxima.edge_map.values, np.array([0.0, 4.0, 0.0], dtype=np.float64)),
        "maximum-oriented contour values",
    )

    require_raises(
        lambda: mmcfilters.ShapeSpaceSaliency.compute_extinction_values(
            tree,
            np.array([1.0, 4.0], dtype=np.float32),
            mmcfilters.ShapeSpaceExtremaPolarity.MINIMA,
        ),
        "short score arrays must be rejected",
    )
    require_raises(
        lambda: mmcfilters.ShapeSpaceSaliency.compute_extinction_values(
            tree,
            np.array([1, 4, 10], dtype=np.int32),
            mmcfilters.ShapeSpaceExtremaPolarity.MINIMA,
        ),
        "integer score arrays must be rejected",
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
