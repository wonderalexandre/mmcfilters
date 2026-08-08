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


def require_raises(fn, message: str):
    try:
        fn()
    except ValueError:
        return
    raise RuntimeError(message)


def make_nested_chain(mmcfilters, radius=1.0):
    parent = [4, 5, 6, 6, 5, 6, 6]
    altitude = np.array([3, 2, 1, 1, 3, 2, 1], dtype=np.uint8)
    return mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
        parent,
        altitude,
        1,
        4,
        mmcfilters.MorphologicalTreeKind.MAX_TREE,
        radius,
    )


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: test_python_shape_space_saliency.py <build-dir>")

    mmcfilters = load_native_module(pathlib.Path(sys.argv[1]).resolve())

    require(hasattr(mmcfilters, "ShapeSpaceExtremaPolarity"), "polarity enum export")
    require(hasattr(mmcfilters, "ShapeSpaceSaliency"), "shape-space saliency export")
    require(
        not hasattr(mmcfilters.ShapeSpaceSaliency, "compute_extinction_values"),
        "snake-case compatibility alias must not be exposed",
    )
    require(
        not hasattr(mmcfilters.ShapeSpaceSaliency, "project_contour_scores"),
        "snake-case compatibility alias must not be exposed",
    )

    tree = make_nested_chain(mmcfilters)
    scores32 = np.array([1.0, 4.0, 10.0], dtype=np.float32)
    extinction = mmcfilters.ShapeSpaceSaliency.computeExtinctionValues(
        tree=tree,
        attribute=scores32,
        polarity=mmcfilters.ShapeSpaceExtremaPolarity.Minima,
    )
    require(set(extinction) == {"extrema", "nodeScores"}, "extinction result keys")
    require(extinction["nodeScores"].dtype == np.float32, "float32 extinction dtype")
    require(
        np.array_equal(extinction["nodeScores"], np.array([9.0, 0.0, 0.0], dtype=np.float32)),
        "float32 extinction values",
    )
    require(len(extinction["extrema"]) == 1, "nested-chain minimum count")
    extremum = extinction["extrema"][0]
    require(
        set(extremum) == {"representative", "birthLevel", "deathLevel", "extinction"},
        "extremum result keys",
    )
    require(extremum["representative"] == 0, "minimum representative")
    require(abs(extremum["birthLevel"] - 1.0) < 1e-6, "minimum birth level")
    require(abs(extremum["deathLevel"] - 10.0) < 1e-6, "minimum death level")
    require(abs(extremum["extinction"] - 9.0) < 1e-6, "minimum extinction")

    projection = mmcfilters.ShapeSpaceSaliency.projectContourScores(
        tree,
        extinction["nodeScores"],
    )
    require(projection["values"].dtype == np.float32, "float32 projection dtype")
    require(
        np.array_equal(projection["values"], np.array([9.0, 0.0, 0.0], dtype=np.float32)),
        "maximum-on-contours projection",
    )

    one_shot = mmcfilters.ShapeSpaceSaliency.compute(
        tree=tree,
        attribute=scores32,
        polarity=mmcfilters.ShapeSpaceExtremaPolarity.Minima,
    )
    require(set(one_shot) == {"extrema", "nodeScores", "edgeMap"}, "one-shot result keys")
    require(np.array_equal(one_shot["nodeScores"], extinction["nodeScores"]), "one-shot node scores")
    require(np.array_equal(one_shot["edgeMap"]["values"], projection["values"]), "one-shot edge map")

    scores64 = np.array([1.0, 5.0, 2.0], dtype=np.float64)
    maxima = mmcfilters.ShapeSpaceSaliency.compute(
        tree,
        scores64,
        mmcfilters.ShapeSpaceExtremaPolarity.Maxima,
        radius=1.0,
    )
    require(maxima["nodeScores"].dtype == np.float64, "float64 extinction dtype")
    require(
        np.array_equal(maxima["nodeScores"], np.array([0.0, 4.0, 0.0], dtype=np.float64)),
        "maximum-oriented extinction values",
    )
    require(
        np.array_equal(maxima["edgeMap"]["values"], np.array([0.0, 4.0, 0.0], dtype=np.float64)),
        "maximum-oriented contour values",
    )

    require_raises(
        lambda: mmcfilters.ShapeSpaceSaliency.computeExtinctionValues(
            tree,
            np.array([1.0, 4.0], dtype=np.float32),
            mmcfilters.ShapeSpaceExtremaPolarity.Minima,
        ),
        "short score arrays must be rejected",
    )
    require_raises(
        lambda: mmcfilters.ShapeSpaceSaliency.computeExtinctionValues(
            tree,
            np.array([1, 4, 10], dtype=np.int32),
            mmcfilters.ShapeSpaceExtremaPolarity.Minima,
        ),
        "integer score arrays must be rejected",
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
