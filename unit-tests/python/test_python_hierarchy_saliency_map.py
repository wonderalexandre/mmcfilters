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


def require_radius_rejected(fn, radius: float, context: str):
    if not np.isfinite(radius):
        expected_reason = "radius must be finite"
    elif radius < 1.0:
        expected_reason = "radius must be at least 1.0"
    else:
        expected_reason = "radius exceeds the supported integer stencil range"
    try:
        fn(radius)
    except ValueError as exc:
        error = str(exc)
        require(context in error, f"{context} radius error must identify the entrypoint: {error}")
        require(expected_reason in error, f"{context} radius error must explain the invalid value: {error}")
        return
    except Exception as exc:
        raise RuntimeError(
            f"{context} invalid radius must raise ValueError, got {type(exc).__name__}: {exc}"
        ) from exc
    raise RuntimeError(f"{context} must reject radius={radius!r}")


def qfz_labels(edge_map, threshold):
    num_vertices = edge_map["numRows"] * edge_map["numCols"]
    parent = list(range(num_vertices))

    def find_root(node):
        root = node
        while parent[root] != root:
            root = parent[root]
        while parent[node] != node:
            nxt = parent[node]
            parent[node] = root
            node = nxt
        return root

    for source, target, value in zip(edge_map["sources"], edge_map["targets"], edge_map["values"]):
        if value < threshold:
            source_root = find_root(int(source))
            target_root = find_root(int(target))
            if source_root != target_root:
                parent[target_root] = source_root

    return np.array([find_root(i) for i in range(num_vertices)], dtype=np.int32)


def make_three_pixel_tree(mmcfilters, kind, radius=1.0):
    parent = [3, 4, 4, 5, 5, 5]
    altitude = np.array([5, 4, 4, 5, 4, 0], dtype=np.uint8)
    return mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
        parent,
        altitude,
        1,
        3,
        kind,
        radius,
    )


def make_three_pixel_min_tree(mmcfilters):
    parent = [3, 4, 4, 5, 5, 5]
    altitude = np.array([0, 1, 1, 0, 1, 5], dtype=np.uint8)
    return mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
        parent,
        altitude,
        1,
        3,
        mmcfilters.MorphologicalTreeKind.MIN_TREE,
        1.0,
    )


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: test_python_hierarchy_saliency_map.py <build-dir>")

    build_dir = pathlib.Path(sys.argv[1]).resolve()
    mmcfilters = load_native_module(build_dir)
    require(
        not hasattr(mmcfilters.HierarchySaliencyMap, "computeAltitudeEdgeMap"),
        "raw altitude edge map must not be exposed in Python",
    )
    require(
        not hasattr(mmcfilters.HierarchySaliencyMap, "computeNodeScoreEdgeMap"),
        "node-score compatibility alias must not be exposed in Python",
    )

    weighted = make_three_pixel_tree(mmcfilters, mmcfilters.MorphologicalTreeKind.MAX_TREE)

    scores32 = np.array([10.0, 20.0, 30.0], dtype=np.float32)
    invalid_radii = (
        float("nan"),
        float("inf"),
        float("-inf"),
        -1.0,
        0.0,
        0.5,
        float(np.nextafter(np.float64(1.0), np.float64(0.0))),
        float(np.finfo(np.float64).max),
    )
    radius_entrypoints = (
        (
            "MorphologicalTreeFactory.createFromHigraParent",
            lambda radius: make_three_pixel_tree(
                mmcfilters,
                mmcfilters.MorphologicalTreeKind.MAX_TREE,
                radius=radius,
            ),
        ),
        (
            "HierarchySaliencyMap.computeSaliencyEdgeMap",
            lambda radius: mmcfilters.HierarchySaliencyMap.computeSaliencyEdgeMap(
                weighted,
                scores32,
                radius=radius,
            ),
        ),
        (
            "HierarchySaliencyMap.computeTopologicalLevelEdgeMap",
            lambda radius: mmcfilters.HierarchySaliencyMap.computeTopologicalLevelEdgeMap(
                weighted,
                radius=radius,
            ),
        ),
        (
            "HierarchySaliencyMap.computeNormalizedAltitudeEdgeMap",
            lambda radius: mmcfilters.HierarchySaliencyMap.computeNormalizedAltitudeEdgeMap(
                weighted,
                radius=radius,
            ),
        ),
        (
            "ComponentTreePartitionHierarchyAdapter.validate",
            lambda radius: mmcfilters.ComponentTreePartitionHierarchyAdapter.validate(
                weighted,
                radius=radius,
            ),
        ),
        (
            "ComponentTreePartitionHierarchyAdapter.computeSaliencyEdgeMap",
            lambda radius: mmcfilters.ComponentTreePartitionHierarchyAdapter.computeSaliencyEdgeMap(
                weighted,
                np.array([1, 1, 2], dtype=np.int32),
                radius=radius,
            ),
        ),
        (
            "HierarchySaliencyMapProjection.nodeContourEdges",
            lambda radius: mmcfilters.HierarchySaliencyMapProjection.nodeContourEdges(
                weighted,
                radius=radius,
            ),
        ),
        (
            "HierarchySaliencyMapProjection.computeIncrementalNodeContours",
            lambda radius: mmcfilters.HierarchySaliencyMapProjection.computeIncrementalNodeContours(
                weighted,
                radius=radius,
            ),
        ),
    )
    for context, entrypoint in radius_entrypoints:
        for radius in invalid_radii:
            require_radius_rejected(entrypoint, radius, context)

    mmcfilters.HierarchySaliencyMapValidation.validateHierarchyValuation(weighted, scores32, strict=True)
    formal_score_map32 = mmcfilters.HierarchySaliencyMap.computeSaliencyEdgeMap(weighted, scores32, strict=True)
    require(formal_score_map32["values"].dtype == np.float32, "formal saliency float32 dtype")
    require(np.array_equal(formal_score_map32["values"], np.array([30.0, 0.0], dtype=np.float32)), "formal saliency float32 values")

    topological_map = mmcfilters.HierarchySaliencyMap.computeTopologicalLevelEdgeMap(weighted)
    require(np.array_equal(topological_map["values"], np.array([1, 0], dtype=topological_map["values"].dtype)), "topological edge map values")
    require(np.array_equal(qfz_labels(topological_map, 0), np.array([0, 1, 2], dtype=np.int32)), "qfz lambda 0 partition")
    require(np.array_equal(qfz_labels(topological_map, 1), np.array([0, 1, 1], dtype=np.int32)), "qfz lambda 1 partition")
    require(np.array_equal(qfz_labels(topological_map, 2), np.array([0, 0, 0], dtype=np.int32)), "qfz lambda 2 partition")

    ranked_scores = mmcfilters.HierarchySaliencyMapValidation.rankHierarchyValuation(weighted, scores32, strict=True)
    require(np.array_equal(ranked_scores, np.array([0, 1, 2], dtype=ranked_scores.dtype)), "ranked strict hierarchy valuation")
    normalized_scores = mmcfilters.HierarchySaliencyMapValidation.computeNormalizedScores(weighted, scores32, strict=True)
    require(normalized_scores.dtype == np.float64, "normalized hierarchy valuation dtype")
    require(np.allclose(normalized_scores, np.array([0.0, 0.5, 1.0], dtype=np.float64)), "normalized strict hierarchy valuation")
    ranked_map = mmcfilters.HierarchySaliencyMap.computeSaliencyEdgeMap(
        weighted,
        ranked_scores,
        strict=True,
    )
    require(np.array_equal(ranked_map["values"], np.array([2, 0], dtype=ranked_map["values"].dtype)), "ranked valuation saliency values")
    canonical_ranked_map = mmcfilters.HierarchySaliencyMap.computeCanonicalRankedSaliencyEdgeMap(
        weighted,
        scores32,
        strict=True,
    )
    require(
        np.array_equal(canonical_ranked_map["values"], np.array([1, 0], dtype=canonical_ranked_map["values"].dtype)),
        "canonical rank uses effective edge levels",
    )
    appearance_levels = np.array([1, 1, 2], dtype=np.int32)
    appearance_map = mmcfilters.HierarchySaliencyMap.computeSaliencyEdgeMap(
        weighted,
        appearance_levels,
        strict=True,
        levelConvention=mmcfilters.HierarchyLevelConvention.PartitionAppearanceLevel,
    )
    require(
        np.array_equal(appearance_map["values"], np.array([1, 0], dtype=appearance_map["values"].dtype)),
        "partition appearance convention applies level(LCA)-1",
    )
    mmcfilters.HierarchySaliencyMapValidation.validateHierarchyConnectivity(weighted)
    mmcfilters.ComponentTreePartitionHierarchyAdapter.validate(weighted)
    adapter_levels = mmcfilters.ComponentTreePartitionHierarchyAdapter.computePartitionAppearanceLevels(weighted)
    require(np.array_equal(adapter_levels, appearance_levels), "component-tree adapter appearance levels")
    adapter_map = mmcfilters.ComponentTreePartitionHierarchyAdapter.computeSaliencyEdgeMap(
        weighted,
        adapter_levels,
        strict=True,
    )
    require(
        np.array_equal(adapter_map["values"], appearance_map["values"]),
        "component-tree adapter applies level(LCA)-1",
    )

    collapsing_scores = np.array([30.0, 30.0, 30.0], dtype=np.float32)
    mmcfilters.HierarchySaliencyMapValidation.validateHierarchyValuation(weighted, collapsing_scores)
    require_raises(
        lambda: mmcfilters.HierarchySaliencyMapValidation.validateHierarchyValuation(weighted, collapsing_scores, strict=True),
        "strict hierarchy valuation must reject equal parent-child levels",
    )
    negative_scores = np.array([-10.0, -5.0, 0.0], dtype=np.float32)
    mmcfilters.HierarchySaliencyMapValidation.validateHierarchyValuation(weighted, negative_scores)
    require_raises(
        lambda: mmcfilters.HierarchySaliencyMap.computeSaliencyEdgeMap(weighted, negative_scores),
        "formal saliency must reject negative valuation values",
    )
    require_raises(
        lambda: mmcfilters.HierarchySaliencyMapValidation.validateHierarchyValuation(weighted, negative_scores, nonnegative=True),
        "non-negative hierarchy valuation policy must reject negative values",
    )
    ranked_negative_scores = mmcfilters.HierarchySaliencyMapValidation.rankHierarchyValuation(weighted, negative_scores)
    require(np.array_equal(ranked_negative_scores, np.array([0, 1, 2], dtype=ranked_negative_scores.dtype)), "ranked negative valuation")
    normalized_negative_scores = mmcfilters.HierarchySaliencyMapValidation.computeNormalizedScores(weighted, negative_scores)
    require(np.allclose(normalized_negative_scores, np.array([0.0, 0.5, 1.0], dtype=np.float64)), "normalized negative valuation")
    require_raises(
        lambda: mmcfilters.HierarchySaliencyMapValidation.computeNormalizedScores(weighted, negative_scores, nonnegative=True),
        "non-negative normalized hierarchy valuation policy must reject negative values",
    )

    invalid_scores = np.array([10.0, 40.0, 30.0], dtype=np.float32)
    require_raises(
        lambda: mmcfilters.HierarchySaliencyMap.computeSaliencyEdgeMap(weighted, invalid_scores),
        "formal saliency must reject valuations that decrease toward the root",
    )
    require_raises(
        lambda: mmcfilters.HierarchySaliencyMapValidation.computeNormalizedScores(weighted, invalid_scores),
        "normalized saliency valuation must reject valuations that decrease toward the root",
    )
    non_finite_scores = np.array([10.0, 20.0, np.nan], dtype=np.float64)
    require_raises(
        lambda: mmcfilters.HierarchySaliencyMapValidation.validateHierarchyValuation(weighted, non_finite_scores),
        "formal saliency valuation must reject non-finite values",
    )

    normalized_max = mmcfilters.HierarchySaliencyMap.computeNormalizedAltitudeEdgeMap(weighted)
    require(normalized_max["values"].dtype == np.float64, "normalized max-tree edge map dtype")
    require(np.allclose(normalized_max["values"], np.array([1.0, 0.0], dtype=np.float64)), "normalized max-tree edge map values")

    max_pixels = mmcfilters.HierarchySaliencyMapProjection.edgeMapToPixelImage(
        normalized_max,
        mmcfilters.EdgeToPixelReducer.Max,
    )
    require(np.allclose(max_pixels, np.array([[1.0, 1.0, 0.0]], dtype=np.float64)), "edge-to-pixel max projection")
    mean_pixels = mmcfilters.HierarchySaliencyMapProjection.edgeMapToPixelImage(
        normalized_max,
        mmcfilters.EdgeToPixelReducer.Mean,
    )
    require(np.allclose(mean_pixels, np.array([[1.0, 0.5, 0.0]], dtype=np.float64)), "edge-to-pixel mean projection")

    high_cut = mmcfilters.HierarchySaliencyMapProjection.thresholdCut(normalized_max, 0.5)
    require(high_cut["numRows"] == 1, "high-threshold contour rows")
    require(high_cut["numCols"] == 3, "high-threshold contour cols")
    require(abs(high_cut["adjacencyRadius"] - 1.0) < 1e-12, "high-threshold contour radius")
    require(np.array_equal(high_cut["sources"], np.array([0], dtype=high_cut["sources"].dtype)), "high-threshold contour sources")
    require(np.array_equal(high_cut["targets"], np.array([1], dtype=high_cut["targets"].dtype)), "high-threshold contour targets")

    low_cut = mmcfilters.HierarchySaliencyMapProjection.thresholdCut(normalized_max, 0.199)
    require(np.array_equal(low_cut["sources"], np.array([0], dtype=low_cut["sources"].dtype)), "low-threshold contour sources")
    require(np.array_equal(low_cut["targets"], np.array([1], dtype=low_cut["targets"].dtype)), "low-threshold contour targets")

    node_contours = mmcfilters.HierarchySaliencyMapProjection.nodeContourEdges(weighted)
    require(np.array_equal(node_contours["sources"], np.array([0], dtype=node_contours["sources"].dtype)), "node contour sources")
    require(np.array_equal(node_contours["targets"], np.array([1], dtype=node_contours["targets"].dtype)), "node contour targets")
    require(np.array_equal(node_contours["nodes"], np.array([2], dtype=node_contours["nodes"].dtype)), "node contour LCA nodes")

    incremental_contours = mmcfilters.HierarchySaliencyMapProjection.computeIncrementalNodeContours(weighted)
    require(incremental_contours["numRows"] == 1, "incremental contour rows")
    require(incremental_contours["numCols"] == 3, "incremental contour cols")
    require(incremental_contours["numNodeSlots"] == weighted.numInternalNodeSlots, "incremental contour node slots")
    require(
        np.array_equal(incremental_contours["offsets"], np.array([0, 0, 0, 1], dtype=incremental_contours["offsets"].dtype)),
        "incremental contour offsets",
    )
    require(np.array_equal(incremental_contours["sources"], np.array([0], dtype=incremental_contours["sources"].dtype)), "incremental contour sources")
    require(np.array_equal(incremental_contours["targets"], np.array([1], dtype=incremental_contours["targets"].dtype)), "incremental contour targets")

    projected_scores = mmcfilters.HierarchySaliencyMapProjection.projectNodeValuation(incremental_contours, scores32)
    require(projected_scores["values"].dtype == np.float32, "incremental projected score dtype")
    require(np.array_equal(projected_scores["sources"], np.array([0], dtype=projected_scores["sources"].dtype)), "incremental projected score sources")
    require(np.array_equal(projected_scores["targets"], np.array([1], dtype=projected_scores["targets"].dtype)), "incremental projected score targets")
    require(np.array_equal(projected_scores["values"], np.array([30.0], dtype=np.float32)), "incremental projected score values")

    score_cut = mmcfilters.HierarchySaliencyMapProjection.thresholdByNodeValuation(incremental_contours, scores32, 25.0)
    require(np.array_equal(score_cut["sources"], np.array([0], dtype=score_cut["sources"].dtype)), "incremental score-cut sources")
    require(np.array_equal(score_cut["targets"], np.array([1], dtype=score_cut["targets"].dtype)), "incremental score-cut targets")

    min_weighted = make_three_pixel_min_tree(mmcfilters)
    normalized_min = mmcfilters.HierarchySaliencyMap.computeNormalizedAltitudeEdgeMap(min_weighted)
    require(np.allclose(normalized_min["values"], np.array([1.0, 0.0], dtype=np.float64)), "normalized min-tree edge map values")

    no_adjacency = make_three_pixel_tree(
        mmcfilters,
        mmcfilters.MorphologicalTreeKind.TREE_OF_SHAPES,
        radius=None,
    )
    require_raises(
        lambda: mmcfilters.HierarchySaliencyMap.computeTopologicalLevelEdgeMap(no_adjacency),
        "topological edge map must require stored or explicit adjacency",
    )
    tos_topological_valuation = np.array([0.0, 0.0, 1.0], dtype=np.float32)
    tos_formal = mmcfilters.HierarchySaliencyMap.computeSaliencyEdgeMap(
        no_adjacency,
        tos_topological_valuation,
        radius=1.0,
        strict=True,
    )
    require(np.array_equal(tos_formal["values"], np.array([1.0, 0.0], dtype=np.float32)), "tree-of-shapes formal saliency values")
    contour_fallback = mmcfilters.HierarchySaliencyMapProjection.nodeContourEdges(no_adjacency, radius=1.0)
    require(np.array_equal(contour_fallback["nodes"], np.array([2], dtype=contour_fallback["nodes"].dtype)), "explicit radius contour fallback nodes")
    incremental_fallback = mmcfilters.HierarchySaliencyMapProjection.computeIncrementalNodeContours(no_adjacency, radius=1.0)
    require(
        np.array_equal(incremental_fallback["offsets"], np.array([0, 0, 0, 1], dtype=incremental_fallback["offsets"].dtype)),
        "explicit radius incremental contour fallback offsets",
    )
    topological_fallback = mmcfilters.HierarchySaliencyMap.computeTopologicalLevelEdgeMap(no_adjacency, radius=1.0)
    require(np.array_equal(topological_fallback["values"], np.array([1, 0], dtype=topological_fallback["values"].dtype)), "explicit radius topological fallback values")
    require_raises(
        lambda: mmcfilters.HierarchySaliencyMap.computeNormalizedAltitudeEdgeMap(no_adjacency, radius=1.0),
        "normalized altitude edge map must reject trees without max/min polarity",
    )

    tos_image = np.array(
        [[1, 2, 1], [2, 3, 2], [1, 2, 1]],
        dtype=np.uint8,
    )
    self_dual_tos = mmcfilters.MorphologicalTreeFactory.createTreeOfShapes(
        tos_image,
        mmcfilters.ToSInterpolation.SelfDual,
    )
    self_dual_stored_map = mmcfilters.HierarchySaliencyMap.computeTopologicalLevelEdgeMap(
        self_dual_tos
    )
    require(
        self_dual_stored_map["values"].size > 0,
        "equal directional ToS stencils must define one stored saliency graph",
    )

    asymmetric_tos = mmcfilters.MorphologicalTreeFactory.createTreeOfShapes(
        tos_image,
        mmcfilters.ToSInterpolation.Min4cMax8c,
    )
    require_raises(
        lambda: mmcfilters.HierarchySaliencyMap.computeTopologicalLevelEdgeMap(asymmetric_tos),
        "distinct directional ToS stencils must require an explicit saliency graph",
    )
    asymmetric_explicit_map = mmcfilters.HierarchySaliencyMap.computeTopologicalLevelEdgeMap(
        asymmetric_tos,
        radius=1.0,
    )
    require(
        asymmetric_explicit_map["values"].size > 0,
        "an explicit radius must resolve directional ToS graph ambiguity",
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
