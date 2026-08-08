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


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: test_python_filter_wrappers.py <build-dir>")

    build_dir = pathlib.Path(sys.argv[1]).resolve()
    mmcfilters = load_native_module(build_dir)
    require(hasattr(mmcfilters, "__version__"), "package import must expose __version__")

    image = np.array(
        [
            [3, 3, 2, 2],
            [3, 4, 4, 2],
            [1, 4, 5, 2],
            [1, 1, 5, 0],
        ],
        dtype=np.uint8,
    )

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
    factory_radius_entrypoints = (
        (
            "MorphologicalTreeFactory.createMaxTree",
            lambda radius: mmcfilters.MorphologicalTreeFactory.createMaxTree(image, radius=radius),
        ),
        (
            "MorphologicalTreeFactory.createMinTree",
            lambda radius: mmcfilters.MorphologicalTreeFactory.createMinTree(image, radius=radius),
        ),
    )
    for context, entrypoint in factory_radius_entrypoints:
        for radius in invalid_radii:
            require_radius_rejected(entrypoint, radius, context)

    weighted = mmcfilters.MorphologicalTreeFactory.createMaxTree(image)
    weighted_reconstruction = weighted.reconstructionImage()

    adjacency = mmcfilters.RegularGridAdjacency2D(4, 4, 1.5)
    require(adjacency.size == 9, "RegularGridAdjacency2D stencil size")
    require(
        set(adjacency.neighborIndices(1, 1))
        == {0, 1, 2, 4, 6, 8, 9, 10},
        "RegularGridAdjacency2D neighbor index traversal",
    )

    keep_all = [True] * weighted.numInternalNodeSlots
    weighted_filters = mmcfilters.AttributeFilters(weighted)
    require(np.array_equal(weighted_filters.filteringDirectRule(keep_all), weighted_reconstruction), "weighted AttributeFilters direct rule keep-all")
    require(np.array_equal(weighted_filters.filteringSubtractiveRule(keep_all), weighted_reconstruction), "weighted AttributeFilters subtractive rule keep-all")
    require(np.array_equal(weighted_filters.filteringByPruningMin(keep_all), weighted_reconstruction), "weighted AttributeFilters pruning min keep-all criterion")
    require(np.array_equal(weighted_filters.filteringByPruningMax(keep_all), weighted_reconstruction), "weighted AttributeFilters pruning max keep-all criterion")
    weighted_box_height = mmcfilters.Attribute.computeSingleAttribute(weighted, mmcfilters.Attribute.BOX_HEIGHT)
    weighted_box_height64 = mmcfilters.Attribute.computeSingleAttribute(weighted, mmcfilters.Attribute.BOX_HEIGHT, dtype=np.float64)
    require(weighted_box_height.dtype == np.float32, "default filter attribute dtype")
    require(weighted_box_height64.dtype == np.float64, "float64 filter attribute dtype")
    require(np.array_equal(weighted_filters.filteringByPruningMin(weighted_box_height, 1.0), weighted_reconstruction), "weighted AttributeFilters pruning min keep-all attribute")
    require(np.array_equal(weighted_filters.filteringByPruningMax(weighted_box_height, 1.0), weighted_reconstruction), "weighted AttributeFilters pruning max keep-all attribute")
    require(np.array_equal(weighted_filters.filteringByPruningMin(weighted_box_height64, 1.0), weighted_reconstruction), "weighted AttributeFilters pruning min float64 attribute")
    require(np.array_equal(weighted_filters.filteringByPruningMax(weighted_box_height64, 1.0), weighted_reconstruction), "weighted AttributeFilters pruning max float64 attribute")
    require_raises(lambda: weighted_filters.filteringByPruningMin(np.array([1.0], dtype=np.float32), 1.0), "weighted filteringByPruningMin must reject short attribute buffer")
    require_raises(lambda: weighted_filters.filteringDirectRule([True]), "weighted filteringDirectRule must reject short criterion")
    require(weighted_filters.getAdaptiveCriterion(keep_all, 2) == [False] * weighted.numInternalNodeSlots, "weighted AttributeFilters adaptive criterion on all-true input")
    require(weighted_filters.getAdaptiveCriterionByDepth(keep_all, 2) == [False] * weighted.numInternalNodeSlots, "weighted AttributeFilters depth adaptive criterion on all-true input")

    viterbi_chain = mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
        [1, 2, 3, 3],
        np.array([5, 5, 3, 0], dtype=np.uint8),
        1,
        1,
        mmcfilters.MorphologicalTreeKind.MAX_TREE,
        1.5,
    )
    viterbi_chain_filters = mmcfilters.AttributeFilters(viterbi_chain)
    chain_remove_attr = np.array([1.1, 0.0, 2.0], dtype=np.float32)
    chain_preserve_attr = np.array([10.0, 0.9, 2.0], dtype=np.float64)
    require(
        np.array_equal(viterbi_chain_filters.filteringByViterbiRule(chain_remove_attr, 1.0), np.array([[0]], dtype=np.uint8)),
        "Python AttributeFilters Viterbi chain remove",
    )
    require(
        np.array_equal(viterbi_chain_filters.filteringByViterbiRule(chain_preserve_attr, 1.0), np.array([[5]], dtype=np.uint8)),
        "Python AttributeFilters Viterbi chain preserve float64",
    )
    require_raises(
        lambda: viterbi_chain_filters.filteringByViterbiRule(np.array([1.0], dtype=np.float32), 1.0),
        "Python Viterbi must reject short attribute buffer",
    )
    require_raises(
        lambda: viterbi_chain_filters.filteringByViterbiRule(np.array([1.0, np.nan, 1.0], dtype=np.float32), 1.0),
        "Python Viterbi must reject NaN attributes",
    )

    viterbi_branch = mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
        [3, 4, 4, 5, 5, 5],
        np.array([5, 4, 4, 5, 4, 0], dtype=np.uint8),
        1,
        3,
        mmcfilters.MorphologicalTreeKind.MAX_TREE,
        1.5,
    )
    branch_attr = np.array([2.0, 0.0, 2.0], dtype=np.float32)
    require(
        np.array_equal(
            mmcfilters.AttributeFilters(viterbi_branch).filteringByViterbiRule(branch_attr, 1.0),
            np.array([[5, 0, 0]], dtype=np.uint8),
        ),
        "Python AttributeFilters Viterbi branch",
    )

    weighted_level_attr = mmcfilters.Attribute.computeSingleAttribute(weighted, mmcfilters.Attribute.LEVEL)
    weighted_level_attr64 = mmcfilters.Attribute.computeSingleAttribute(weighted, mmcfilters.Attribute.LEVEL, dtype=np.float64)
    keep_all_policy = mmcfilters.ExtinctionSelectionPolicy.byTopK(1024)
    keep_one_policy = mmcfilters.ExtinctionSelectionPolicy.byTopK(1)
    rank_contours = mmcfilters.ExtinctionContourScorePolicy.RankScore
    require(np.array_equal(weighted_filters.filteringByExtinction(weighted_level_attr, keep_all_policy), weighted_reconstruction), "weighted AttributeFilters extinction filtering keep-all")
    require(np.array_equal(weighted_filters.filteringByExtinction(weighted_level_attr64, keep_all_policy), weighted_reconstruction), "weighted AttributeFilters extinction filtering float64")
    require(weighted_filters.contourMapByExtinction(weighted_level_attr, keep_all_policy, rank_contours).shape == (4, 4), "weighted AttributeFilters extinction contour shape")
    saliency64 = weighted_filters.contourMapByExtinction(weighted_level_attr64, keep_all_policy, rank_contours)
    require(saliency64.shape == (4, 4), "weighted AttributeFilters extinction float64 saliency shape")
    require(saliency64.dtype == np.float64, "weighted AttributeFilters extinction float64 saliency dtype")
    require_raises(lambda: weighted_filters.filteringByExtinction(np.array([1.0], dtype=np.float32), keep_one_policy), "weighted filteringByExtinction must reject short attribute buffer")
    require_raises(lambda: weighted_filters.contourMapByExtinction(np.array([1.0], dtype=np.float32), keep_one_policy, rank_contours), "weighted contourMapByExtinction must reject short attribute buffer")
    require_raises(lambda: weighted_filters.filteringByExtinction(weighted_level_attr, mmcfilters.ExtinctionSelectionPolicy.byTopK(-1)), "weighted filteringByExtinction must reject negative keep count")
    require_raises(lambda: weighted_filters.contourMapByExtinction(weighted_level_attr, mmcfilters.ExtinctionSelectionPolicy.byTopK(-1), rank_contours), "weighted contourMapByExtinction must reject negative keep count")
    require(
        not hasattr(mmcfilters, "AttributeOpeningPrimitivesFamily"),
        "AttributeOpeningPrimitivesFamily must not be exported in the Python API",
    )

    weighted_uao = mmcfilters.UltimateAttributeOpening(weighted, weighted_box_height)
    weighted_uao64 = mmcfilters.UltimateAttributeOpening(weighted, weighted_box_height64)
    weighted_uao.execute(int(weighted.numRows))
    weighted_uao64.execute(float(weighted.numRows))
    require(weighted_uao.getMaxContrastImage().shape == (4, 4), "weighted UltimateAttributeOpening max contrast shape")
    require(np.array_equal(weighted_uao64.getMaxContrastImage(), weighted_uao.getMaxContrastImage()), "weighted UltimateAttributeOpening float64 max contrast")
    require(weighted_uao.getAssociatedImage().shape == (4, 4), "weighted UltimateAttributeOpening associated image shape")
    require(np.array_equal(weighted_uao64.getAssociatedImage(), weighted_uao.getAssociatedImage()), "weighted UltimateAttributeOpening float64 associated image")
    require(weighted_uao.getAssociatedColoredImage().shape == (4, 12), "weighted UltimateAttributeOpening associated color image shape")
    weighted_uao.executeWithMSER(int(weighted.numRows), 1)
    weighted_uao64.executeWithMSER(float(weighted.numRows), 1)
    require(weighted_uao.getMaxContrastImage().shape == (4, 4), "weighted UltimateAttributeOpening MSER execute shape")
    require(weighted_uao64.getMaxContrastImage().shape == (4, 4), "weighted UltimateAttributeOpening float64 MSER execute shape")
    weighted_uao.executeWithDepthStability(int(weighted.numRows), 1)
    weighted_uao64.executeWithDepthStability(float(weighted.numRows), 1)
    require(weighted_uao.getMaxContrastImage().shape == (4, 4), "weighted UltimateAttributeOpening depth stability execute shape")
    require(weighted_uao64.getMaxContrastImage().shape == (4, 4), "weighted UltimateAttributeOpening float64 depth stability execute shape")

    tos = mmcfilters.MorphologicalTreeFactory.createTreeOfShapes(image)
    tos_area = mmcfilters.Attribute.computeSingleTopologyAttribute(tos, mmcfilters.Attribute.AREA)
    tos_filters = mmcfilters.AttributeFilters(tos)
    tos_keep_all = [True] * tos.numInternalNodeSlots
    require(tos_filters.getAdaptiveCriterionByDepth(tos_keep_all, 1) == [False] * tos.numInternalNodeSlots, "ToS depth adaptive criterion on all-true input")
    depth_computer = mmcfilters.DepthStableRegionComputer(tos)
    require_raises(lambda: depth_computer.getVariation(0), "DepthStableRegionComputer variation requires computation")
    require_raises(lambda: depth_computer.getVariations(), "DepthStableRegionComputer variation buffer requires computation")
    require_raises(lambda: depth_computer.getNumNodes(), "DepthStableRegionComputer selected count requires computation")
    require_raises(
        lambda: depth_computer.nodeWithMinimumVariationInWindow(0),
        "DepthStableRegionComputer minimum requires computation",
    )
    require_raises(
        lambda: depth_computer.ascendantInStabilityWindow(0),
        "DepthStableRegionComputer ascendant requires computation",
    )
    require_raises(
        lambda: depth_computer.descendantInStabilityWindow(0),
        "DepthStableRegionComputer descendant requires computation",
    )
    depth_mask = depth_computer.computeByDepth(1)
    require(depth_mask.shape == (tos.numInternalNodeSlots,), "DepthStableRegionComputer mask shape")
    require(depth_mask.dtype == np.uint8, "DepthStableRegionComputer mask dtype")
    depth_variations = depth_computer.getVariations()
    require(depth_variations.shape == (tos.numInternalNodeSlots,), "DepthStableRegionComputer variation shape")
    require(depth_variations.dtype == np.float32, "DepthStableRegionComputer default variation dtype")
    for node_id in tos.aliveNodeIds:
        require(
            np.isclose(
                depth_computer.getVariation(node_id),
                float(depth_variations[node_id]),
                atol=1e-6,
                equal_nan=True,
            ),
            "DepthStableRegionComputer getVariation",
        )
    require_raises(lambda: depth_computer.computeByDepth(0), "DepthStableRegionComputer must reject zero depth delta")
    depth_computer64 = mmcfilters.DepthStableRegionComputer(tos, tos_area.astype(np.float64))
    depth_computer64.computeByDepth(1)
    require(depth_computer64.getVariations().dtype == np.float64, "DepthStableRegionComputer float64 variation dtype")
    require_raises(lambda: mmcfilters.ExtinctionValues(tos, tos_area), "ToS ExtinctionValues must be rejected")
    require_raises(lambda: tos_filters.filteringByExtinction(tos_area, keep_one_policy), "ToS AttributeFilters extinction filtering must be rejected")
    tos_uao = mmcfilters.UltimateAttributeOpening(tos, tos_area)
    require_raises(lambda: tos_uao.executeWithMSER(int(tos.numRows), 1), "ToS UltimateAttributeOpening must reject altitude MSER")
    tos_uao.executeWithDepthStability(int(tos.numRows), 1)
    require(tos_uao.getMaxContrastImage().shape == (4, 4), "ToS UltimateAttributeOpening depth stability execute shape")
    require_raises(
        lambda: mmcfilters.UltimateAttributeOpening(weighted, np.array([1.0], dtype=np.float32)),
        "weighted UltimateAttributeOpening must reject short attribute buffer",
    )
    nonsquare = np.array([[3, 3, 2], [1, 4, 5]], dtype=np.uint8)
    nonsquare_weighted = mmcfilters.MorphologicalTreeFactory.createMaxTree(nonsquare)
    nonsquare_box_height = mmcfilters.Attribute.computeSingleAttribute(nonsquare_weighted, mmcfilters.Attribute.BOX_HEIGHT)
    nonsquare_uao = mmcfilters.UltimateAttributeOpening(nonsquare_weighted, nonsquare_box_height)
    nonsquare_uao.execute(int(nonsquare_weighted.numRows))
    require(nonsquare_uao.getMaxContrastImage().shape == (2, 3), "non-square UltimateAttributeOpening max contrast shape")
    require(nonsquare_uao.getAssociatedImage().shape == (2, 3), "non-square UltimateAttributeOpening associated image shape")
    require(nonsquare_uao.getAssociatedColoredImage().shape == (2, 9), "non-square UltimateAttributeOpening associated color image shape")

    weighted_extinction = mmcfilters.ExtinctionValues(weighted, weighted_level_attr)
    weighted_extinction64 = mmcfilters.ExtinctionValues(weighted, weighted_level_attr64)
    require(np.array_equal(weighted_extinction.filtering(keep_all_policy), weighted_reconstruction), "weighted ExtinctionValues filtering keep-all")
    require(np.array_equal(weighted_extinction64.filtering(keep_all_policy), weighted_reconstruction), "weighted ExtinctionValues float64 filtering keep-all")
    dominant_threshold = float(weighted_extinction.getRegionalExtrema()[0][2])
    below_all_policy = mmcfilters.ExtinctionSelectionPolicy.byThreshold(-1.0)
    dominant_policy = mmcfilters.ExtinctionSelectionPolicy.byThreshold(dominant_threshold)
    require(
        np.array_equal(weighted_extinction.filtering(below_all_policy), weighted_reconstruction),
        "weighted ExtinctionValues threshold below all extinctions keeps all extrema",
    )
    require(
        np.array_equal(weighted_extinction.filtering(dominant_policy), weighted_extinction.filtering(keep_one_policy)),
        "weighted ExtinctionValues threshold at dominant extinction keeps strongest extremum",
    )
    require(
        np.array_equal(weighted_filters.filteringByExtinction(weighted_level_attr, dominant_policy), weighted_extinction.filtering(keep_one_policy)),
        "weighted AttributeFilters threshold extinction filtering",
    )
    require(weighted_extinction.contourMap(keep_all_policy, rank_contours).shape == (4, 4), "weighted ExtinctionValues contour shape")
    require(weighted_extinction64.contourMap(keep_all_policy, rank_contours).dtype == np.float64, "weighted ExtinctionValues float64 contour dtype")
    extinction_attribute = weighted_extinction.getExtinctionValueAttribute()
    require(extinction_attribute.dtype == np.float32, "weighted ExtinctionValues extinction attribute dtype")
    require(extinction_attribute.shape == (weighted.numInternalNodeSlots,), "weighted ExtinctionValues extinction attribute shape")
    for leaf, _cutoff, value in weighted_extinction.getRegionalExtrema():
        require(
            float(extinction_attribute[int(leaf)]) == float(np.float32(value)),
            "weighted ExtinctionValues leaf must receive its extinction value",
        )
    mmcfilters.HierarchySaliencyMapValidation.validateHierarchyValuation(weighted, extinction_attribute, nonnegative=True)
    ranked_extinction_valuation = weighted_extinction.computeRankedExtinctionValueAttribute()
    require(ranked_extinction_valuation.shape == (weighted.numInternalNodeSlots,), "weighted ExtinctionValues ranked extinction attribute shape")
    mmcfilters.HierarchySaliencyMapValidation.validateHierarchyValuation(weighted, ranked_extinction_valuation, nonnegative=True)
    formal_extinction_map = weighted_extinction.computeFormalSaliencyEdgeMap(ranked=True)
    monotone_extinction_map = weighted_extinction.computeMonotoneExtinctionProjection(ranked=True)
    direct_extinction_map = mmcfilters.HierarchySaliencyMap.computeCanonicalRankedSaliencyEdgeMap(weighted, extinction_attribute)
    require(np.array_equal(monotone_extinction_map["sources"], direct_extinction_map["sources"]), "weighted monotone extinction projection sources")
    require(np.array_equal(monotone_extinction_map["targets"], direct_extinction_map["targets"]), "weighted monotone extinction projection targets")
    require(np.array_equal(monotone_extinction_map["values"], direct_extinction_map["values"]), "weighted monotone extinction projection values")
    require(np.all(formal_extinction_map["values"] >= 0), "weighted persistence saliency ranks are non-negative")
    formal_extinction64_map = weighted_extinction64.computeFormalSaliencyEdgeMap(radius=1.5)
    require(formal_extinction64_map["values"].dtype == np.float64, "weighted formal extinction saliency float64 dtype")
    require(abs(formal_extinction64_map["adjacencyRadius"] - 1.5) < 1e-12, "weighted formal extinction explicit radius")
    for radius in invalid_radii:
        require_radius_rejected(
            lambda invalid_radius: weighted_extinction.computeFormalSaliencyEdgeMap(
                radius=invalid_radius,
                ranked=True,
            ),
            radius,
            "ExtinctionValues.computeFormalSaliencyEdgeMap",
        )
        require_radius_rejected(
            lambda invalid_radius: weighted_extinction.computeMonotoneExtinctionProjection(
                radius=invalid_radius,
                ranked=True,
            ),
            radius,
            "ExtinctionValues.computeMonotoneExtinctionProjection",
        )
    require_raises(lambda: weighted_extinction.filtering(mmcfilters.ExtinctionSelectionPolicy.byTopK(-1)), "weighted ExtinctionValues filtering must reject negative keep count")
    require_raises(lambda: weighted_extinction.filtering(mmcfilters.ExtinctionSelectionPolicy.byThreshold(np.nan)), "weighted ExtinctionValues threshold filtering must reject NaN")
    require_raises(lambda: weighted_filters.filteringByExtinction(weighted_level_attr, mmcfilters.ExtinctionSelectionPolicy.byThreshold(np.nan)), "weighted AttributeFilters threshold extinction filtering must reject NaN")
    require_raises(lambda: weighted_extinction.contourMap(mmcfilters.ExtinctionSelectionPolicy.byTopK(-1), rank_contours), "weighted ExtinctionValues contour must reject negative keep count")
    require_raises(
        lambda: mmcfilters.ExtinctionValues(weighted, np.array([1.0], dtype=np.float32)),
        "weighted ExtinctionValues must reject short attribute buffer",
    )
    require_raises(
        lambda: weighted_filters.filteringByExtinction(np.array([1.0], dtype=np.float32), dominant_policy),
        "weighted AttributeFilters threshold extinction filtering must reject short attribute buffer",
    )

    stale_weighted = mmcfilters.MorphologicalTreeFactory.createMaxTree(image)
    stale_keep_all = [True] * stale_weighted.numInternalNodeSlots
    stale_level_attr = mmcfilters.Attribute.computeSingleAttribute(stale_weighted, mmcfilters.Attribute.LEVEL)
    stale_filters = mmcfilters.AttributeFilters(stale_weighted)
    stale_extinction = mmcfilters.ExtinctionValues(stale_weighted, stale_level_attr)
    stale_uao = mmcfilters.UltimateAttributeOpening(stale_weighted, stale_level_attr)
    stale_weighted.mergeNodeIntoParent(4)
    require_raises(lambda: stale_filters.filteringDirectRule(stale_keep_all), "AttributeFilters must reject use after topology mutation")
    require_raises(lambda: stale_extinction.filtering(keep_all_policy), "ExtinctionValues must reject filtering after topology mutation")
    require_raises(lambda: stale_extinction.contourMap(keep_all_policy, rank_contours), "ExtinctionValues must reject contour after topology mutation")
    require_raises(lambda: stale_extinction.getExtinctionValueAttribute(), "ExtinctionValues must reject extinction attribute after topology mutation")
    require_raises(lambda: stale_extinction.computeFormalSaliencyEdgeMap(), "ExtinctionValues must reject formal saliency after topology mutation")
    require_raises(lambda: stale_uao.execute(4), "UltimateAttributeOpening must reject execute after topology mutation")

    casf = mmcfilters.CasfComponentTrees(image, mmcfilters.CasfComponentTreesAttribute.AREA)
    require_raises(
        lambda: mmcfilters.CasfComponentTrees(image.astype(np.int32), mmcfilters.CasfComponentTreesAttribute.AREA),
        "CASF must reject non-uint8 integer arrays",
    )
    require_raises(
        lambda: mmcfilters.CasfComponentTrees(image.astype(np.int64), mmcfilters.CasfComponentTreesAttribute.AREA),
        "CASF must reject int64 arrays",
    )
    require_raises(
        lambda: mmcfilters.CasfComponentTrees(image.astype(np.float32), mmcfilters.CasfComponentTreesAttribute.AREA),
        "CASF must reject float arrays",
    )
    require_raises(
        lambda: mmcfilters.CasfComponentTrees(image.astype(np.float64), mmcfilters.CasfComponentTreesAttribute.AREA),
        "CASF must reject float64 arrays",
    )
    require_raises(
        lambda: mmcfilters.CasfComponentTrees(image.astype(bool), mmcfilters.CasfComponentTreesAttribute.AREA),
        "CASF must reject bool arrays",
    )
    require_raises(
        lambda: mmcfilters.CasfComponentTrees(image.astype(object), mmcfilters.CasfComponentTreesAttribute.AREA),
        "CASF must reject object arrays",
    )
    require_raises(
        lambda: mmcfilters.CasfComponentTrees(image[:, ::-1], mmcfilters.CasfComponentTreesAttribute.AREA),
        "CASF must reject non-contiguous uint8 arrays",
    )
    require(np.array_equal(casf.filter([]), image), "CASF empty threshold sequence must preserve image")
    require(np.array_equal(casf.filter([0.0]), image), "CASF zero threshold must preserve image")
    filtered = casf.filter([2.0])
    require(filtered.shape == image.shape, "CASF positive threshold image shape")
    min_parent, min_altitude = casf.exportMinTree()
    max_parent, max_altitude = casf.exportMaxTree()
    require(len(min_parent) == len(min_altitude), "CASF exported min-tree shape")
    require(len(max_parent) == len(max_altitude), "CASF exported max-tree shape")
    require(casf.minTree.reconstructionImage().shape == image.shape, "CASF minTree property")
    require(casf.maxTree.reconstructionImage().shape == image.shape, "CASF maxTree property")

    bbox_casf = mmcfilters.CasfComponentTrees(image, mmcfilters.CasfComponentTreesAttribute.BOUNDING_BOX_DIAGONAL)
    require(bbox_casf.filter([2.0]).shape == image.shape, "CASF bounding-box path")

    min_for_adjust = mmcfilters.MorphologicalTreeFactory.createMinTree(image)
    max_for_adjust = mmcfilters.MorphologicalTreeFactory.createMaxTree(image)
    adjust = mmcfilters.DualMinMaxTreeIncrementalFilter(min_for_adjust, max_for_adjust)
    max_candidates = [
        node_id
        for node_id in max_for_adjust.getAliveNodeIds()
        if node_id != max_for_adjust.getRoot() and max_for_adjust.getNumProperParts(node_id) <= 1
    ]
    adjust.pruneMaxTreeAndUpdateMinTree(max_candidates[:1])
    require(adjust.minTree.reconstructionImage().shape == image.shape, "adjust minTree property after max prune")
    require(adjust.maxTree.reconstructionImage().shape == image.shape, "adjust maxTree property after max prune")

    min_candidates = [
        node_id
        for node_id in adjust.minTree.getAliveNodeIds()
        if node_id != adjust.minTree.getRoot() and adjust.minTree.getNumProperParts(node_id) <= 1
    ]
    adjust.pruneMinTreeAndUpdateMaxTree(min_candidates[:1])
    require(adjust.minTree.reconstructionImage().shape == image.shape, "adjust minTree property after min prune")
    require(adjust.maxTree.reconstructionImage().shape == image.shape, "adjust maxTree property after min prune")

    print("python filter wrappers ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
