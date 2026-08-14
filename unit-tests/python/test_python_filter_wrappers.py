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
            "MorphologicalTreeFactory.create_max_tree",
            lambda radius: mmcfilters.MorphologicalTreeFactory.create_max_tree(image, radius=radius),
        ),
        (
            "MorphologicalTreeFactory.create_min_tree",
            lambda radius: mmcfilters.MorphologicalTreeFactory.create_min_tree(image, radius=radius),
        ),
    )
    for context, entrypoint in factory_radius_entrypoints:
        for radius in invalid_radii:
            require_radius_rejected(entrypoint, radius, context)

    valued_tree = mmcfilters.MorphologicalTreeFactory.create_max_tree(image)
    valued_tree_reconstruction = valued_tree.reconstruct_from_node_altitudes()

    adjacency = mmcfilters.RegularGridAdjacency2D(4, 4, 1.5)
    require(adjacency.size == 9, "RegularGridAdjacency2D stencil size")
    require(
        set(adjacency.neighbor_indices(1, 1))
        == {0, 1, 2, 4, 6, 8, 9, 10},
        "RegularGridAdjacency2D neighbor index traversal",
    )

    keep_all = [True] * valued_tree.num_internal_node_slots
    keep_all_mask = mmcfilters.NodePreservationMask(keep_all)
    direct_filter = mmcfilters.DirectAttributeFilter(valued_tree)
    hard_subtractive_filter = mmcfilters.HardSubtractiveAttributeFilter(valued_tree)
    soft_subtractive_filter = mmcfilters.SoftSubtractiveAttributeFilter(valued_tree)
    valued_tree_filters = mmcfilters.AttributeFilters(valued_tree)
    require(
        np.array_equal(direct_filter.apply_direct_attribute_filter(keep_all_mask), valued_tree_reconstruction),
        "valued_tree DirectAttributeFilter keep-all",
    )
    require(
        np.array_equal(
            hard_subtractive_filter.apply_hard_subtractive_attribute_filter(keep_all_mask),
            valued_tree_reconstruction.astype(np.int64),
        ),
        "valued_tree HardSubtractiveAttributeFilter keep-all",
    )
    require(
        np.array_equal(
            soft_subtractive_filter.apply_soft_subtractive_attribute_filter(
                np.ones(valued_tree.num_internal_node_slots, dtype=np.float32)
            ),
            valued_tree_reconstruction.astype(np.float32),
        ),
        "valued_tree SoftSubtractiveAttributeFilter unit scores",
    )
    reject_all_mask = mmcfilters.NodePreservationMask([False] * valued_tree.num_internal_node_slots)
    require(
        np.array_equal(
            hard_subtractive_filter.apply_hard_subtractive_attribute_filter(reject_all_mask),
            np.zeros_like(valued_tree_reconstruction, dtype=np.int64),
        ),
        "all-false hard subtractive mask must produce zero",
    )
    zero_scores = np.zeros(valued_tree.num_internal_node_slots, dtype=np.float64)
    require(
        np.array_equal(
            soft_subtractive_filter.apply_soft_subtractive_attribute_filter(zero_scores),
            np.zeros_like(valued_tree_reconstruction, dtype=np.float64),
        ),
        "all-zero soft subtractive scores must produce zero",
    )
    require(
        np.array_equal(valued_tree.reconstruct_from_node_contributions(zero_scores), np.zeros_like(valued_tree_reconstruction, dtype=np.float64)),
        "general zero-baseline node-contribution reconstruction",
    )
    prune_none_mask = mmcfilters.to_node_pruning_mask(keep_all_mask)
    require(prune_none_mask.to_list() == [False] * valued_tree.num_internal_node_slots, "explicit preservation-to-pruning conversion")
    require(mmcfilters.to_node_preservation_mask(prune_none_mask).to_list() == keep_all, "explicit pruning-to-preservation conversion")
    require(np.array_equal(valued_tree_filters.filtering_by_pruning_min(keep_all_mask), valued_tree_reconstruction), "valued_tree AttributeFilters pruning min keep-all mask")
    require(np.array_equal(valued_tree_filters.filtering_by_pruning_max(keep_all_mask), valued_tree_reconstruction), "valued_tree AttributeFilters pruning max keep-all mask")
    valued_tree_box_height = mmcfilters.Attribute.compute_single_attribute(valued_tree, mmcfilters.Attribute.BOUNDING_BOX_HEIGHT)
    valued_tree_box_height64 = mmcfilters.Attribute.compute_single_attribute(valued_tree, mmcfilters.Attribute.BOUNDING_BOX_HEIGHT, dtype=np.float64)
    require(valued_tree_box_height.dtype == np.float32, "default filter attribute dtype")
    require(valued_tree_box_height64.dtype == np.float64, "float64 filter attribute dtype")
    require(np.array_equal(valued_tree_filters.filtering_by_pruning_min(valued_tree_box_height, 1.0), valued_tree_reconstruction), "valued_tree AttributeFilters pruning min keep-all attribute")
    require(np.array_equal(valued_tree_filters.filtering_by_pruning_max(valued_tree_box_height, 1.0), valued_tree_reconstruction), "valued_tree AttributeFilters pruning max keep-all attribute")
    require(np.array_equal(valued_tree_filters.filtering_by_pruning_min(valued_tree_box_height64, 1.0), valued_tree_reconstruction), "valued_tree AttributeFilters pruning min float64 attribute")
    require(np.array_equal(valued_tree_filters.filtering_by_pruning_max(valued_tree_box_height64, 1.0), valued_tree_reconstruction), "valued_tree AttributeFilters pruning max float64 attribute")
    require_raises(lambda: valued_tree_filters.filtering_by_pruning_min(np.array([1.0], dtype=np.float32), 1.0), "valued_tree filtering_by_pruning_min must reject short attribute buffer")
    require_raises(
        lambda: direct_filter.apply_direct_attribute_filter(mmcfilters.NodePreservationMask([True])),
        "DirectAttributeFilter must reject a short preservation mask",
    )
    reject_root = keep_all.copy()
    reject_root[valued_tree.root] = False
    require_raises(
        lambda: direct_filter.apply_direct_attribute_filter(mmcfilters.NodePreservationMask(reject_root)),
        "DirectAttributeFilter must reject root suppression",
    )
    require_raises(
        lambda: soft_subtractive_filter.apply_soft_subtractive_attribute_filter(
            np.full(valued_tree.num_internal_node_slots, np.nan, dtype=np.float32)
        ),
        "SoftSubtractiveAttributeFilter must reject NaN scores",
    )
    require_raises(
        lambda: soft_subtractive_filter.apply_soft_subtractive_attribute_filter(
            np.full(valued_tree.num_internal_node_slots, -0.01, dtype=np.float32)
        ),
        "SoftSubtractiveAttributeFilter must reject scores below zero",
    )
    require_raises(
        lambda: soft_subtractive_filter.apply_soft_subtractive_attribute_filter(
            np.full(valued_tree.num_internal_node_slots, 1.01, dtype=np.float64)
        ),
        "SoftSubtractiveAttributeFilter must reject scores above one",
    )
    valued_tree_all_preserved = mmcfilters.NodePreservationMask(keep_all)
    altitude_adjusted = mmcfilters.adjust_node_preservation_mask_by_altitude_stability(
        valued_tree,
        valued_tree_all_preserved,
        2,
        mmcfilters.IncompleteStabilityWindowPolicy.PRESERVE_INPUT_DECISION,
    )
    depth_adjusted = mmcfilters.adjust_node_preservation_mask_by_depth_stability(
        valued_tree,
        valued_tree_all_preserved,
        2,
        mmcfilters.IncompleteStabilityWindowPolicy.PRESERVE_INPUT_DECISION,
    )
    require(altitude_adjusted.to_list() == keep_all, "altitude stability preserves an all-preserved input mask")
    require(depth_adjusted.to_list() == keep_all, "depth stability preserves an all-preserved input mask")
    threshold_mask = mmcfilters.compute_node_preservation_mask(valued_tree_box_height, 1.0)
    require(isinstance(threshold_mask, mmcfilters.NodePreservationMask), "thresholding returns NodePreservationMask")
    require(threshold_mask.to_list() == keep_all, "thresholding uses preservation polarity")

    viterbi_chain = mmcfilters.MorphologicalTreeFactory.create_from_higra_parent(
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
        np.array_equal(viterbi_chain_filters.filtering_by_viterbi_rule(chain_remove_attr, 1.0), np.array([[0]], dtype=np.uint8)),
        "Python AttributeFilters Viterbi chain remove",
    )
    require(
        np.array_equal(viterbi_chain_filters.filtering_by_viterbi_rule(chain_preserve_attr, 1.0), np.array([[5]], dtype=np.uint8)),
        "Python AttributeFilters Viterbi chain preserve float64",
    )
    require_raises(
        lambda: viterbi_chain_filters.filtering_by_viterbi_rule(np.array([1.0], dtype=np.float32), 1.0),
        "Python Viterbi must reject short attribute buffer",
    )
    require_raises(
        lambda: viterbi_chain_filters.filtering_by_viterbi_rule(np.array([1.0, np.nan, 1.0], dtype=np.float32), 1.0),
        "Python Viterbi must reject NaN attributes",
    )

    viterbi_branch = mmcfilters.MorphologicalTreeFactory.create_from_higra_parent(
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
            mmcfilters.AttributeFilters(viterbi_branch).filtering_by_viterbi_rule(branch_attr, 1.0),
            np.array([[5, 0, 0]], dtype=np.uint8),
        ),
        "Python AttributeFilters Viterbi branch",
    )

    valued_tree_gray_level_height = mmcfilters.Attribute.compute_single_attribute(valued_tree, mmcfilters.Attribute.GRAY_LEVEL_HEIGHT)
    valued_tree_gray_level_height64 = mmcfilters.Attribute.compute_single_attribute(valued_tree, mmcfilters.Attribute.GRAY_LEVEL_HEIGHT, dtype=np.float64)
    keep_all_policy = mmcfilters.ExtinctionSelectionPolicy.by_top_k(1024)
    keep_one_policy = mmcfilters.ExtinctionSelectionPolicy.by_top_k(1)
    rank_contours = mmcfilters.ExtinctionContourScorePolicy.RANK_SCORE
    require(np.array_equal(valued_tree_filters.filtering_by_extinction(valued_tree_gray_level_height, keep_all_policy), valued_tree_reconstruction), "valued_tree AttributeFilters extinction filtering keep-all")
    require(np.array_equal(valued_tree_filters.filtering_by_extinction(valued_tree_gray_level_height64, keep_all_policy), valued_tree_reconstruction), "valued_tree AttributeFilters extinction filtering float64")
    require(valued_tree_filters.contour_map_by_extinction(valued_tree_gray_level_height, keep_all_policy, rank_contours).shape == (4, 4), "valued_tree AttributeFilters extinction contour shape")
    saliency64 = valued_tree_filters.contour_map_by_extinction(valued_tree_gray_level_height64, keep_all_policy, rank_contours)
    require(saliency64.shape == (4, 4), "valued_tree AttributeFilters extinction float64 saliency shape")
    require(saliency64.dtype == np.float64, "valued_tree AttributeFilters extinction float64 saliency dtype")
    require_raises(lambda: valued_tree_filters.filtering_by_extinction(np.array([1.0], dtype=np.float32), keep_one_policy), "valued_tree filtering_by_extinction must reject short attribute buffer")
    require_raises(lambda: valued_tree_filters.contour_map_by_extinction(np.array([1.0], dtype=np.float32), keep_one_policy, rank_contours), "valued_tree contour_map_by_extinction must reject short attribute buffer")
    require_raises(lambda: valued_tree_filters.filtering_by_extinction(valued_tree_gray_level_height, mmcfilters.ExtinctionSelectionPolicy.by_top_k(-1)), "valued_tree filtering_by_extinction must reject negative keep count")
    require_raises(lambda: valued_tree_filters.contour_map_by_extinction(valued_tree_gray_level_height, mmcfilters.ExtinctionSelectionPolicy.by_top_k(-1), rank_contours), "valued_tree contour_map_by_extinction must reject negative keep count")
    require(
        not hasattr(mmcfilters, "AttributeOpeningPrimitivesFamily"),
        "AttributeOpeningPrimitivesFamily must not be exported in the Python API",
    )

    valued_tree_uao = mmcfilters.UltimateAttributeOpening(valued_tree, valued_tree_box_height)
    valued_tree_uao64 = mmcfilters.UltimateAttributeOpening(valued_tree, valued_tree_box_height64)
    valued_tree_uao.execute(int(valued_tree.num_rows))
    valued_tree_uao64.execute(float(valued_tree.num_rows))
    require(valued_tree_uao.get_max_contrast_image().shape == (4, 4), "valued_tree UltimateAttributeOpening max contrast shape")
    require(np.array_equal(valued_tree_uao64.get_max_contrast_image(), valued_tree_uao.get_max_contrast_image()), "valued_tree UltimateAttributeOpening float64 max contrast")
    require(valued_tree_uao.get_associated_image().shape == (4, 4), "valued_tree UltimateAttributeOpening associated image shape")
    require(np.array_equal(valued_tree_uao64.get_associated_image(), valued_tree_uao.get_associated_image()), "valued_tree UltimateAttributeOpening float64 associated image")
    require(valued_tree_uao.get_associated_colored_image().shape == (4, 12), "valued_tree UltimateAttributeOpening associated color image shape")
    valued_tree_uao.execute_with_mser(int(valued_tree.num_rows), 1)
    valued_tree_uao64.execute_with_mser(float(valued_tree.num_rows), 1)
    require(valued_tree_uao.get_max_contrast_image().shape == (4, 4), "valued_tree UltimateAttributeOpening MSER execute shape")
    require(valued_tree_uao64.get_max_contrast_image().shape == (4, 4), "valued_tree UltimateAttributeOpening float64 MSER execute shape")
    valued_tree_uao.execute_with_depth_stability(int(valued_tree.num_rows), 1)
    valued_tree_uao64.execute_with_depth_stability(float(valued_tree.num_rows), 1)
    require(valued_tree_uao.get_max_contrast_image().shape == (4, 4), "valued_tree UltimateAttributeOpening depth stability execute shape")
    require(valued_tree_uao64.get_max_contrast_image().shape == (4, 4), "valued_tree UltimateAttributeOpening float64 depth stability execute shape")

    tos = mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(image)
    tos_area = mmcfilters.Attribute.compute_single_topology_attribute(tos, mmcfilters.Attribute.AREA)
    tos_keep_all = [True] * tos.num_internal_node_slots
    tos_depth_adjusted = mmcfilters.adjust_node_preservation_mask_by_depth_stability(
        tos,
        mmcfilters.NodePreservationMask(tos_keep_all),
        1,
    )
    require(tos_depth_adjusted.to_list() == tos_keep_all, "ToS depth stability preserves an all-preserved input mask")
    require_raises(
        lambda: mmcfilters.adjust_node_preservation_mask_by_altitude_stability(
            tos,
            mmcfilters.NodePreservationMask(tos_keep_all),
            1,
        ),
        "ToS altitude stability must reject unconstrained altitude order",
    )
    depth_computer = mmcfilters.DepthStableRegionComputer(tos)
    require_raises(lambda: depth_computer.get_variation(0), "DepthStableRegionComputer variation requires computation")
    require_raises(lambda: depth_computer.get_variations(), "DepthStableRegionComputer variation buffer requires computation")
    require_raises(lambda: depth_computer.num_nodes(), "DepthStableRegionComputer selected count requires computation")
    require_raises(
        lambda: depth_computer.node_with_minimum_variation_in_window(0),
        "DepthStableRegionComputer minimum requires computation",
    )
    require_raises(
        lambda: depth_computer.ancestor_in_stability_window(0),
        "DepthStableRegionComputer ancestor requires computation",
    )
    require_raises(
        lambda: depth_computer.descendant_in_stability_window(0),
        "DepthStableRegionComputer descendant requires computation",
    )
    depth_mask = depth_computer.compute_by_depth(1)
    require(depth_mask.shape == (tos.num_internal_node_slots,), "DepthStableRegionComputer mask shape")
    require(depth_mask.dtype == np.uint8, "DepthStableRegionComputer mask dtype")
    depth_variations = depth_computer.get_variations()
    require(depth_variations.shape == (tos.num_internal_node_slots,), "DepthStableRegionComputer variation shape")
    require(depth_variations.dtype == np.float32, "DepthStableRegionComputer default variation dtype")
    for node_id in tos.alive_node_ids:
        require(
            np.isclose(
                depth_computer.get_variation(node_id),
                float(depth_variations[node_id]),
                atol=1e-6,
                equal_nan=True,
            ),
            "DepthStableRegionComputer get_variation",
        )
    require_raises(lambda: depth_computer.compute_by_depth(0), "DepthStableRegionComputer must reject a zero depth-window radius")
    depth_computer64 = mmcfilters.DepthStableRegionComputer(tos, tos_area.astype(np.float64))
    depth_computer64.compute_by_depth(1)
    require(depth_computer64.get_variations().dtype == np.float64, "DepthStableRegionComputer float64 variation dtype")
    require_raises(lambda: mmcfilters.ExtinctionValues(tos, tos_area), "ToS ExtinctionValues must be rejected")
    require_raises(lambda: tos_filters.filtering_by_extinction(tos_area, keep_one_policy), "ToS AttributeFilters extinction filtering must be rejected")
    tos_uao = mmcfilters.UltimateAttributeOpening(tos, tos_area)
    require_raises(lambda: tos_uao.execute_with_mser(int(tos.num_rows), 1), "ToS UltimateAttributeOpening must reject altitude MSER")
    tos_uao.execute_with_depth_stability(int(tos.num_rows), 1)
    require(tos_uao.get_max_contrast_image().shape == (4, 4), "ToS UltimateAttributeOpening depth stability execute shape")
    require_raises(
        lambda: mmcfilters.UltimateAttributeOpening(valued_tree, np.array([1.0], dtype=np.float32)),
        "valued_tree UltimateAttributeOpening must reject short attribute buffer",
    )
    nonsquare = np.array([[3, 3, 2], [1, 4, 5]], dtype=np.uint8)
    nonsquare_valued_tree = mmcfilters.MorphologicalTreeFactory.create_max_tree(nonsquare)
    nonsquare_box_height = mmcfilters.Attribute.compute_single_attribute(nonsquare_valued_tree, mmcfilters.Attribute.BOUNDING_BOX_HEIGHT)
    nonsquare_uao = mmcfilters.UltimateAttributeOpening(nonsquare_valued_tree, nonsquare_box_height)
    nonsquare_uao.execute(int(nonsquare_valued_tree.num_rows))
    require(nonsquare_uao.get_max_contrast_image().shape == (2, 3), "non-square UltimateAttributeOpening max contrast shape")
    require(nonsquare_uao.get_associated_image().shape == (2, 3), "non-square UltimateAttributeOpening associated image shape")
    require(nonsquare_uao.get_associated_colored_image().shape == (2, 9), "non-square UltimateAttributeOpening associated color image shape")

    valued_tree_extinction = mmcfilters.ExtinctionValues(valued_tree, valued_tree_gray_level_height)
    valued_tree_extinction64 = mmcfilters.ExtinctionValues(valued_tree, valued_tree_gray_level_height64)
    require(np.array_equal(valued_tree_extinction.filtering(keep_all_policy), valued_tree_reconstruction), "valued_tree ExtinctionValues filtering keep-all")
    require(np.array_equal(valued_tree_extinction64.filtering(keep_all_policy), valued_tree_reconstruction), "valued_tree ExtinctionValues float64 filtering keep-all")
    dominant_threshold = float(valued_tree_extinction.get_regional_extrema()[0][2])
    below_all_policy = mmcfilters.ExtinctionSelectionPolicy.by_threshold(-1.0)
    dominant_policy = mmcfilters.ExtinctionSelectionPolicy.by_threshold(dominant_threshold)
    require(
        np.array_equal(valued_tree_extinction.filtering(below_all_policy), valued_tree_reconstruction),
        "valued_tree ExtinctionValues threshold below all extinctions keeps all extrema",
    )
    require(
        np.array_equal(valued_tree_extinction.filtering(dominant_policy), valued_tree_extinction.filtering(keep_one_policy)),
        "valued_tree ExtinctionValues threshold at dominant extinction keeps strongest extremum",
    )
    require(
        np.array_equal(valued_tree_filters.filtering_by_extinction(valued_tree_gray_level_height, dominant_policy), valued_tree_extinction.filtering(keep_one_policy)),
        "valued_tree AttributeFilters threshold extinction filtering",
    )
    require(valued_tree_extinction.contour_map(keep_all_policy, rank_contours).shape == (4, 4), "valued_tree ExtinctionValues contour shape")
    require(valued_tree_extinction64.contour_map(keep_all_policy, rank_contours).dtype == np.float64, "valued_tree ExtinctionValues float64 contour dtype")
    extinction_attribute = valued_tree_extinction.get_extinction_value_attribute()
    require(extinction_attribute.dtype == np.float32, "valued_tree ExtinctionValues extinction attribute dtype")
    require(extinction_attribute.shape == (valued_tree.num_internal_node_slots,), "valued_tree ExtinctionValues extinction attribute shape")
    for leaf, _cutoff, value in valued_tree_extinction.get_regional_extrema():
        require(
            float(extinction_attribute[int(leaf)]) == float(np.float32(value)),
            "valued_tree ExtinctionValues leaf must receive its extinction value",
        )
    mmcfilters.HierarchySaliencyMapValidation.validate_hierarchy_valuation(valued_tree, extinction_attribute, nonnegative=True)
    ranked_extinction_valuation = valued_tree_extinction.compute_ranked_extinction_value_attribute()
    require(ranked_extinction_valuation.shape == (valued_tree.num_internal_node_slots,), "valued_tree ExtinctionValues ranked extinction attribute shape")
    mmcfilters.HierarchySaliencyMapValidation.validate_hierarchy_valuation(valued_tree, ranked_extinction_valuation, nonnegative=True)
    formal_extinction_map = valued_tree_extinction.compute_formal_saliency_edge_map(ranked=True)
    monotone_extinction_map = valued_tree_extinction.compute_monotone_extinction_projection(ranked=True)
    direct_extinction_map = mmcfilters.HierarchySaliencyMap.compute_canonical_ranked_saliency_edge_map(valued_tree, extinction_attribute)
    require(np.array_equal(monotone_extinction_map["sources"], direct_extinction_map["sources"]), "valued_tree monotone extinction projection sources")
    require(np.array_equal(monotone_extinction_map["targets"], direct_extinction_map["targets"]), "valued_tree monotone extinction projection targets")
    require(np.array_equal(monotone_extinction_map["values"], direct_extinction_map["values"]), "valued_tree monotone extinction projection values")
    require(np.all(formal_extinction_map["values"] >= 0), "valued_tree persistence saliency ranks are non-negative")
    formal_extinction64_map = valued_tree_extinction64.compute_formal_saliency_edge_map(radius=1.5)
    require(formal_extinction64_map["values"].dtype == np.float64, "valued_tree formal extinction saliency float64 dtype")
    require(abs(formal_extinction64_map["adjacency_radius"] - 1.5) < 1e-12, "valued_tree formal extinction explicit radius")
    for radius in invalid_radii:
        require_radius_rejected(
            lambda invalid_radius: valued_tree_extinction.compute_formal_saliency_edge_map(
                radius=invalid_radius,
                ranked=True,
            ),
            radius,
            "ExtinctionValues.compute_formal_saliency_edge_map",
        )
        require_radius_rejected(
            lambda invalid_radius: valued_tree_extinction.compute_monotone_extinction_projection(
                radius=invalid_radius,
                ranked=True,
            ),
            radius,
            "ExtinctionValues.compute_monotone_extinction_projection",
        )
    require_raises(lambda: valued_tree_extinction.filtering(mmcfilters.ExtinctionSelectionPolicy.by_top_k(-1)), "valued_tree ExtinctionValues filtering must reject negative keep count")
    require_raises(lambda: valued_tree_extinction.filtering(mmcfilters.ExtinctionSelectionPolicy.by_threshold(np.nan)), "valued_tree ExtinctionValues threshold filtering must reject NaN")
    require_raises(lambda: valued_tree_filters.filtering_by_extinction(valued_tree_gray_level_height, mmcfilters.ExtinctionSelectionPolicy.by_threshold(np.nan)), "valued_tree AttributeFilters threshold extinction filtering must reject NaN")
    require_raises(lambda: valued_tree_extinction.contour_map(mmcfilters.ExtinctionSelectionPolicy.by_top_k(-1), rank_contours), "valued_tree ExtinctionValues contour must reject negative keep count")
    require_raises(
        lambda: mmcfilters.ExtinctionValues(valued_tree, np.array([1.0], dtype=np.float32)),
        "valued_tree ExtinctionValues must reject short attribute buffer",
    )
    require_raises(
        lambda: valued_tree_filters.filtering_by_extinction(np.array([1.0], dtype=np.float32), dominant_policy),
        "valued_tree AttributeFilters threshold extinction filtering must reject short attribute buffer",
    )

    stale_valued_tree = mmcfilters.MorphologicalTreeFactory.create_max_tree(image)
    stale_keep_all = [True] * stale_valued_tree.num_internal_node_slots
    stale_attribute = mmcfilters.Attribute.compute_single_attribute(stale_valued_tree, mmcfilters.Attribute.GRAY_LEVEL_HEIGHT)
    stale_filter = mmcfilters.DirectAttributeFilter(stale_valued_tree)
    stale_extinction = mmcfilters.ExtinctionValues(stale_valued_tree, stale_attribute)
    stale_uao = mmcfilters.UltimateAttributeOpening(stale_valued_tree, stale_attribute)
    stale_valued_tree.merge_node_into_parent(4)
    require_raises(
        lambda: stale_filter.apply_direct_attribute_filter(mmcfilters.NodePreservationMask(stale_keep_all)),
        "DirectAttributeFilter must reject use after topology mutation",
    )
    require_raises(lambda: stale_extinction.filtering(keep_all_policy), "ExtinctionValues must reject filtering after topology mutation")
    require_raises(lambda: stale_extinction.contour_map(keep_all_policy, rank_contours), "ExtinctionValues must reject contour after topology mutation")
    require_raises(lambda: stale_extinction.get_extinction_value_attribute(), "ExtinctionValues must reject extinction attribute after topology mutation")
    require_raises(lambda: stale_extinction.compute_formal_saliency_edge_map(), "ExtinctionValues must reject formal saliency after topology mutation")
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
    min_parent, min_altitude = casf.export_min_tree()
    max_parent, max_altitude = casf.export_max_tree()
    require(len(min_parent) == len(min_altitude), "CASF exported min-tree shape")
    require(len(max_parent) == len(max_altitude), "CASF exported max-tree shape")
    require(casf.min_tree.reconstruct_from_node_altitudes().shape == image.shape, "CASF min_tree property")
    require(casf.max_tree.reconstruct_from_node_altitudes().shape == image.shape, "CASF max_tree property")

    bbox_casf = mmcfilters.CasfComponentTrees(image, mmcfilters.CasfComponentTreesAttribute.BOUNDING_BOX_DIAGONAL)
    require(bbox_casf.filter([2.0]).shape == image.shape, "CASF bounding-box path")

    min_for_adjust = mmcfilters.MorphologicalTreeFactory.create_min_tree(image)
    max_for_adjust = mmcfilters.MorphologicalTreeFactory.create_max_tree(image)
    adjust = mmcfilters.DualMinMaxTreeIncrementalFilter(min_for_adjust, max_for_adjust)
    max_candidates = [
        node_id
        for node_id in max_for_adjust.alive_node_ids
        if node_id != max_for_adjust.root and max_for_adjust.proper_part_cardinality(node_id) <= 1
    ]
    adjust.prune_max_tree_and_update_min_tree(max_candidates[:1])
    require(adjust.min_tree.reconstruct_from_node_altitudes().shape == image.shape, "adjust min_tree property after max prune")
    require(adjust.max_tree.reconstruct_from_node_altitudes().shape == image.shape, "adjust max_tree property after max prune")

    min_candidates = [
        node_id
        for node_id in adjust.min_tree.alive_node_ids
        if node_id != adjust.min_tree.root and adjust.min_tree.proper_part_cardinality(node_id) <= 1
    ]
    adjust.prune_min_tree_and_update_max_tree(min_candidates[:1])
    require(adjust.min_tree.reconstruct_from_node_altitudes().shape == image.shape, "adjust min_tree property after min prune")
    require(adjust.max_tree.reconstruct_from_node_altitudes().shape == image.shape, "adjust max_tree property after min prune")

    print("python filter wrappers ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
