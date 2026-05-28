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

    weighted = mmcfilters.MorphologicalTreeFactory.createMaxTree(image)
    weighted_reconstruction = weighted.reconstructionImage()

    adjacency = mmcfilters.AdjacencyRelation(4, 4, 1.5)
    require(adjacency.size == 9, "AdjacencyRelation stencil size")
    require(type(adjacency.getAdjPixels(1, 1)).__name__ == "AdjacencyRelation", "AdjacencyRelation getAdjPixels return type")

    require_raises(
        lambda: mmcfilters.MorphologicalTree(image, True),
        "direct MorphologicalTree construction should not be public",
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

    weighted_level_attr = mmcfilters.Attribute.computeSingleAttribute(weighted, mmcfilters.Attribute.LEVEL)
    weighted_level_attr64 = mmcfilters.Attribute.computeSingleAttribute(weighted, mmcfilters.Attribute.LEVEL, dtype=np.float64)
    require(np.array_equal(weighted_filters.filteringByExtinction(weighted_level_attr, 1024), weighted_reconstruction), "weighted AttributeFilters extinction filtering keep-all")
    require(np.array_equal(weighted_filters.filteringByExtinction(weighted_level_attr64, 1024), weighted_reconstruction), "weighted AttributeFilters extinction filtering float64")
    require(weighted_filters.saliencyMapByExtinction(weighted_level_attr, 1024).shape == (4, 4), "weighted AttributeFilters extinction saliency shape")
    saliency64 = weighted_filters.saliencyMapByExtinction(weighted_level_attr64, 1024)
    require(saliency64.shape == (4, 4), "weighted AttributeFilters extinction float64 saliency shape")
    require(saliency64.dtype == np.float64, "weighted AttributeFilters extinction float64 saliency dtype")
    require_raises(lambda: weighted_filters.filteringByExtinction(np.array([1.0], dtype=np.float32), 1), "weighted filteringByExtinction must reject short attribute buffer")
    require_raises(lambda: weighted_filters.saliencyMapByExtinction(np.array([1.0], dtype=np.float32), 1), "weighted saliencyMapByExtinction must reject short attribute buffer")
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
    require(np.array_equal(weighted_extinction.filtering(1024), weighted_reconstruction), "weighted ExtinctionValues filtering keep-all")
    require(np.array_equal(weighted_extinction64.filtering(1024), weighted_reconstruction), "weighted ExtinctionValues float64 filtering keep-all")
    require(weighted_extinction.saliencyMap(1024).shape == (4, 4), "weighted ExtinctionValues saliency shape")
    require(weighted_extinction64.saliencyMap(1024).dtype == np.float64, "weighted ExtinctionValues float64 saliency dtype")
    require_raises(
        lambda: mmcfilters.ExtinctionValues(weighted, np.array([1.0], dtype=np.float32)),
        "weighted ExtinctionValues must reject short attribute buffer",
    )

    stale_weighted = mmcfilters.MorphologicalTreeFactory.createMaxTree(image)
    stale_keep_all = [True] * stale_weighted.numInternalNodeSlots
    stale_level_attr = mmcfilters.Attribute.computeSingleAttribute(stale_weighted, mmcfilters.Attribute.LEVEL)
    stale_filters = mmcfilters.AttributeFilters(stale_weighted)
    stale_extinction = mmcfilters.ExtinctionValues(stale_weighted, stale_level_attr)
    stale_uao = mmcfilters.UltimateAttributeOpening(stale_weighted, stale_level_attr)
    stale_weighted.mergeNodeIntoParent(4)
    require_raises(lambda: stale_filters.filteringDirectRule(stale_keep_all), "AttributeFilters must reject use after topology mutation")
    require_raises(lambda: stale_extinction.filtering(1024), "ExtinctionValues must reject filtering after topology mutation")
    require_raises(lambda: stale_extinction.saliencyMap(1024), "ExtinctionValues must reject saliency after topology mutation")
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
