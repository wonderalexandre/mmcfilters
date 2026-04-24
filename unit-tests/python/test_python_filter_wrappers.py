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

    weighted = mmcfilters.WeightedMorphologicalTree(image, True)
    weighted_reconstruction = weighted.reconstructionImage()

    adjacency = mmcfilters.AdjacencyRelation(4, 4, 1.5)
    require(adjacency.size == 9, "AdjacencyRelation stencil size")
    require(type(adjacency.getAdjPixels(1, 1)).__name__ == "AdjacencyRelation", "AdjacencyRelation getAdjPixels return type")

    tree = mmcfilters.MorphologicalTree(image, True)
    keep_all = [True] * weighted.numInternalNodeSlots
    require_raises(
        lambda: mmcfilters.AttributeFilters(tree).filteringDirectRule(keep_all),
        "topology-only MorphologicalTree must reject weighted filtering operations",
    )

    weighted_filters = mmcfilters.AttributeFilters(weighted)
    require(np.array_equal(weighted_filters.filteringDirectRule(keep_all), weighted_reconstruction), "weighted AttributeFilters direct rule keep-all")
    require(np.array_equal(weighted_filters.filteringSubtractiveRule(keep_all), weighted_reconstruction), "weighted AttributeFilters subtractive rule keep-all")
    require(np.array_equal(weighted_filters.filteringMin(keep_all), weighted_reconstruction), "weighted AttributeFilters pruning min keep-all criterion")
    require(np.array_equal(weighted_filters.filteringMax(keep_all), weighted_reconstruction), "weighted AttributeFilters pruning max keep-all criterion")
    weighted_box_height = mmcfilters.Attribute.computeSingleAttribute(weighted, mmcfilters.Attribute.BOX_HEIGHT)
    require(np.array_equal(weighted_filters.filteringMin(weighted_box_height, 1.0), weighted_reconstruction), "weighted AttributeFilters pruning min keep-all attribute")
    require(np.array_equal(weighted_filters.filteringMax(weighted_box_height, 1.0), weighted_reconstruction), "weighted AttributeFilters pruning max keep-all attribute")
    require_raises(lambda: weighted_filters.filteringMin(np.array([1.0], dtype=np.float32), 1.0), "weighted filteringMin must reject short attribute buffer")
    require_raises(lambda: weighted_filters.filteringDirectRule([True]), "weighted filteringDirectRule must reject short criterion")
    require(weighted_filters.getAdaptiveCriterion(keep_all, 2) == [False] * weighted.numInternalNodeSlots, "weighted AttributeFilters adaptive criterion on all-true input")

    weighted_level_attr = mmcfilters.Attribute.computeSingleAttribute(weighted, mmcfilters.Attribute.LEVEL)
    require(np.array_equal(weighted_filters.filteringByExtinction(weighted_level_attr, 1024), weighted_reconstruction), "weighted AttributeFilters extinction filtering keep-all")
    require(weighted_filters.saliencyMapByExtinction(weighted_level_attr, 1024).shape == (4, 4), "weighted AttributeFilters extinction saliency shape")
    require_raises(lambda: weighted_filters.filteringByExtinction(np.array([1.0], dtype=np.float32), 1), "weighted filteringByExtinction must reject short attribute buffer")
    require_raises(lambda: weighted_filters.saliencyMapByExtinction(np.array([1.0], dtype=np.float32), 1), "weighted saliencyMapByExtinction must reject short attribute buffer")

    weighted_primitives = mmcfilters.AttributeOpeningPrimitivesFamily(weighted, weighted_box_height, float(weighted.numRows))
    require(weighted_primitives.numPrimitives >= 1, "weighted AttributeOpeningPrimitivesFamily numPrimitives")
    require(len(weighted_primitives.getThresholdsPrimitive()) == weighted_primitives.numPrimitives, "weighted AttributeOpeningPrimitivesFamily threshold count")
    require(weighted_primitives.getNodesWithMaximumCriterium() == [4], "weighted AttributeOpeningPrimitivesFamily nodes with maximum criterion")
    require(weighted_primitives.getPrimitive(float(weighted.numRows)).shape == (4, 4), "weighted AttributeOpeningPrimitivesFamily primitive image shape")
    require(weighted_primitives.restOfImage.shape == (4, 4), "weighted AttributeOpeningPrimitivesFamily rest image shape")
    require_raises(
        lambda: mmcfilters.AttributeOpeningPrimitivesFamily(weighted, np.array([1.0], dtype=np.float32), float(weighted.numRows)),
        "weighted AttributeOpeningPrimitivesFamily must reject short attribute buffer",
    )

    weighted_uao = mmcfilters.UltimateAttributeOpening(weighted, weighted_box_height)
    weighted_uao.execute(int(weighted.numRows))
    require(weighted_uao.getMaxContrastImage().shape == (4, 4), "weighted UltimateAttributeOpening max contrast shape")
    require(weighted_uao.getAssociatedImage().shape == (4, 4), "weighted UltimateAttributeOpening associated image shape")
    require(weighted_uao.getAssociatedColoredImage().shape == (4, 12), "weighted UltimateAttributeOpening associated color image shape")
    weighted_uao.executeWithMSER(int(weighted.numRows), 1)
    require(weighted_uao.getMaxContrastImage().shape == (4, 4), "weighted UltimateAttributeOpening MSER execute shape")
    require_raises(
        lambda: mmcfilters.UltimateAttributeOpening(weighted, np.array([1.0], dtype=np.float32)),
        "weighted UltimateAttributeOpening must reject short attribute buffer",
    )

    weighted_extinction = mmcfilters.ExtinctionValues(weighted, weighted_level_attr)
    require(np.array_equal(weighted_extinction.filtering(1024), weighted_reconstruction), "weighted ExtinctionValues filtering keep-all")
    require(weighted_extinction.saliencyMap(1024).shape == (4, 4), "weighted ExtinctionValues saliency shape")
    require_raises(
        lambda: mmcfilters.ExtinctionValues(weighted, np.array([1.0], dtype=np.float32)),
        "weighted ExtinctionValues must reject short attribute buffer",
    )

    print("python filter wrappers ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
