#!/usr/bin/env python3

import gc
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


def require_raises(callback, message: str):
    try:
        callback()
    except (TypeError, ValueError):
        return
    raise RuntimeError(message)


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: test_python_contour_traces.py <build-dir>")

    mmcfilters = load_native_module(pathlib.Path(sys.argv[1]).resolve())
    require(hasattr(mmcfilters, "ContourTraceComputation"), "package must expose ContourTraceComputation")
    require(hasattr(mmcfilters, "ContourLoopKind"), "package must expose ContourLoopKind")
    require(hasattr(mmcfilters, "ContourTraceSide"), "package must expose ContourTraceSide")
    require_raises(
        lambda: mmcfilters.ContourTraceComputation.extraction(None),
        "contour trace extraction must reject None",
    )

    image = np.array(
        [
            [2, 2, 2],
            [2, 1, 2],
            [2, 2, 2],
        ],
        dtype=np.uint8,
    )
    tree = mmcfilters.MorphologicalTreeFactory.createMaxTree(image)
    area = mmcfilters.Attribute.computeSingleTopologyAttribute(tree, mmcfilters.Attribute.AREA)
    ring_nodes = [node for node, value in enumerate(area.tolist()) if int(value) == 8]
    require(len(ring_nodes) == 1, "ring fixture must expose one area-8 node")

    tree_reference_count = sys.getrefcount(tree)
    traces = mmcfilters.ContourTraceComputation.extraction(tree)
    require(
        sys.getrefcount(tree) > tree_reference_count,
        "ContourTraces must retain a strong reference to its source tree",
    )
    del tree
    gc.collect()

    # ContourTraces borrows the topology internally, so its Python binding must
    # keep the source tree alive for every lazy query.
    require(traces.isMaterialized is False, "traces must start without global materialization")
    edges = traces.getEdges(ring_nodes[0])
    require(len(edges) == 16, "ring trace edge count")
    require(traces.isEdgeMaterialized(ring_nodes[0]) is True, "getEdges must materialize node edges")
    require(traces.isNodeTraced(ring_nodes[0]) is False, "getEdges must not trace loops")

    loops = traces.getLoops(ring_nodes[0])
    require(len(loops) == 2, "ring trace loop count")
    require(traces.isNodeTraced(ring_nodes[0]) is True, "getLoops must trace node loops")
    external = [loop for loop in loops if loop.kind == mmcfilters.ContourLoopKind.External]
    internal = [loop for loop in loops if loop.kind == mmcfilters.ContourLoopKind.Internal]
    require(len(external) == 1, "ring must have one external loop")
    require(len(internal) == 1, "ring must have one internal loop")
    require(external[0].edge_count == 12, "ring external loop edge count")
    require(internal[0].edge_count == 4, "ring internal loop edge count")
    require(external[0].signed_area2 > 0, "external loop signed area")
    require(internal[0].signed_area2 < 0, "internal loop signed area")
    require(len(traces.getLoopEdges(external[0])) == 12, "external loop edge range")
    require(len(traces.getLoopEdges(internal[0])) == 4, "internal loop edge range")

    traces.materializeAll()
    require(traces.isMaterialized is True, "materializeAll must trace all nodes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
