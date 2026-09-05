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


def trace_signature(trace):
    signature = []
    for boundary in trace.boundaries():
        edges = tuple((edge.pixel, str(edge.side)) for edge in trace.boundary_edges(boundary))
        first_edge = min(range(len(edges)), key=edges.__getitem__)
        normalized_edges = edges[first_edge:] + edges[:first_edge]
        signature.append(
            (
                boundary.doubled_signed_area < 0,
                boundary.edge_count,
                boundary.doubled_signed_area,
                normalized_edges,
            )
        )
    return tuple(sorted(signature))


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: test_python_contour_traces.py <build-dir>")

    mmcfilters = load_native_module(pathlib.Path(sys.argv[1]).resolve())
    require(hasattr(mmcfilters, "ContourTraceComputation"), "package must expose ContourTraceComputation")
    require(hasattr(mmcfilters, "ContourTrace"), "package must expose ContourTrace")
    require(hasattr(mmcfilters, "ContourBoundaryKind"), "package must expose ContourBoundaryKind")
    require(hasattr(mmcfilters, "ContourSide"), "package must expose ContourSide")
    require_raises(
        lambda: mmcfilters.ContourTraceComputation(None),
        "contour trace construction must reject None",
    )

    image = np.array(
        [
            [2, 2, 2],
            [2, 1, 2],
            [2, 2, 2],
        ],
        dtype=np.uint8,
    )
    tree = mmcfilters.MorphologicalTreeFactory.create_max_tree(image)
    area = mmcfilters.Attribute.compute_single_topology_attribute(tree, mmcfilters.Attribute.AREA)
    ring_nodes = [node for node, value in enumerate(area.tolist()) if int(value) == 8]
    require(len(ring_nodes) == 1, "ring fixture must expose one area-8 node")

    tree_reference_count = sys.getrefcount(tree)
    traces = mmcfilters.ContourTraceComputation(tree)
    require(isinstance(traces, mmcfilters.ContourTraceComputation), "constructor must return the public computation type")
    require(
        sys.getrefcount(tree) > tree_reference_count,
        "ContourTraceComputation must retain a strong reference to its source tree",
    )
    num_nodes = tree.num_nodes
    del tree
    gc.collect()

    trace = traces.trace(ring_nodes[0])
    require(isinstance(trace, mmcfilters.ContourTrace), "trace must return an owned contour trace")
    require(len(trace.edges()) == 16, "ring trace edge count")
    boundaries = trace.boundaries()
    require(len(boundaries) == 2, "ring trace boundary count")
    external = [boundary for boundary in boundaries if boundary.kind == mmcfilters.ContourBoundaryKind.EXTERNAL]
    internal = [boundary for boundary in boundaries if boundary.kind == mmcfilters.ContourBoundaryKind.INTERNAL]
    require(len(external) == 1, "ring must have one external boundary")
    require(len(internal) == 1, "ring must have one internal boundary")
    require(external[0].edge_count == 12, "ring external boundary edge count")
    require(internal[0].edge_count == 4, "ring internal boundary edge count")
    require(external[0].doubled_signed_area > 0, "external boundary signed area")
    require(internal[0].doubled_signed_area < 0, "internal boundary signed area")
    require(len(trace.boundary_edges(external[0])) == 12, "external boundary edge range")
    require(len(trace.boundary_edges(internal[0])) == 4, "internal boundary edge range")
    direct_external = trace.external_boundary()
    require(direct_external.kind == mmcfilters.ContourBoundaryKind.EXTERNAL, "direct external boundary kind")
    ordered_edges = trace.boundary_edges(direct_external)
    ordered_pixels = trace.boundary_pixels(direct_external)
    require(ordered_pixels == [edge.pixel for edge in ordered_edges], "boundary pixel projection must preserve edge order")
    require(len(set(ordered_pixels)) < len(ordered_pixels), "boundary pixel projection must retain repeated pixels")

    iteration_signatures = {}
    retained_trace = None
    for node, current_trace in traces:
        iteration_signatures[node] = trace_signature(current_trace)
        if retained_trace is None:
            retained_trace = current_trace
    require(len(iteration_signatures) == num_nodes, "iteration must visit every live node")
    require(trace_signature(retained_trace), "an iterated trace must retain its data after traversal")
    for node, signature in iteration_signatures.items():
        require(trace_signature(traces.trace(node)) == signature, "iteration and node query must agree")

    callback_signatures = {}
    traces.for_each_trace(lambda node, current_trace: callback_signatures.update({node: trace_signature(current_trace)}))
    require(callback_signatures == iteration_signatures, "callback and iterator traversal must agree")

    # The center zero reaches the exterior only under 8-connected background.
    # A repeated grid vertex must not prematurely close a 4/8 contour.
    image = np.array([[1, 0, 2], [0, 2, 0], [0, 0, 0]], dtype=np.uint8)
    for radius, expected in ((1.0, [(16, 14)]), (1.5, [(4, -2), (12, 16)])):
        tree = mmcfilters.MorphologicalTreeFactory.create_min_tree(image, radius)
        area = mmcfilters.Attribute.compute_single_topology_attribute(tree, mmcfilters.Attribute.AREA)
        node = next(node for node, value in enumerate(area.tolist()) if int(value) == 7)
        trace = mmcfilters.ContourTraceComputation(tree).trace(node)
        actual = sorted((boundary.edge_count, boundary.doubled_signed_area) for boundary in trace.boundaries())
        require(actual == expected, f"complementary contour pairing for radius={radius}: {actual}")

    # A valued view supplies the polarity needed to trace complementary shapes.
    image = np.ones((5, 5), dtype=np.uint8)
    image[1, 1] = image[2, 2] = 2
    tree = mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(image)
    area = mmcfilters.Attribute.compute_single_topology_attribute(tree, mmcfilters.Attribute.AREA)
    node = next(node for node, value in enumerate(area.tolist()) if int(value) == 2)
    trace = mmcfilters.ContourTraceComputation(tree).trace(node)
    boundaries = trace.boundaries()
    require(len(boundaries) == 1 and boundaries[0].edge_count == 8, "upper 8-connected shape must have one self-touching cycle")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
