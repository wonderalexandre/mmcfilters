#!/usr/bin/env python3
"""Benchmark representative calls through the public Python API."""

from __future__ import annotations

import argparse
import hashlib
import json
import statistics
import time
from dataclasses import dataclass
from typing import Any, Callable

import mmcfilters
import numpy as np


@dataclass(frozen=True)
class Scenario:
    name: str
    median_ms: float
    mad_ms: float
    minimum_ms: float
    maximum_ms: float
    checksum: str


def make_input(pattern: str, rows: int, cols: int) -> np.ndarray:
    row, col = np.indices((rows, cols), dtype=np.int64)
    if pattern == "structured":
        radial = (row - rows // 2) ** 2 + (col - cols // 2) ** 2
        values = radial + 17 * row + 31 * col
    elif pattern == "noise":
        values = (row * cols + col).astype(np.uint32) + np.uint32(0x9E3779B9)
        values ^= values >> np.uint32(16)
        values *= np.uint32(0x7FEB352D)
        values ^= values >> np.uint32(15)
        values *= np.uint32(0x846CA68B)
        values ^= values >> np.uint32(16)
    elif pattern == "ramp":
        values = 3 * row + 5 * col
    elif pattern == "geometric":
        values = np.full((rows, cols), 24, dtype=np.int64)
        radius = max(2, min(rows, cols) // 4)
        values[(row - rows // 2) ** 2 + (col - cols // 2) ** 2 <= radius**2] = 180
        values[rows // 8 : rows // 3, cols // 7 : cols // 2] = 92
        values[rows // 2 : 7 * rows // 8, 5 * cols // 8 : 7 * cols // 8] = 232
    else:
        values = np.full((rows, cols), 127, dtype=np.int64)
    return np.ascontiguousarray(values.astype(np.uint8))


def append_array(digest: Any, value: Any) -> None:
    array = np.ascontiguousarray(value)
    digest.update(str(array.dtype).encode())
    digest.update(repr(array.shape).encode())
    digest.update(array.tobytes())


def array_checksum(value: Any) -> str:
    digest = hashlib.blake2b(digest_size=8)
    append_array(digest, value)
    return digest.hexdigest()


def tree_checksum(tree: Any) -> str:
    parent, altitude = tree.exportHigraHierarchy()
    digest = hashlib.blake2b(digest_size=8)
    append_array(digest, parent)
    append_array(digest, altitude)
    return digest.hexdigest()


def casf_checksum(casf: Any, image: Any | None = None) -> str:
    digest = hashlib.blake2b(digest_size=8)
    if image is not None:
        append_array(digest, image)
    for hierarchy in (casf.exportMinTree(), casf.exportMaxTree()):
        append_array(digest, hierarchy[0])
        append_array(digest, hierarchy[1])
    return digest.hexdigest()


def canonical_attributes_checksum(names: Any, values: Any) -> str:
    matrix = np.asarray(values)
    if matrix.ndim == 1:
        matrix = matrix[:, None]
    if isinstance(names, dict):
        columns = [(name, matrix[:, index]) for name, index in names.items()]
    else:
        columns = list(zip(names, matrix.T, strict=True))
    digest = hashlib.blake2b(digest_size=8)
    for name, column in sorted(columns, key=lambda item: item[0]):
        digest.update(name.encode())
        append_array(digest, column)
    return digest.hexdigest()


def benchmark(name: str, repetitions: int, operation: Callable[[], Any], checksum_of: Callable[[Any], str]) -> Scenario:
    expected = checksum_of(operation())
    samples: list[float] = []
    for _ in range(repetitions):
        start = time.perf_counter_ns()
        result = operation()
        finish = time.perf_counter_ns()
        checksum = checksum_of(result)
        if checksum != expected:
            raise RuntimeError(f"{name} produced a non-deterministic result")
        samples.append((finish - start) / 1_000_000.0)
    center = statistics.median(samples)
    return Scenario(name, center, statistics.median(abs(value - center) for value in samples), min(samples), max(samples), expected)


def benchmark_prepared(
    name: str,
    repetitions: int,
    prepare: Callable[[], Any],
    operation: Callable[[Any], Any],
    checksum_of: Callable[[Any, Any], str],
) -> Scenario:
    state = prepare()
    expected = checksum_of(state, operation(state))
    samples: list[float] = []
    for _ in range(repetitions):
        state = prepare()
        start = time.perf_counter_ns()
        result = operation(state)
        finish = time.perf_counter_ns()
        checksum = checksum_of(state, result)
        if checksum != expected:
            raise RuntimeError(f"{name} produced a non-deterministic result")
        samples.append((finish - start) / 1_000_000.0)
    center = statistics.median(samples)
    return Scenario(name, center, statistics.median(abs(value - center) for value in samples), min(samples), max(samples), expected)


def group_requests(profile: str) -> list[tuple[str, Any]]:
    groups = [
        ("gray_level", mmcfilters.Attribute.Group.GRAY_LEVEL),
        ("tree_topology", mmcfilters.Attribute.Group.TREE_TOPOLOGY),
    ]
    if profile != "smoke":
        groups[1:1] = [
            ("shape", mmcfilters.Attribute.Group.SHAPE),
            ("moments", mmcfilters.Attribute.Group.MOMENTS),
            ("boundary", mmcfilters.Attribute.Group.BOUNDARY),
        ]
        groups.append(("all", mmcfilters.Attribute.ALL))
    return groups


def pipeline_area_direct(image: np.ndarray) -> tuple[Any, Any, np.ndarray, np.ndarray]:
    tree = mmcfilters.MorphologicalTreeFactory.createMaxTree(image, radius=1.5)
    names, attributes = mmcfilters.Attribute.computeAttributes(tree, [mmcfilters.Attribute.Group.SHAPE])
    area = attributes[:, names["AREA"]]
    criterion = (area >= image.size / 16.0).tolist()
    criterion[tree.getRoot()] = True
    filtered = mmcfilters.AttributeFilters(tree).filteringDirectRule(criterion)
    return tree, names, attributes, filtered


def pipeline_checksum(result: tuple[Any, Any, np.ndarray, np.ndarray]) -> str:
    tree, names, attributes, image = result
    digest = hashlib.blake2b(digest_size=8)
    digest.update(tree_checksum(tree).encode())
    digest.update(canonical_attributes_checksum(names, attributes).encode())
    append_array(digest, image)
    return digest.hexdigest()


def parse_quantiles(value: str) -> tuple[float, float, float]:
    try:
        quantiles = tuple(float(item) for item in value.split(","))
    except ValueError as error:
        raise argparse.ArgumentTypeError("CASF quantiles must be numeric") from error
    if len(quantiles) != 3 or any(not 0.0 < item < 1.0 for item in quantiles) or tuple(sorted(set(quantiles))) != quantiles:
        raise argparse.ArgumentTypeError("CASF quantiles must contain three strictly increasing values between zero and one")
    return quantiles


def casf_area_thresholds(image: np.ndarray, quantiles: tuple[float, float, float]) -> list[float]:
    values: list[float] = []
    trees = (
        mmcfilters.MorphologicalTreeFactory.createMinTree(image, radius=1.5),
        mmcfilters.MorphologicalTreeFactory.createMaxTree(image, radius=1.5),
    )
    for tree in trees:
        area = mmcfilters.Attribute.computeSingleAttribute(tree, mmcfilters.Attribute.AREA, dtype=np.float64)
        root = tree.getRoot()
        values.extend(float(area[node]) for node in tree.getAliveNodeIds() if node != root)
    if not values:
        return [0.0, 0.0, 0.0]
    ordered = np.sort(np.asarray(values, dtype=np.float64))
    thresholds: list[float] = []
    for quantile in quantiles:
        index = int(np.floor(quantile * (ordered.size - 1)))
        if thresholds and ordered[index] <= thresholds[-1]:
            index = int(np.searchsorted(ordered, thresholds[-1], side="right"))
            index = min(index, ordered.size - 1)
        thresholds.append(float(ordered[index]))
    return thresholds


def run(options: argparse.Namespace) -> list[Scenario]:
    image = make_input(options.input, options.rows, options.cols)
    tree = mmcfilters.MorphologicalTreeFactory.createMaxTree(image, radius=1.5)
    scenarios = [
        benchmark("python.construction.end_to_end.max_tree", options.repetitions, lambda: mmcfilters.MorphologicalTreeFactory.createMaxTree(image, radius=1.5), tree_checksum),
        benchmark("python.construction.end_to_end.min_tree", options.repetitions, lambda: mmcfilters.MorphologicalTreeFactory.createMinTree(image, radius=1.5), tree_checksum),
    ]

    for attribute_name in ("AREA", "MAX_DIST"):
        attribute = getattr(mmcfilters.Attribute, attribute_name)
        scenarios.append(
            benchmark(
                f"python.attributes.established_input.scalar_{attribute_name.lower()}",
                options.repetitions,
                lambda attribute=attribute: mmcfilters.Attribute.computeSingleAttribute(tree, attribute, dtype=np.float64),
                array_checksum,
            )
        )

    for group_name, group in group_requests(options.profile):
        grouped = benchmark(
            f"python.attribute_groups.established_input.{group_name}",
            options.repetitions,
            lambda group=group: mmcfilters.Attribute.computeAttributes(tree, [group], dtype=np.float64),
            lambda result: canonical_attributes_checksum(result[0], result[1]),
        )
        scenarios.append(grouped)
        if options.profile != "smoke":
            layout, _ = mmcfilters.Attribute.computeAttributes(tree, [group], dtype=np.float64)
            names = list(layout)

            def sequential(names: list[str] = names) -> tuple[list[str], np.ndarray]:
                columns = [
                    mmcfilters.Attribute.computeSingleAttribute(tree, getattr(mmcfilters.Attribute, name), dtype=np.float64)
                    for name in names
                ]
                return names, np.column_stack(columns)

            scalar = benchmark(
                f"python.attribute_groups.established_input.{group_name}_sequential_scalars",
                options.repetitions,
                sequential,
                lambda result: canonical_attributes_checksum(result[0], result[1]),
            )
            if scalar.checksum != grouped.checksum:
                raise RuntimeError(f"grouped and sequential Python results differ for {group_name}")
            scenarios.append(scalar)

    area = mmcfilters.Attribute.computeSingleAttribute(tree, mmcfilters.Attribute.AREA, dtype=np.float64)
    scenarios.append(
        benchmark(
            "python.filters.established_input.pruning_max_area",
            options.repetitions,
            lambda: mmcfilters.AttributeFilters(tree).filteringByPruningMax(area, image.size / 16.0),
            array_checksum,
        )
    )
    thresholds = casf_area_thresholds(image, options.casf_quantiles)
    scenarios.append(
        benchmark(
            "python.casf.end_to_end.construction_area",
            options.repetitions,
            lambda: mmcfilters.CasfComponentTrees(image, mmcfilters.CasfComponentTreesAttribute.AREA, radius=1.5),
            casf_checksum,
        )
    )
    scenarios.append(
        benchmark_prepared(
            "python.casf.established_input.incremental_area_sequence",
            options.repetitions,
            lambda: mmcfilters.CasfComponentTrees(image, mmcfilters.CasfComponentTreesAttribute.AREA, radius=1.5),
            lambda casf: casf.filter(thresholds),
            lambda casf, result: casf_checksum(casf, result),
        )
    )
    scenarios.append(benchmark("python.pipelines.end_to_end.max_tree_shape_group_direct_area", options.repetitions, lambda: pipeline_area_direct(image), pipeline_checksum))
    return scenarios


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=("smoke", "core", "publication"), default="smoke")
    parser.add_argument("--input", choices=("structured", "noise", "ramp", "geometric", "flat"), default="structured")
    parser.add_argument("--rows", type=int)
    parser.add_argument("--cols", type=int)
    parser.add_argument("--repetitions", type=int)
    parser.add_argument("--casf-quantiles", type=parse_quantiles, default=(0.1, 0.5, 0.9))
    parser.add_argument("--format", choices=("key-value", "jsonl"), default="key-value")
    options = parser.parse_args()
    defaults = {"smoke": (48, 48, 2), "core": (192, 192, 5), "publication": (512, 512, 15)}[options.profile]
    options.rows = options.rows or defaults[0]
    options.cols = options.cols or defaults[1]
    options.repetitions = options.repetitions or defaults[2]
    if min(options.rows, options.cols, options.repetitions) <= 0:
        parser.error("rows, cols, and repetitions must be positive")
    return options


def main() -> int:
    options = parse_arguments()
    scenarios = run(options)
    metadata = {
        "benchmark": "python_scientific_api",
        "contract_mode": "PYTHON_FACADE",
        "profile": options.profile,
        "input": options.input,
        "rows": options.rows,
        "cols": options.cols,
        "repetitions": options.repetitions,
        "casf_quantiles": ",".join(f"{value:.6f}" for value in options.casf_quantiles),
        "scenario_count": len(scenarios),
    }
    if options.format == "jsonl":
        print(json.dumps({"record": "metadata", **metadata}, separators=(",", ":")))
        for scenario in scenarios:
            print(json.dumps({"record": "scenario", **scenario.__dict__}, separators=(",", ":")))
    else:
        for key, value in metadata.items():
            print(f"{key}={value}")
        for scenario in scenarios:
            print(f"{scenario.name}_median_ms={scenario.median_ms:.6f}")
            print(f"{scenario.name}_mad_ms={scenario.mad_ms:.6f}")
            print(f"{scenario.name}_minimum_ms={scenario.minimum_ms:.6f}")
            print(f"{scenario.name}_maximum_ms={scenario.maximum_ms:.6f}")
            print(f"{scenario.name}_checksum={scenario.checksum}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
