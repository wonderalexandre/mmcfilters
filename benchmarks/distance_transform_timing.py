#!/usr/bin/env python3
"""Benchmark DIST_TRANSF versus DIST_TRANSF_EXACT and their maximum attributes."""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import itertools
import json
import math
import platform
import statistics
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from PIL import Image, __version__ as PILLOW_VERSION

import mmcfilters

from distance_transform_sensitivity import build_tree, discover_images, git_metadata, parse_tree_kinds, sha256_file


RESOLUTIONS = {
    "480p": (480, 854),
    "720p": (720, 1280),
    "1080p": (1080, 1920),
}
MODES = ("dist_transf", "dist_transf_exact", "max_dist", "max_dist_exact")
EXACT_SUFFIX = "_EXACT"
OBSERVATION_FIELDS = (
    "resolution",
    "image_index",
    "image",
    "tree",
    "repetition",
    "order_position",
    "mode",
    "rows",
    "columns",
    "pixels",
    "nodes",
    "seconds",
    "milliseconds",
    "nanoseconds_per_pixel",
    "nanoseconds_per_node",
    "megapixels_per_second",
    "meganodes_per_second",
    "output_columns",
    "checksum",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("image_directory", type=Path)
    parser.add_argument("output_directory", type=Path)
    parser.add_argument("--count", type=int, default=10)
    parser.add_argument(
        "--trees",
        default="max,min,tos,residual_unrestricted,residual_saturated",
        help="Comma-separated tree kinds (default: all five studied kinds).",
    )
    parser.add_argument("--resolutions", default="480p,720p,1080p")
    parser.add_argument("--modes", default=",".join(MODES))
    parser.add_argument("--repetitions", type=int, default=1)
    parser.add_argument("--radius", type=float, default=1.5)
    parser.add_argument("--infinity-pixel", type=int, default=0)
    parser.add_argument("--bootstrap-replicates", type=int, default=20_000)
    parser.add_argument("--expected-contract-mode", choices=("CHECKED", "UNCHECKED"))
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def parse_choices(value: str, allowed: Iterable[str], label: str) -> list[str]:
    result = [item.strip().lower() for item in value.split(",") if item.strip()]
    allowed_set = set(allowed)
    if not result or any(item not in allowed_set for item in result):
        raise ValueError(f"{label} must contain values from: {','.join(allowed)}")
    if len(set(result)) != len(result):
        raise ValueError(f"{label} must not contain duplicates")
    return result


def resized_image(path: Path, rows: int, columns: int) -> np.ndarray:
    with Image.open(path) as source:
        gray = source.convert("L")
        resized = gray.resize((columns, rows), resample=Image.Resampling.LANCZOS)
        return np.ascontiguousarray(np.asarray(resized, dtype=np.uint8))


def array_sha256(array: np.ndarray) -> str:
    digest = hashlib.sha256()
    digest.update(f"{array.shape[0]}x{array.shape[1]}:{array.dtype}:".encode())
    digest.update(memoryview(np.ascontiguousarray(array)))
    return digest.hexdigest()


def file_sha256(path: Path) -> str:
    return sha256_file(path)


def module_metadata() -> dict[str, Any]:
    module_path = Path(mmcfilters._native.__file__).resolve()
    cache_path: Path | None = None
    for parent in module_path.parents:
        candidate = parent / "CMakeCache.txt"
        if candidate.is_file():
            cache_path = candidate
            break

    cache_values: dict[str, str] = {}
    if cache_path is not None:
        wanted = {
            "CMAKE_BUILD_TYPE",
            "CMAKE_CXX_COMPILER",
            "CMAKE_CXX_FLAGS_RELEASE",
            "CMAKE_OSX_ARCHITECTURES",
            "MMCFILTERS_CONTRACT_MODE",
        }
        for line in cache_path.read_text(encoding="utf-8", errors="replace").splitlines():
            if ":" not in line or "=" not in line:
                continue
            key = line.split(":", 1)[0]
            if key in wanted:
                cache_values[key] = line.split("=", 1)[1]
    return {
        "native_module": str(module_path),
        "native_module_sha256": file_sha256(module_path),
        "cmake_cache": str(cache_path) if cache_path is not None else None,
        "cmake": cache_values,
    }


def cpu_name() -> str:
    try:
        result = subprocess.run(
            ["sysctl", "-n", "machdep.cpu.brand_string"],
            check=True,
            capture_output=True,
            text=True,
        )
        value = result.stdout.strip()
        if value:
            return value
    except (OSError, subprocess.CalledProcessError):
        pass
    return platform.processor() or platform.machine()


def discover_requests() -> tuple[dict[str, Any], dict[str, list[str]]]:
    tiny = np.arange(12 * 16, dtype=np.uint8).reshape(12, 16)
    tree = mmcfilters.MorphologicalTreeFactory.create_max_tree(tiny, radius=1.5)
    approximate_layout, _ = mmcfilters.Attribute.compute_topology_attributes(
        tree,
        [mmcfilters.Attribute.Group.DIST_TRANSF],
        dtype=np.float32,
    )
    exact_layout, _ = mmcfilters.Attribute.compute_topology_attributes(
        tree,
        [mmcfilters.Attribute.Group.DIST_TRANSF_EXACT],
        dtype=np.float32,
    )
    approximate_names = list(approximate_layout)
    exact_names = list(exact_layout)
    if len(approximate_names) != 29 or any(name.endswith(EXACT_SUFFIX) for name in approximate_names):
        raise RuntimeError("DIST_TRANSF must expose exactly 29 unsuffixed approximate attributes")
    if len(exact_names) != 29 or any(not name.endswith(EXACT_SUFFIX) for name in exact_names):
        raise RuntimeError("DIST_TRANSF_EXACT must expose exactly 29 exact attributes")
    if {name + EXACT_SUFFIX for name in approximate_names} != set(exact_names):
        raise RuntimeError("DIST_TRANSF and DIST_TRANSF_EXACT must form aligned attribute pairs")
    requests = {
        "dist_transf": [mmcfilters.Attribute.Group.DIST_TRANSF],
        "dist_transf_exact": [mmcfilters.Attribute.Group.DIST_TRANSF_EXACT],
        "max_dist": [mmcfilters.Attribute.MAX_DIST],
        "max_dist_exact": [mmcfilters.Attribute.MAX_DIST_EXACT],
    }
    names = {
        "dist_transf": approximate_names,
        "dist_transf_exact": exact_names,
        "max_dist": ["MAX_DIST"],
        "max_dist_exact": ["MAX_DIST_EXACT"],
    }
    return requests, names


def warm_up(tree_kinds: list[str], modes: list[str], requests: dict[str, Any], radius: float, infinity_pixel: int) -> None:
    row, column = np.indices((24, 32), dtype=np.int64)
    image = np.ascontiguousarray(((17 * row + 29 * column + row * column) % 256).astype(np.uint8))
    for tree_kind in tree_kinds:
        tree = build_tree(image, tree_kind, radius, infinity_pixel)
        for mode in modes:
            mmcfilters.Attribute.compute_topology_attributes(tree, requests[mode], dtype=np.float32)


def result_checksum(mode: str, layout: Any, values: np.ndarray, expected_names: list[str]) -> str:
    names = list(layout)
    if names != expected_names:
        raise RuntimeError(f"{mode} returned an unexpected attribute layout")
    matrix = np.ascontiguousarray(values)
    digest = hashlib.blake2b(digest_size=16)
    digest.update(mode.encode())
    digest.update("\0".join(names).encode())
    digest.update(str(matrix.dtype).encode())
    digest.update(repr(matrix.shape).encode())
    digest.update(matrix.tobytes(order="C"))
    return digest.hexdigest()


def append_observation(path: Path, row: dict[str, Any]) -> None:
    write_header = not path.exists() or path.stat().st_size == 0
    with path.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=OBSERVATION_FIELDS)
        if write_header:
            writer.writeheader()
        writer.writerow({field: row[field] for field in OBSERVATION_FIELDS})
        stream.flush()


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: Iterable[str] | None = None) -> None:
    if not rows:
        if fieldnames is None:
            raise ValueError(f"cannot infer columns for empty CSV: {path}")
        with path.open("w", newline="", encoding="utf-8") as stream:
            csv.DictWriter(stream, fieldnames=list(fieldnames)).writeheader()
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def observation_key(row: dict[str, Any]) -> tuple[str, int, str, int, str]:
    return (
        str(row["resolution"]),
        int(row["image_index"]),
        str(row["tree"]),
        int(row["repetition"]),
        str(row["mode"]),
    )


def geometric_mean(values: Iterable[float]) -> float:
    array = np.asarray(list(values), dtype=np.float64)
    if array.size == 0 or np.any(array <= 0.0):
        return math.nan
    return float(np.exp(np.mean(np.log(array))))


def bootstrap_geometric_mean(values: list[float], replicates: int, seed: int) -> tuple[float, float]:
    if len(values) < 2 or replicates <= 0:
        return math.nan, math.nan
    logs = np.log(np.asarray(values, dtype=np.float64))
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, logs.size, size=(replicates, logs.size))
    samples = np.exp(np.mean(logs[indices], axis=1))
    return float(np.quantile(samples, 0.025)), float(np.quantile(samples, 0.975))


def summarize_observations(rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[(row["resolution"], row["tree"], row["mode"])].append(row)
    result: list[dict[str, Any]] = []
    resolution_order = {name: index for index, name in enumerate(RESOLUTIONS)}
    for (resolution, tree, mode), group in sorted(
        groups.items(), key=lambda item: (resolution_order[item[0][0]], item[0][1], MODES.index(item[0][2]))
    ):
        seconds = [float(row["seconds"]) for row in group]
        ns_pixel = [float(row["nanoseconds_per_pixel"]) for row in group]
        ns_node = [float(row["nanoseconds_per_node"]) for row in group]
        result.append({
            "resolution": resolution,
            "tree": tree,
            "mode": mode,
            "images": len({row["image"] for row in group}),
            "samples": len(group),
            "median_seconds": statistics.median(seconds),
            "mean_seconds": statistics.fmean(seconds),
            "stdev_seconds": statistics.pstdev(seconds),
            "minimum_seconds": min(seconds),
            "maximum_seconds": max(seconds),
            "median_nanoseconds_per_pixel": statistics.median(ns_pixel),
            "median_nanoseconds_per_node": statistics.median(ns_node),
            "median_megapixels_per_second": statistics.median(float(row["megapixels_per_second"]) for row in group),
            "median_meganodes_per_second": statistics.median(float(row["meganodes_per_second"]) for row in group),
        })
    return result


def summarize_builds(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[(str(row["resolution"]), str(row["tree"]))].append(row)
    result: list[dict[str, Any]] = []
    resolution_order = {name: index for index, name in enumerate(RESOLUTIONS)}
    for (resolution, tree), group in sorted(groups.items(), key=lambda item: (resolution_order[item[0][0]], item[0][1])):
        values = [float(row["build_seconds"]) for row in group]
        result.append({
            "resolution": resolution,
            "tree": tree,
            "images": len(group),
            "median_build_seconds": statistics.median(values),
            "mean_build_seconds": statistics.fmean(values),
            "minimum_build_seconds": min(values),
            "maximum_build_seconds": max(values),
            "median_nodes": statistics.median(int(row["nodes"]) for row in group),
        })
    return result


def paired_rows(observations: list[dict[str, str]], bootstrap_replicates: int) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    by_case: dict[tuple[str, int, str, int], dict[str, float]] = defaultdict(dict)
    image_names: dict[tuple[str, int, str, int], str] = {}
    for row in observations:
        key = (row["resolution"], int(row["image_index"]), row["tree"], int(row["repetition"]))
        by_case[key][row["mode"]] = float(row["seconds"])
        image_names[key] = row["image"]

    raw_pairs: list[dict[str, Any]] = []
    for (resolution, image_index, tree, repetition), values in sorted(by_case.items()):
        if set(values) != set(MODES):
            raise RuntimeError(f"incomplete timing modes for {resolution}/{image_index}/{tree}/{repetition}")
        raw_pairs.append({
            "resolution": resolution,
            "image_index": image_index,
            "image": image_names[(resolution, image_index, tree, repetition)],
            "tree": tree,
            "repetition": repetition,
            "dist_transf_exact_over_dist_transf": values["dist_transf_exact"] / values["dist_transf"],
            "max_dist_exact_over_max_dist": values["max_dist_exact"] / values["max_dist"],
        })

    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in raw_pairs:
        grouped[(row["resolution"], row["tree"])].append(row)

    summaries: list[dict[str, Any]] = []
    resolution_order = {name: index for index, name in enumerate(RESOLUTIONS)}
    for (resolution, tree), group in sorted(grouped.items(), key=lambda item: (resolution_order[item[0][0]], item[0][1])):
        by_image: dict[int, list[dict[str, Any]]] = defaultdict(list)
        for row in group:
            by_image[int(row["image_index"])].append(row)
        family_ratios = [geometric_mean(item["dist_transf_exact_over_dist_transf"] for item in image_group) for image_group in by_image.values()]
        maximum_ratios = [geometric_mean(item["max_dist_exact_over_max_dist"] for item in image_group) for image_group in by_image.values()]
        seed_material = f"{resolution}:{tree}".encode()
        seed = int.from_bytes(hashlib.sha256(seed_material).digest()[:8], "little")
        family_low, family_high = bootstrap_geometric_mean(family_ratios, bootstrap_replicates, seed)
        maximum_low, maximum_high = bootstrap_geometric_mean(maximum_ratios, bootstrap_replicates, seed ^ 0xA5A5A5A5)
        summaries.append({
            "resolution": resolution,
            "tree": tree,
            "images": len(by_image),
            "dist_transf_exact_over_dist_transf_geomean": geometric_mean(family_ratios),
            "dist_transf_exact_over_dist_transf_median": statistics.median(family_ratios),
            "dist_transf_exact_over_dist_transf_ci95_low": family_low,
            "dist_transf_exact_over_dist_transf_ci95_high": family_high,
            "max_dist_exact_over_max_dist_geomean": geometric_mean(maximum_ratios),
            "max_dist_exact_over_max_dist_median": statistics.median(maximum_ratios),
            "max_dist_exact_over_max_dist_ci95_low": maximum_low,
            "max_dist_exact_over_max_dist_ci95_high": maximum_high,
        })
    return raw_pairs, summaries


def scaling_rows(summary: list[dict[str, Any]]) -> list[dict[str, Any]]:
    lookup = {(row["resolution"], row["tree"], row["mode"]): row for row in summary}
    result: list[dict[str, Any]] = []
    pixel_480 = RESOLUTIONS["480p"][0] * RESOLUTIONS["480p"][1]
    pixel_720 = RESOLUTIONS["720p"][0] * RESOLUTIONS["720p"][1]
    pixel_1080 = RESOLUTIONS["1080p"][0] * RESOLUTIONS["1080p"][1]
    trees = sorted({row["tree"] for row in summary})
    modes = [mode for mode in MODES if any(row["mode"] == mode for row in summary)]
    for tree in trees:
        for mode in modes:
            keys = [(resolution, tree, mode) for resolution in RESOLUTIONS]
            if not all(key in lookup for key in keys):
                continue
            time_480 = float(lookup[("480p", tree, mode)]["median_seconds"])
            time_720 = float(lookup[("720p", tree, mode)]["median_seconds"])
            time_1080 = float(lookup[("1080p", tree, mode)]["median_seconds"])
            result.append({
                "tree": tree,
                "mode": mode,
                "time_ratio_720p_over_480p": time_720 / time_480,
                "time_ratio_1080p_over_720p": time_1080 / time_720,
                "time_ratio_1080p_over_480p": time_1080 / time_480,
                "pixel_ratio_720p_over_480p": pixel_720 / pixel_480,
                "pixel_ratio_1080p_over_720p": pixel_1080 / pixel_720,
                "pixel_ratio_1080p_over_480p": pixel_1080 / pixel_480,
                "scaling_exponent_480p_to_1080p": math.log(time_1080 / time_480) / math.log(pixel_1080 / pixel_480),
            })
    return result


def format_number(value: Any) -> str:
    number = float(value)
    if not math.isfinite(number):
        return "nan"
    if number == 0.0:
        return "0"
    if abs(number) >= 1000.0 or abs(number) < 0.001:
        return f"{number:.3e}"
    return f"{number:.3f}"


def build_report(
    metadata: dict[str, Any],
    summary: list[dict[str, Any]],
    pair_summary: list[dict[str, Any]],
    scaling: list[dict[str, Any]],
    builds: list[dict[str, Any]],
) -> str:
    timing_lookup = {(row["resolution"], row["tree"], row["mode"]): row for row in summary}
    pair_lookup = {(row["resolution"], row["tree"]): row for row in pair_summary}
    resolution_description = ", ".join(
        f"{name} ({metadata['configuration']['resolution_shapes'][name][1]}×{metadata['configuration']['resolution_shapes'][name][0]})"
        for name in metadata["configuration"]["resolutions"]
    )
    tos_protocol = []
    if "tos" in metadata["configuration"]["tree_kinds"]:
        tos_protocol.append(
            "- ToS: default minimum-4/maximum-8 complementary grid, no padding, infinity pixel zero, `uint8` altitudes."
        )
    lines = [
        "# Distance-transform group and maximum execution time on ICDAR",
        "",
        "## Protocol",
        "",
        f"- Inputs: first {metadata['configuration']['image_count']} ICDAR files in lexicographic order.",
        f"- Resolutions: {resolution_description}; Lanczos resizing is outside the timed region.",
        f"- Trees: `{', '.join(metadata['configuration']['tree_kinds'])}`; tree construction is reported separately.",
        *tos_protocol,
        f"- Contract mode: `{metadata['configuration']['contract_mode']}`; the runner rejects a module with a different mode.",
        "- Timed call: public `compute_topology_attributes(..., dtype=float32)` over an established tree.",
        "- Modes: `DIST_TRANSF` (29 approximate attributes), `DIST_TRANSF_EXACT` (29 exact attributes), `MAX_DIST`, and `MAX_DIST_EXACT`.",
        "- Every complete result matrix is hashed after the timed region.",
        f"- Repetitions per image/tree/mode: {metadata['configuration']['repetitions']}; mode order is deterministically balanced.",
        "- Resolution order is permuted by image and tree order is rotated by image/resolution to limit thermal-order confounding.",
        "- Ratios are paired by image. Confidence intervals use a deterministic image-level bootstrap.",
        "",
        f"Images are the experimental units (n={metadata['configuration']['image_count']}); individual tree nodes and repeated calls are not treated as independent replicates.",
        "",
    ]

    comparisons = (
        ("Complete groups", "dist_transf", "dist_transf_exact", "dist_transf_exact_over_dist_transf"),
        ("Maximum attributes", "max_dist", "max_dist_exact", "max_dist_exact_over_max_dist"),
    )
    for title, approximate_mode, exact_mode, ratio_prefix in comparisons:
        lines.extend([
            "",
            f"## {title}",
            "",
            "| Resolution | Tree | Approximate median (s) | Exact median (s) | Exact / approximate | 95% CI |",
            "| --- | --- | ---: | ---: | ---: | ---: |",
        ])
        for resolution in metadata["configuration"]["resolutions"]:
            for tree in metadata["configuration"]["tree_kinds"]:
                approximate = timing_lookup[(resolution, tree, approximate_mode)]
                exact = timing_lookup[(resolution, tree, exact_mode)]
                pair = pair_lookup[(resolution, tree)]
                lines.append(
                    f"| {resolution} | `{tree}` | {format_number(approximate['median_seconds'])} | "
                    f"{format_number(exact['median_seconds'])} | {format_number(pair[ratio_prefix + '_geomean'])}× | "
                    f"[{format_number(pair[ratio_prefix + '_ci95_low'])}, {format_number(pair[ratio_prefix + '_ci95_high'])}] |"
                )

    family_ratios = [float(row["dist_transf_exact_over_dist_transf_geomean"]) for row in pair_summary]
    maximum_ratios = [float(row["max_dist_exact_over_max_dist_geomean"]) for row in pair_summary]
    lines.extend([
        "",
        "## Ratio ranges",
        "",
        f"- `DIST_TRANSF_EXACT / DIST_TRANSF`: {min(family_ratios):.3f}× to {max(family_ratios):.3f}×.",
        f"- `MAX_DIST_EXACT / MAX_DIST`: {min(maximum_ratios):.3f}× to {max(maximum_ratios):.3f}×.",
    ])

    lines.extend([
        "",
        "## Resolution scaling",
        "",
        "| Tree | Mode | 720p / 480p | 1080p / 720p | 1080p / 480p | Scaling exponent |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ])
    for row in scaling:
        lines.append(
            f"| `{row['tree']}` | {row['mode']} | {format_number(row['time_ratio_720p_over_480p'])}× | "
            f"{format_number(row['time_ratio_1080p_over_720p'])}× | {format_number(row['time_ratio_1080p_over_480p'])}× | "
            f"{format_number(row['scaling_exponent_480p_to_1080p'])} |"
        )

    lines.extend([
        "",
        "## Tree construction outside the attribute timer",
        "",
        "| Resolution | Tree | Median construction (s) | Median nodes |",
        "| --- | --- | ---: | ---: |",
    ])
    for row in builds:
        lines.append(
            f"| {row['resolution']} | `{row['tree']}` | {format_number(row['median_build_seconds'])} | "
            f"{int(float(row['median_nodes'])):,} |"
        )

    lines.extend([
        "",
        "## Interpretation constraints",
        "",
        "- Timings include output allocation and Python-visible matrix materialization, but exclude resizing and tree construction.",
        "- Lanczos upsampling changes both pixel count and gray-level/tree complexity; resolution ratios are empirical, not pure asymptotic estimates.",
        "- Group and scalar comparisons are independent paired workloads; no combined exact/approximate call is timed.",
        "- Absolute times describe the recorded Release build and host. Paired ratios are more portable than raw seconds.",
        "",
        "## Artifacts",
        "",
        "- `experiment.json`: configuration, software, build, hardware, and input metadata.",
        "- `observations.csv`: raw timed calls in actual execution order.",
        "- `summary.csv`: descriptive time and throughput summaries.",
        "- `pairs.csv` and `pair-summary.csv`: paired exact/approximate ratios for complete groups and maximum attributes.",
        "- `scaling.csv`: empirical resolution ratios and log-log scaling exponent.",
        "- `tree-cases.csv` and `build-summary.csv`: node counts and untimed construction costs.",
        "",
    ])
    return "\n".join(lines)


def main() -> int:
    arguments = parse_args()
    if arguments.repetitions <= 0:
        raise ValueError("--repetitions must be positive")
    if arguments.bootstrap_replicates < 0:
        raise ValueError("--bootstrap-replicates must be non-negative")

    tree_kinds = parse_tree_kinds(arguments.trees)
    resolutions = parse_choices(arguments.resolutions, RESOLUTIONS, "--resolutions")
    modes = parse_choices(arguments.modes, MODES, "--modes")
    if set(modes) != set(MODES):
        raise ValueError("the analysis requires DIST_TRANSF, DIST_TRANSF_EXACT, MAX_DIST, and MAX_DIST_EXACT modes")
    images = discover_images(arguments.image_directory.resolve(), arguments.count)
    output_directory = arguments.output_directory.resolve()
    output_directory.mkdir(parents=True, exist_ok=True)
    observations_path = output_directory / "observations.csv"
    metadata_path = output_directory / "experiment.json"

    requests, attribute_names = discover_requests()
    repository = Path(__file__).resolve().parents[1]
    active_module = module_metadata()
    active_contract_mode = active_module["cmake"].get("MMCFILTERS_CONTRACT_MODE")
    if active_contract_mode not in {"CHECKED", "UNCHECKED"}:
        raise RuntimeError("the native module CMake cache does not identify a CHECKED or UNCHECKED contract mode")
    if arguments.expected_contract_mode is not None and active_contract_mode != arguments.expected_contract_mode:
        raise RuntimeError(
            f"expected {arguments.expected_contract_mode} native module, received {active_contract_mode}"
        )
    configuration = {
        "image_directory": str(arguments.image_directory.resolve()),
        "images": [path.name for path in images],
        "image_count": len(images),
        "selection_order": "lexicographic filename order",
        "tree_kinds": tree_kinds,
        "resolutions": resolutions,
        "resolution_shapes": {name: list(RESOLUTIONS[name]) for name in resolutions},
        "resize_filter": "Pillow LANCZOS",
        "modes": modes,
        "repetitions": arguments.repetitions,
        "radius": arguments.radius,
        "infinity_pixel": arguments.infinity_pixel,
        "dtype": "float32",
        "bootstrap_replicates": arguments.bootstrap_replicates,
        "contract_mode": active_contract_mode,
        "result_checksum": "BLAKE2b-128 over the complete layout and float32 result matrix",
        "attribute_names": attribute_names,
        "tos_default": "canonical Min4Max8, domain extension None, infinity pixel 0, UInt8",
        "execution_schedule": "resolution permutations by image; cyclic tree rotations by image/resolution; mode permutations by case/repetition",
    }
    if metadata_path.exists():
        previous = json.loads(metadata_path.read_text(encoding="utf-8"))
        if not arguments.resume:
            raise ValueError(f"output directory already contains an experiment; pass --resume: {output_directory}")
        if previous.get("configuration") != configuration:
            raise ValueError("--resume configuration does not match the existing experiment")
        if previous.get("module", {}).get("native_module_sha256") != active_module["native_module_sha256"]:
            raise ValueError("--resume native module does not match the existing experiment")
        metadata = previous
    else:
        metadata = {
            "schema_version": 3,
            "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "command": sys.argv,
            "configuration": configuration,
            "python": sys.version,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "cpu": cpu_name(),
            "numpy": np.__version__,
            "pillow": PILLOW_VERSION,
            "mmcfilters": getattr(mmcfilters, "__version__", "unknown"),
            "module": active_module,
            "git": git_metadata(repository),
        }
        metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    existing_observations = read_csv(observations_path)
    completed = {observation_key(row) for row in existing_observations}
    input_rows_by_key: dict[tuple[str, int], dict[str, Any]] = {}
    for row in read_csv(output_directory / "resized-inputs.csv"):
        input_rows_by_key[(row["resolution"], int(row["image_index"]))] = row
    tree_rows_by_key: dict[tuple[str, int, str], dict[str, Any]] = {}
    for row in read_csv(output_directory / "tree-cases.csv"):
        tree_rows_by_key[(row["resolution"], int(row["image_index"]), row["tree"])] = row

    warm_up(tree_kinds, modes, requests, arguments.radius, arguments.infinity_pixel)
    campaign_start = time.perf_counter()
    mode_orders = list(itertools.permutations(modes))
    total_tree_cases = len(resolutions) * len(images) * len(tree_kinds)
    tree_case_index = 0
    resolution_orders = list(itertools.permutations(resolutions))
    for image_index, path in enumerate(images):
        resolution_order = resolution_orders[image_index % len(resolution_orders)]
        for resolution in resolution_order:
            resolution_index = resolutions.index(resolution)
            rows, columns = RESOLUTIONS[resolution]
            image = resized_image(path, rows, columns)
            input_rows_by_key[(resolution, image_index)] = {
                "resolution": resolution,
                "image_index": image_index,
                "image": path.name,
                "source_sha256": file_sha256(path),
                "resized_sha256": array_sha256(image),
                "rows": rows,
                "columns": columns,
                "pixels": image.size,
            }
            write_csv(output_directory / "resized-inputs.csv", list(input_rows_by_key.values()))

            tree_shift = (image_index + resolution_index) % len(tree_kinds)
            tree_order = tree_kinds[tree_shift:] + tree_kinds[:tree_shift]
            for tree_kind in tree_order:
                tree_index = tree_kinds.index(tree_kind)
                tree_case_index += 1
                required_keys = {
                    (resolution, image_index, tree_kind, repetition, mode)
                    for repetition in range(arguments.repetitions)
                    for mode in modes
                }
                if required_keys <= completed:
                    print(f"[{tree_case_index:03d}/{total_tree_cases:03d}] resume-skip {resolution} {path.name} {tree_kind}", flush=True)
                    continue

                build_start = time.perf_counter_ns()
                tree = build_tree(image, tree_kind, arguments.radius, arguments.infinity_pixel)
                build_seconds = (time.perf_counter_ns() - build_start) / 1_000_000_000.0
                nodes = len(tree.alive_node_ids)
                for repetition in range(arguments.repetitions):
                    order_index = (resolution_index * len(images) * len(tree_kinds) + image_index * len(tree_kinds) + tree_index + repetition) % len(mode_orders)
                    order = mode_orders[order_index]
                    for order_position, mode in enumerate(order):
                        key = (resolution, image_index, tree_kind, repetition, mode)
                        if key in completed:
                            continue
                        gc.collect()
                        start = time.perf_counter_ns()
                        layout, values = mmcfilters.Attribute.compute_topology_attributes(
                            tree,
                            requests[mode],
                            dtype=np.float32,
                        )
                        elapsed_ns = time.perf_counter_ns() - start
                        checksum = result_checksum(mode, layout, values, attribute_names[mode])
                        seconds = elapsed_ns / 1_000_000_000.0
                        observation = {
                            "resolution": resolution,
                            "image_index": image_index,
                            "image": path.name,
                            "tree": tree_kind,
                            "repetition": repetition,
                            "order_position": order_position,
                            "mode": mode,
                            "rows": rows,
                            "columns": columns,
                            "pixels": image.size,
                            "nodes": nodes,
                            "seconds": seconds,
                            "milliseconds": elapsed_ns / 1_000_000.0,
                            "nanoseconds_per_pixel": elapsed_ns / image.size,
                            "nanoseconds_per_node": elapsed_ns / nodes,
                            "megapixels_per_second": (image.size / 1_000_000.0) / seconds,
                            "meganodes_per_second": (nodes / 1_000_000.0) / seconds,
                            "output_columns": len(layout),
                            "checksum": checksum,
                        }
                        append_observation(observations_path, observation)
                        completed.add(key)
                        del values, layout
                        print(
                            f"[{tree_case_index:03d}/{total_tree_cases:03d}] {resolution} {path.name} {tree_kind} "
                            f"rep={repetition + 1} mode={mode} nodes={nodes} time={seconds:.3f}s",
                            flush=True,
                        )

                tree_rows_by_key[(resolution, image_index, tree_kind)] = {
                    "resolution": resolution,
                    "image_index": image_index,
                    "image": path.name,
                    "tree": tree_kind,
                    "rows": rows,
                    "columns": columns,
                    "pixels": image.size,
                    "nodes": nodes,
                    "build_seconds": build_seconds,
                }
                write_csv(output_directory / "tree-cases.csv", list(tree_rows_by_key.values()))
                del tree

    observations = read_csv(observations_path)
    expected_observations = len(resolutions) * len(images) * len(tree_kinds) * arguments.repetitions * len(modes)
    if len(observations) != expected_observations:
        raise RuntimeError(f"expected {expected_observations} observations, found {len(observations)}")
    tree_rows = list(tree_rows_by_key.values())
    summary = summarize_observations(observations)
    builds = summarize_builds(tree_rows)
    raw_pairs, pair_summary = paired_rows(observations, arguments.bootstrap_replicates)
    scaling = scaling_rows(summary)
    write_csv(output_directory / "summary.csv", summary)
    write_csv(output_directory / "build-summary.csv", builds)
    write_csv(output_directory / "pairs.csv", raw_pairs)
    write_csv(output_directory / "pair-summary.csv", pair_summary)
    write_csv(
        output_directory / "scaling.csv",
        scaling,
        fieldnames=(
            "tree",
            "mode",
            "time_ratio_720p_over_480p",
            "time_ratio_1080p_over_720p",
            "time_ratio_1080p_over_480p",
            "pixel_ratio_720p_over_480p",
            "pixel_ratio_1080p_over_720p",
            "pixel_ratio_1080p_over_480p",
            "scaling_exponent_480p_to_1080p",
        ),
    )
    (output_directory / "analysis.md").write_text(
        build_report(metadata, summary, pair_summary, scaling, builds),
        encoding="utf-8",
    )
    metadata["completed_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    invocation_wall_seconds = time.perf_counter() - campaign_start
    previous_wall_seconds = float(
        metadata.get("campaign_wall_seconds_total", metadata.get("campaign_wall_seconds_this_invocation", 0.0))
    )
    metadata["campaign_wall_seconds_this_invocation"] = invocation_wall_seconds
    metadata["campaign_wall_seconds_total"] = previous_wall_seconds + invocation_wall_seconds
    metadata["observation_count"] = len(observations)
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"artifacts={output_directory}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
