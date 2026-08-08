#!/usr/bin/env python3
"""Validate and summarize the ICDAR tree-construction benchmark."""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path


METHODS = (
    "tos_max4c_min8c",
    "tos_self_dual",
    "max_tree_8c",
    "min_tree_8c",
    "residual_unrestricted_8c",
    "residual_saturated_8c",
)

METHOD_LABELS = {
    "tos_max4c_min8c": "ToS Max4cMin8c",
    "tos_self_dual": "ToS self-dual",
    "max_tree_8c": "Max-tree 8c",
    "min_tree_8c": "Min-tree 8c",
    "residual_unrestricted_8c": "Residual irrestrita 8c",
    "residual_saturated_8c": "Residual saturada 8c",
}

RESOLUTIONS = ("480p", "720p", "1080p")
EXPECTED_DIMENSIONS = {
    "480p": (480, 853, 409_440),
    "720p": (720, 1280, 921_600),
    "1080p": (1080, 1920, 2_073_600),
}
DEFAULT_PREFIX = "tree-construction-comparison-icdar-first10-2026-08-04"

PAIRWISE_COMPARISONS = (
    ("tos_self_dual", "tos_max4c_min8c"),
    ("min_tree_8c", "max_tree_8c"),
    ("residual_saturated_8c", "residual_unrestricted_8c"),
    ("residual_unrestricted_8c", "tos_self_dual"),
    ("residual_saturated_8c", "tos_self_dual"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("raw_csv", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--bootstrap-samples", type=int, default=20_000)
    parser.add_argument("--prefix", default=DEFAULT_PREFIX)
    return parser.parse_args()


def percentile(sorted_values: list[float], probability: float) -> float:
    position = probability * (len(sorted_values) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return sorted_values[lower]
    fraction = position - lower
    return (
        sorted_values[lower] * (1.0 - fraction)
        + sorted_values[upper] * fraction
    )


def bootstrap_mean_interval(
    values: list[float],
    samples: int,
    seed: int,
) -> tuple[float, float]:
    generator = random.Random(seed)
    size = len(values)
    means = [
        sum(values[generator.randrange(size)] for _ in range(size)) / size
        for _ in range(samples)
    ]
    means.sort()
    return percentile(means, 0.025), percentile(means, 0.975)


def geometric_mean(values: list[float]) -> float:
    return math.exp(sum(math.log(value) for value in values) / len(values))


def regression_exponent(points: list[tuple[int, float]]) -> float:
    x_values = [math.log(pixels) for pixels, _ in points]
    y_values = [math.log(milliseconds) for _, milliseconds in points]
    x_mean = statistics.mean(x_values)
    y_mean = statistics.mean(y_values)
    numerator = sum(
        (x_value - x_mean) * (y_value - y_mean)
        for x_value, y_value in zip(x_values, y_values)
    )
    denominator = sum((x_value - x_mean) ** 2 for x_value in x_values)
    return numerator / denominator


def load_and_validate(path: Path) -> list[dict[str, object]]:
    with path.open(newline="", encoding="utf-8") as stream:
        raw_rows = list(csv.DictReader(stream))

    expected_rows = len(RESOLUTIONS) * 10 * len(METHODS) * 5
    if len(raw_rows) != expected_rows:
        raise ValueError(f"expected {expected_rows} rows, found {len(raw_rows)}")

    rows: list[dict[str, object]] = []
    groups: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    schedules: dict[tuple[str, str, int], list[dict[str, object]]] = defaultdict(list)
    for raw in raw_rows:
        row: dict[str, object] = {
            "resolution": raw["resolution"],
            "image": raw["image"],
            "rows": int(raw["rows"]),
            "cols": int(raw["cols"]),
            "pixels": int(raw["pixels"]),
            "repetition": int(raw["repetition"]),
            "position": int(raw["position"]),
            "algorithm": raw["algorithm"],
            "construction_ms": float(raw["construction_ms"]),
            "nodes": int(raw["nodes"]),
        }
        resolution = str(row["resolution"])
        algorithm = str(row["algorithm"])
        if resolution not in RESOLUTIONS:
            raise ValueError(f"unexpected resolution: {resolution}")
        if algorithm not in METHODS:
            raise ValueError(f"unexpected algorithm: {algorithm}")
        dimensions = (row["rows"], row["cols"], row["pixels"])
        if dimensions != EXPECTED_DIMENSIONS[resolution]:
            raise ValueError(f"unexpected dimensions for {resolution}: {dimensions}")
        if float(row["construction_ms"]) <= 0.0:
            raise ValueError("construction times must be positive")
        groups[(resolution, str(row["image"]), algorithm)].append(row)
        schedules[(resolution, str(row["image"]), int(row["repetition"]))].append(row)
        rows.append(row)

    expected_groups = len(RESOLUTIONS) * 10 * len(METHODS)
    if len(groups) != expected_groups:
        raise ValueError(f"expected {expected_groups} groups, found {len(groups)}")
    for key, group in groups.items():
        if len(group) != 5:
            raise ValueError(f"expected five repetitions for {key}, found {len(group)}")
        if {int(row["repetition"]) for row in group} != set(range(5)):
            raise ValueError(f"invalid repetitions for {key}")
        if len({int(row["nodes"]) for row in group}) != 1:
            raise ValueError(f"node count changed between repetitions for {key}")

    for key, schedule in schedules.items():
        if {str(row["algorithm"]) for row in schedule} != set(METHODS):
            raise ValueError(f"invalid algorithm schedule for {key}")
        if {int(row["position"]) for row in schedule} != set(range(6)):
            raise ValueError(f"invalid position schedule for {key}")
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"cannot write an empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = load_and_validate(args.raw_csv)

    groups: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        groups[(
            str(row["resolution"]),
            str(row["image"]),
            str(row["algorithm"]),
        )].append(row)

    per_image: list[dict[str, object]] = []
    per_image_lookup: dict[tuple[str, str, str], float] = {}
    for resolution in RESOLUTIONS:
        for image_index in range(10):
            image = f"test_{image_index:03d}.png"
            for method in METHODS:
                group = groups[(resolution, image, method)]
                times = [float(row["construction_ms"]) for row in group]
                median_ms = statistics.median(times)
                relative_range = (max(times) - min(times)) / median_ms
                per_image_lookup[(resolution, image, method)] = median_ms
                per_image.append({
                    "resolution": resolution,
                    "image": image,
                    "algorithm": method,
                    "algorithm_label": METHOD_LABELS[method],
                    "median_ms": round(median_ms, 6),
                    "min_ms": round(min(times), 6),
                    "max_ms": round(max(times), 6),
                    "relative_range_pct": round(relative_range * 100.0, 4),
                    "nodes": int(group[0]["nodes"]),
                })

    aggregates: list[dict[str, object]] = []
    aggregate_lookup: dict[tuple[str, str], dict[str, object]] = {}
    for resolution_index, resolution in enumerate(RESOLUTIONS):
        pixels = EXPECTED_DIMENSIONS[resolution][2]
        resolution_rows: list[dict[str, object]] = []
        for method_index, method in enumerate(METHODS):
            values = [
                per_image_lookup[(resolution, f"test_{index:03d}.png", method)]
                for index in range(10)
            ]
            nodes = [
                int(groups[(resolution, f"test_{index:03d}.png", method)][0]["nodes"])
                for index in range(10)
            ]
            relative_ranges = [
                float(row["relative_range_pct"])
                for row in per_image
                if row["resolution"] == resolution
                and row["algorithm"] == method
            ]
            ci_low, ci_high = bootstrap_mean_interval(
                values,
                args.bootstrap_samples,
                20260804 + resolution_index * 100 + method_index,
            )
            result: dict[str, object] = {
                "resolution": resolution,
                "pixels": pixels,
                "megapixels": round(pixels / 1_000_000.0, 6),
                "algorithm": method,
                "algorithm_label": METHOD_LABELS[method],
                "mean_ms": round(statistics.mean(values), 6),
                "mean_ci95_low_ms": round(ci_low, 6),
                "mean_ci95_high_ms": round(ci_high, 6),
                "median_ms": round(statistics.median(values), 6),
                "standard_deviation_ms": round(statistics.stdev(values), 6),
                "mean_ms_per_megapixel": round(
                    statistics.mean(values) / (pixels / 1_000_000.0), 6
                ),
                "mean_nodes": round(statistics.mean(nodes), 1),
                "median_repetition_relative_range_pct": round(
                    statistics.median(relative_ranges), 4
                ),
            }
            resolution_rows.append(result)
            aggregate_lookup[(resolution, method)] = result
        fastest = min(float(row["mean_ms"]) for row in resolution_rows)
        for rank, result in enumerate(
            sorted(resolution_rows, key=lambda row: float(row["mean_ms"])),
            start=1,
        ):
            result["rank"] = rank
            result["relative_to_fastest"] = round(
                float(result["mean_ms"]) / fastest, 6
            )
        aggregates.extend(resolution_rows)

    scaling: list[dict[str, object]] = []
    for method in METHODS:
        points = [
            (
                EXPECTED_DIMENSIONS[resolution][2],
                float(aggregate_lookup[(resolution, method)]["mean_ms"]),
            )
            for resolution in RESOLUTIONS
        ]
        scaling.append({
            "algorithm": method,
            "algorithm_label": METHOD_LABELS[method],
            "empirical_exponent": round(regression_exponent(points), 6),
            "time_ratio_1080p_over_480p": round(points[-1][1] / points[0][1], 6),
            "pixel_ratio_1080p_over_480p": round(points[-1][0] / points[0][0], 6),
        })

    pairwise: list[dict[str, object]] = []
    for resolution in RESOLUTIONS:
        for numerator, denominator in PAIRWISE_COMPARISONS:
            ratios = [
                per_image_lookup[(resolution, f"test_{index:03d}.png", numerator)]
                / per_image_lookup[(resolution, f"test_{index:03d}.png", denominator)]
                for index in range(10)
            ]
            pairwise.append({
                "resolution": resolution,
                "comparison": f"{numerator}/{denominator}",
                "numerator_label": METHOD_LABELS[numerator],
                "denominator_label": METHOD_LABELS[denominator],
                "ratio_of_means": round(
                    float(aggregate_lookup[(resolution, numerator)]["mean_ms"])
                    / float(aggregate_lookup[(resolution, denominator)]["mean_ms"]),
                    6,
                ),
                "geometric_mean_paired_ratio": round(geometric_mean(ratios), 6),
                "median_paired_ratio": round(statistics.median(ratios), 6),
                "min_paired_ratio": round(min(ratios), 6),
                "max_paired_ratio": round(max(ratios), 6),
            })

    prefix = args.prefix
    write_csv(args.output_dir / f"{prefix}-per-image.csv", per_image)
    write_csv(args.output_dir / f"{prefix}-summary.csv", aggregates)
    write_csv(args.output_dir / f"{prefix}-scaling.csv", scaling)
    write_csv(args.output_dir / f"{prefix}-paired-ratios.csv", pairwise)
    validation = {
        "status": "ready_to_share",
        "raw_rows": len(rows),
        "groups": len(groups),
        "repetitions_per_group": 5,
        "images_per_resolution": 10,
        "resolutions": list(RESOLUTIONS),
        "algorithms": list(METHODS),
        "exact_reconstruction_checked_during_warmup": True,
        "deterministic_node_counts": True,
        "positive_times": True,
        "complete_rotated_schedules": True,
        "measurement_scope": "external_complete_public_construction_call",
        "production_internal_timing": False,
        "bootstrap_samples": args.bootstrap_samples,
        "bootstrap_seed_base": 20260804,
    }
    with (args.output_dir / f"{prefix}-validation.json").open(
        "w", encoding="utf-8"
    ) as stream:
        json.dump(validation, stream, indent=2)
        stream.write("\n")


if __name__ == "__main__":
    main()
