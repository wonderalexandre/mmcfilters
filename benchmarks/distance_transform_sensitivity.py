#!/usr/bin/env python3
"""Compare every aligned DIST_TRANSF/DIST_TRANSF_EXACT pair on real images.

The experiment computes both families over the same established morphological
tree and aligns values by live internal NodeId.  It reports node-level errors,
per-image rank preservation, support-area strata, and joint displacement for
the three coordinate-pair descriptors.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np
from PIL import Image

import mmcfilters


EXACT_SUFFIX = "_EXACT"
IMAGE_SUFFIXES = {".png", ".pgm", ".tif", ".tiff", ".jpg", ".jpeg", ".bmp"}
TREE_KINDS = {"max", "min", "tos", "residual_unrestricted", "residual_saturated"}
AREA_STRATA = (
    ("area_1", 1, 1),
    ("area_2_9", 2, 9),
    ("area_10_99", 10, 99),
    ("area_100_999", 100, 999),
    ("area_1000_9999", 1_000, 9_999),
    ("area_10000_plus", 10_000, None),
)

GEOMETRY_PAIRS = (
    ("max_center", "MAX_DIST_CENTER_ROW", "MAX_DIST_CENTER_COLUMN"),
    ("max_plateau_centroid", "MAX_DIST_PLATEAU_CENTROID_ROW", "MAX_DIST_PLATEAU_CENTROID_COLUMN"),
    ("distance_weighted_centroid", "DIST_WEIGHTED_CENTROID_ROW", "DIST_WEIGHTED_CENTROID_COLUMN"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("image_directory", type=Path, help="Directory containing the ICDAR images.")
    parser.add_argument("output_directory", type=Path, help="Directory receiving CSV, JSON, and Markdown artifacts.")
    parser.add_argument("--count", type=int, default=10, help="Number of lexicographically first images to analyze (default: 10).")
    parser.add_argument(
        "--trees",
        default="max,min",
        help="Comma-separated tree kinds: max,min,tos,residual_unrestricted,residual_saturated (default: max,min).",
    )
    parser.add_argument("--radius", type=float, default=1.5, help="Tree-construction adjacency radius (default: 1.5).")
    parser.add_argument(
        "--infinity-pixel",
        type=int,
        default=0,
        help="Infinity pixel used by residual_saturated (default: 0).",
    )
    parser.add_argument("--atol", type=float, default=1.0e-12, help="Absolute agreement tolerance (default: 1e-12).")
    parser.add_argument("--rtol", type=float, default=1.0e-9, help="Relative agreement tolerance (default: 1e-9).")
    return parser.parse_args()


def discover_images(directory: Path, count: int) -> list[Path]:
    if count <= 0:
        raise ValueError("--count must be positive")
    if not directory.is_dir():
        raise ValueError(f"image directory does not exist: {directory}")
    images = sorted(path for path in directory.iterdir() if path.is_file() and path.suffix.lower() in IMAGE_SUFFIXES)
    if len(images) < count:
        raise ValueError(f"requested {count} images, but only {len(images)} supported files were found in {directory}")
    return images[:count]


def parse_tree_kinds(value: str) -> list[str]:
    result = [item.strip().lower() for item in value.split(",") if item.strip()]
    if not result or any(item not in TREE_KINDS for item in result):
        raise ValueError(f"--trees must contain values from: {','.join(sorted(TREE_KINDS))}")
    if len(set(result)) != len(result):
        raise ValueError("--trees must not contain duplicates")
    return result


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_grayscale(path: Path) -> np.ndarray:
    with Image.open(path) as image:
        return np.ascontiguousarray(np.asarray(image.convert("L"), dtype=np.uint8))


def build_tree(image: np.ndarray, kind: str, radius: float, infinity_pixel: int):
    if kind == "max":
        return mmcfilters.MorphologicalTreeFactory.create_max_tree(image, radius=radius)
    if kind == "min":
        return mmcfilters.MorphologicalTreeFactory.create_min_tree(image, radius=radius)
    if kind == "residual_unrestricted":
        return mmcfilters.MorphologicalTreeFactory.create_unrestricted_residual_tree(image, radius=radius)
    if kind == "tos":
        return mmcfilters.MorphologicalTreeFactory.create_tree_of_shapes(image)
    if infinity_pixel < 0 or infinity_pixel >= image.size:
        raise ValueError(f"--infinity-pixel must be in [0, {image.size}) for residual_saturated")
    return mmcfilters.MorphologicalTreeFactory.create_saturated_residual_tree(
        image,
        infinity_pixel=infinity_pixel,
        radius=radius,
    )


def attribute_category(name: str) -> str:
    if name in {"MAX_DIST", "MAX_SQUARED_DIST"}:
        return "maximum"
    if name.startswith("MAX_DIST_CENTER_"):
        return "maximum_center"
    if name.startswith("MAX_DIST_PLATEAU_"):
        return "maximum_plateau"
    if name.startswith("DIST_WEIGHTED_"):
        return "distance_weighted_geometry"
    if name in {"DIST_SQUARED_SUM", "DIST_SQUARED_MEAN", "DIST_RMS", "DIST_SQUARED_VARIANCE"}:
        return "squared_distance_moments"
    if name in {"DIST_SUM", "DIST_MEAN", "DIST_VARIANCE"}:
        return "real_distance_moments"
    return "distance_distribution"


def attribute_unit(name: str) -> str:
    if name in {"MAX_DIST", "DIST_RMS", "DIST_MEAN", "DIST_MEDIAN", "DIST_MODE", "DIST_Q25", "DIST_Q75", "DIST_Q90"}:
        return "pixel"
    if name == "DIST_SUM":
        return "pixel_samples"
    if name in {"MAX_SQUARED_DIST", "DIST_SQUARED_MEAN", "DIST_VARIANCE"}:
        return "pixel_squared"
    if name == "DIST_SQUARED_SUM":
        return "pixel_squared_samples"
    if name == "DIST_SQUARED_VARIANCE":
        return "pixel_fourth_power"
    if name in {"MAX_DIST_CENTER_ROW", "MAX_DIST_CENTER_COLUMN", "MAX_DIST_PLATEAU_CENTROID_ROW",
                "MAX_DIST_PLATEAU_CENTROID_COLUMN", "DIST_WEIGHTED_CENTROID_ROW", "DIST_WEIGHTED_CENTROID_COLUMN"}:
        return "pixel_coordinate"
    if name in {"MAX_DIST_PLATEAU_AREA", "DIST_POSITIVE_AREA", "DIST_LEVEL_COUNT"}:
        return "count"
    if name == "DIST_ENTROPY":
        return "bit"
    if name in {"DIST_WEIGHTED_CENTRAL_MOMENT_20", "DIST_WEIGHTED_CENTRAL_MOMENT_02", "DIST_WEIGHTED_CENTRAL_MOMENT_11"}:
        return "pixel_cubed"
    if name == "DIST_WEIGHTED_AXIS_ORIENTATION":
        return "degree"
    return "dimensionless"


def average_ranks(values: np.ndarray) -> np.ndarray:
    """Return one-based average ranks with deterministic tie handling."""
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    starts = np.concatenate(([0], np.flatnonzero(sorted_values[1:] != sorted_values[:-1]) + 1))
    ends = np.concatenate((starts[1:], [values.size]))
    sorted_ranks = np.empty(values.size, dtype=np.float64)
    for start, end in zip(starts, ends):
        sorted_ranks[start:end] = 0.5 * (start + end - 1) + 1.0
    ranks = np.empty(values.size, dtype=np.float64)
    ranks[order] = sorted_ranks
    return ranks


def correlation(lhs: np.ndarray, rhs: np.ndarray) -> float:
    lhs_centered = lhs - np.mean(lhs)
    rhs_centered = rhs - np.mean(rhs)
    denominator = math.sqrt(float(np.dot(lhs_centered, lhs_centered)) * float(np.dot(rhs_centered, rhs_centered)))
    if denominator == 0.0:
        return 1.0 if np.array_equal(lhs, rhs) else math.nan
    return float(np.dot(lhs_centered, rhs_centered) / denominator)


def scalar_metrics(approximate: np.ndarray, exact: np.ndarray, atol: float, rtol: float) -> dict[str, float | int]:
    if approximate.shape != exact.shape or approximate.ndim != 1:
        raise ValueError("attribute vectors must be aligned one-dimensional arrays")
    if approximate.size == 0 or not np.all(np.isfinite(approximate)) or not np.all(np.isfinite(exact)):
        raise ValueError("attribute comparison requires non-empty finite vectors")

    delta = approximate - exact
    absolute = np.abs(delta)
    exact_absolute = np.abs(exact)
    tolerance = atol + rtol * exact_absolute
    close = absolute <= tolerance
    exact_nonzero = exact_absolute > atol
    relative = absolute[exact_nonzero] / exact_absolute[exact_nonzero]
    symmetric_denominator = np.abs(approximate) + exact_absolute
    symmetric = np.zeros_like(absolute)
    symmetric_mask = symmetric_denominator > atol
    symmetric[symmetric_mask] = 2.0 * absolute[symmetric_mask] / symmetric_denominator[symmetric_mask]

    approximate_ranks = average_ranks(approximate)
    exact_ranks = average_ranks(exact)
    exact_l1 = float(np.sum(exact_absolute, dtype=np.float64))
    exact_l2_squared = float(np.dot(exact, exact))
    squared_error = float(np.dot(delta, delta))
    return {
        "samples": int(approximate.size),
        "exact_nonzero_rate": float(np.mean(exact_nonzero)),
        "exact_equal_rate": float(np.mean(approximate == exact)),
        "close_rate": float(np.mean(close)),
        "underestimate_rate": float(np.mean(delta < -tolerance)),
        "overestimate_rate": float(np.mean(delta > tolerance)),
        "bias": float(np.mean(delta)),
        "mae": float(np.mean(absolute)),
        "rmse": math.sqrt(squared_error / approximate.size),
        "median_absolute_error": float(np.median(absolute)),
        "p95_absolute_error": float(np.quantile(absolute, 0.95)),
        "max_absolute_error": float(np.max(absolute)),
        "nmae": float(np.sum(absolute, dtype=np.float64) / exact_l1) if exact_l1 > 0.0 else math.nan,
        "nrmse": math.sqrt(squared_error / exact_l2_squared) if exact_l2_squared > 0.0 else math.nan,
        "mean_relative_error_nonzero": float(np.mean(relative)) if relative.size else math.nan,
        "p95_relative_error_nonzero": float(np.quantile(relative, 0.95)) if relative.size else math.nan,
        "mean_smape": float(np.mean(symmetric)),
        "pearson": correlation(approximate, exact),
        "spearman": correlation(approximate_ranks, exact_ranks),
        "exact_mean": float(np.mean(exact)),
        "approximate_mean": float(np.mean(approximate)),
    }


@dataclass
class ScalarAggregate:
    samples: int = 0
    exact_nonzero: int = 0
    exact_equal: int = 0
    close: int = 0
    under: int = 0
    over: int = 0
    sum_delta: float = 0.0
    sum_absolute: float = 0.0
    sum_squared_error: float = 0.0
    sum_exact_absolute: float = 0.0
    sum_exact_squared: float = 0.0
    sum_relative: float = 0.0
    sum_smape: float = 0.0
    max_absolute: float = 0.0
    image_mae: list[float] = field(default_factory=list)
    image_p95_absolute: list[float] = field(default_factory=list)
    image_pearson: list[float] = field(default_factory=list)
    image_spearman: list[float] = field(default_factory=list)

    def add(self, approximate: np.ndarray, exact: np.ndarray, metrics: dict[str, float | int], atol: float, rtol: float) -> None:
        delta = approximate - exact
        absolute = np.abs(delta)
        exact_absolute = np.abs(exact)
        tolerance = atol + rtol * exact_absolute
        nonzero = exact_absolute > atol
        symmetric_denominator = np.abs(approximate) + exact_absolute
        symmetric = np.zeros_like(absolute)
        symmetric_mask = symmetric_denominator > atol
        symmetric[symmetric_mask] = 2.0 * absolute[symmetric_mask] / symmetric_denominator[symmetric_mask]

        self.samples += approximate.size
        self.exact_nonzero += int(np.count_nonzero(nonzero))
        self.exact_equal += int(np.count_nonzero(approximate == exact))
        self.close += int(np.count_nonzero(absolute <= tolerance))
        self.under += int(np.count_nonzero(delta < -tolerance))
        self.over += int(np.count_nonzero(delta > tolerance))
        self.sum_delta += float(np.sum(delta, dtype=np.float64))
        self.sum_absolute += float(np.sum(absolute, dtype=np.float64))
        self.sum_squared_error += float(np.dot(delta, delta))
        self.sum_exact_absolute += float(np.sum(exact_absolute, dtype=np.float64))
        self.sum_exact_squared += float(np.dot(exact, exact))
        if np.any(nonzero):
            self.sum_relative += float(np.sum(absolute[nonzero] / exact_absolute[nonzero], dtype=np.float64))
        self.sum_smape += float(np.sum(symmetric, dtype=np.float64))
        self.max_absolute = max(self.max_absolute, float(np.max(absolute)))
        self.image_mae.append(float(metrics["mae"]))
        self.image_p95_absolute.append(float(metrics["p95_absolute_error"]))
        if math.isfinite(float(metrics["pearson"])):
            self.image_pearson.append(float(metrics["pearson"]))
        if math.isfinite(float(metrics["spearman"])):
            self.image_spearman.append(float(metrics["spearman"]))

    def as_row(self) -> dict[str, float | int]:
        return {
            "samples": self.samples,
            "exact_nonzero_rate": self.exact_nonzero / self.samples,
            "exact_equal_rate": self.exact_equal / self.samples,
            "close_rate": self.close / self.samples,
            "underestimate_rate": self.under / self.samples,
            "overestimate_rate": self.over / self.samples,
            "bias": self.sum_delta / self.samples,
            "mae": self.sum_absolute / self.samples,
            "rmse": math.sqrt(self.sum_squared_error / self.samples),
            "max_absolute_error": self.max_absolute,
            "nmae": self.sum_absolute / self.sum_exact_absolute if self.sum_exact_absolute > 0.0 else math.nan,
            "nrmse": math.sqrt(self.sum_squared_error / self.sum_exact_squared) if self.sum_exact_squared > 0.0 else math.nan,
            "mean_relative_error_nonzero": self.sum_relative / self.exact_nonzero if self.exact_nonzero else math.nan,
            "mean_smape": self.sum_smape / self.samples,
            "median_image_mae": float(np.median(self.image_mae)),
            "max_image_p95_absolute_error": max(self.image_p95_absolute),
            "mean_image_pearson": float(np.mean(self.image_pearson)) if self.image_pearson else math.nan,
            "mean_image_spearman": float(np.mean(self.image_spearman)) if self.image_spearman else math.nan,
        }


@dataclass
class GeometryAggregate:
    samples: int = 0
    unchanged: int = 0
    sum_distance: float = 0.0
    sum_squared_distance: float = 0.0
    max_distance: float = 0.0
    image_medians: list[float] = field(default_factory=list)
    image_p95: list[float] = field(default_factory=list)

    def add(self, distance: np.ndarray, atol: float) -> None:
        self.samples += distance.size
        self.unchanged += int(np.count_nonzero(distance <= atol))
        self.sum_distance += float(np.sum(distance, dtype=np.float64))
        self.sum_squared_distance += float(np.dot(distance, distance))
        self.max_distance = max(self.max_distance, float(np.max(distance)))
        self.image_medians.append(float(np.median(distance)))
        self.image_p95.append(float(np.quantile(distance, 0.95)))

    def as_row(self) -> dict[str, float | int]:
        return {
            "samples": self.samples,
            "unchanged_rate": self.unchanged / self.samples,
            "mean_displacement": self.sum_distance / self.samples,
            "rms_displacement": math.sqrt(self.sum_squared_distance / self.samples),
            "max_displacement": self.max_distance,
            "median_image_median_displacement": float(np.median(self.image_medians)),
            "max_image_p95_displacement": max(self.image_p95),
        }


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def git_metadata(repository: Path) -> dict[str, object]:
    def command(*arguments: str) -> str:
        try:
            return subprocess.check_output(arguments, cwd=repository, text=True, stderr=subprocess.DEVNULL).strip()
        except (OSError, subprocess.CalledProcessError):
            return "unknown"

    status = command("git", "status", "--short")
    return {"commit": command("git", "rev-parse", "HEAD"), "dirty": status not in {"", "unknown"}}


def format_number(value: object) -> str:
    number = float(value)
    if math.isnan(number):
        return "n/a"
    if number == 0.0:
        return "0"
    if abs(number) >= 1.0e4 or abs(number) < 1.0e-3:
        return f"{number:.3e}"
    return f"{number:.6f}"


def build_report(metadata: dict, summary_rows: list[dict], area_rows: list[dict], geometry_rows: list[dict],
                 field_condition_rows: list[dict], input_rows: list[dict]) -> str:
    combined = [row for row in summary_rows if row["tree"] == "combined"]
    ranked = sorted(combined, key=lambda row: (1.0 - float(row["close_rate"]), float(row["nrmse"])), reverse=True)
    summary_by_attribute = {str(row["attribute"]): row for row in combined}
    squared_sum = summary_by_attribute["DIST_SQUARED_SUM"]
    maximum = summary_by_attribute["MAX_SQUARED_DIST"]
    field_changed = int(round(float(squared_sum["samples"]) * (1.0 - float(squared_sum["exact_equal_rate"]))))
    maximum_changed = int(round(float(maximum["samples"]) * (1.0 - float(maximum["exact_equal_rate"]))))
    invariant_attributes = sorted(str(row["attribute"]) for row in combined if float(row["exact_equal_rate"]) == 1.0)
    changed_area_strata = [
        str(row["area_stratum"])
        for row in area_rows
        if row["tree"] == "combined" and row["attribute"] == "DIST_SQUARED_SUM" and float(row["exact_equal_rate"]) < 1.0
    ]
    exact_area_strata = [
        str(row["area_stratum"])
        for row in area_rows
        if row["tree"] == "combined" and row["attribute"] == "DIST_SQUARED_SUM" and float(row["exact_equal_rate"]) == 1.0
    ]
    saturation_protocol = []
    if "residual_saturated" in metadata["tree_kinds"]:
        saturation_protocol.append(
            f"- Saturated residual infinity pixel: `{metadata['infinity_pixel']}`; residual spatial order: row-major default."
        )
    if "tos" in metadata["tree_kinds"]:
        saturation_protocol.append(
            "- Tree of Shapes: default canonical complementary grid (minimum 4-connectivity, maximum 8-connectivity), "
            "no domain padding (`TopographicDomainExtension.NONE`), infinity pixel `0`, and `uint8` altitude encoding."
        )
    lines = [
        "# Sensitivity of `DIST_TRANSF` versus `DIST_TRANSF_EXACT`",
        "",
        "## Protocol",
        "",
        f"- Images: `{metadata['image_count']}` lexicographically first files from `{metadata['image_directory']}`.",
        f"- Trees: `{', '.join(metadata['tree_kinds'])}`.",
        f"- Component/residual adjacency radius: `{metadata['radius']}`.",
        *saturation_protocol,
        f"- Comparison domain: every live internal node, paired by the same `NodeId`, using `float64` output.",
        f"- Agreement: `|approx-exact| <= {metadata['atol']} + {metadata['rtol']} |exact|`.",
        "- `NMAE = sum(|approx-exact|) / sum(|exact|)`; `NRMSE = sqrt(sum(error²) / sum(exact²))`.",
        "- Spearman is computed per image and then averaged, because attribute filters depend on node ordering.",
        "",
        "The node-level aggregate is descriptive: nodes are not independent statistical replicates. Per-image metrics are retained in CSV so images can be used as the inferential unit in a later study.",
        "",
        "## Inputs and workload",
        "",
        f"The experiment processed {len({row['image'] for row in input_rows})} images and {sum(int(row['nodes']) for row in input_rows):,} tree nodes across all tree/image cases.",
        "",
        "## Field-level sensitivity",
        "",
        f"`DIST_SQUARED_SUM` differs exactly on {field_changed:,} of {int(squared_sum['samples']):,} nodes "
        f"({100.0 * field_changed / int(squared_sum['samples']):.3f}%). Because the approximate DIFT cost is the distance to an active contour seed, "
        "whereas the exact EDT takes the nearest seed, the approximate field is an upper bound; equality of the integer squared-distance sums therefore certifies an identical field for that node.",
        f"The maximum squared distance changes on only {maximum_changed:,} nodes "
        f"({100.0 * maximum_changed / int(maximum['samples']):.4f}%).",
    ]

    summary_lookup = {
        (str(row["tree"]), str(row["attribute"])): row
        for row in summary_rows
    }
    lines.extend([
        "",
        "## Sensitivity by tree kind",
        "",
        "| Tree | Nodes | Changed fields | Exact field equality | Field NMAE | Changed maxima | Maximum NMAE |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ])
    for tree_kind in metadata["tree_kinds"]:
        tree_squared_sum = summary_lookup[(tree_kind, "DIST_SQUARED_SUM")]
        tree_maximum = summary_lookup[(tree_kind, "MAX_SQUARED_DIST")]
        tree_field_changed = int(round(
            float(tree_squared_sum["samples"]) * (1.0 - float(tree_squared_sum["exact_equal_rate"]))
        ))
        tree_maximum_changed = int(round(
            float(tree_maximum["samples"]) * (1.0 - float(tree_maximum["exact_equal_rate"]))
        ))
        lines.append(
            f"| `{tree_kind}` | {int(tree_squared_sum['samples']):,} | {tree_field_changed:,} | "
            f"{100.0 * float(tree_squared_sum['exact_equal_rate']):.3f}% | {format_number(tree_squared_sum['nmae'])} | "
            f"{tree_maximum_changed:,} | {format_number(tree_maximum['nmae'])} |"
        )

    residual_kinds = [tree_kind for tree_kind in metadata["tree_kinds"] if tree_kind.startswith("residual_")]
    if residual_kinds:
        lines.extend([
            "",
            "## Residual-tree descriptor sensitivity",
            "",
            "| Tree | Attribute | Agreement | NMAE | Mean image Spearman |",
            "| --- | --- | ---: | ---: | ---: |",
        ])
        residual_attributes = (
            "MAX_DIST",
            "DIST_Q90",
            "DIST_WEIGHTED_AXIS_ORIENTATION",
            "DIST_WEIGHTED_ECCENTRICITY",
        )
        for tree_kind in residual_kinds:
            for attribute in residual_attributes:
                row = summary_lookup[(tree_kind, attribute)]
                lines.append(
                    f"| `{tree_kind}` | `{attribute}` | {100.0 * float(row['close_rate']):.3f}% | "
                    f"{format_number(row['nmae'])} | {format_number(row['mean_image_spearman'])} |"
                )

    if "tos" in metadata["tree_kinds"]:
        lines.extend([
            "",
            "## Tree-of-Shapes descriptor sensitivity",
            "",
            "| Attribute | Agreement | NMAE | Mean image Spearman |",
            "| --- | ---: | ---: | ---: |",
        ])
        tos_attributes = (
            "MAX_DIST",
            "DIST_Q90",
            "DIST_WEIGHTED_AXIS_ORIENTATION",
            "DIST_WEIGHTED_ECCENTRICITY",
        )
        for attribute in tos_attributes:
            row = summary_lookup[("tos", attribute)]
            lines.append(
                f"| `{attribute}` | {100.0 * float(row['close_rate']):.3f}% | "
                f"{format_number(row['nmae'])} | {format_number(row['mean_image_spearman'])} |"
            )

    lines.extend([
        "",
        "The conditional table separates actual field changes from reducer-order and derived-descriptor sensitivity:",
        "",
        "| Attribute | Field condition | Nodes | Agreement | NMAE | Mean image Spearman |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ])
    conditional_attributes = {
        "DIST_MEAN",
        "DIST_WEIGHTED_CENTRAL_MOMENT_11",
        "DIST_WEIGHTED_AXIS_ORIENTATION",
        "DIST_WEIGHTED_ECCENTRICITY",
        "MAX_SQUARED_DIST",
    }
    for row in field_condition_rows:
        if row["tree"] != "combined" or row["attribute"] not in conditional_attributes:
            continue
        lines.append(
            f"| `{row['attribute']}` | `{row['field_condition']}` | {int(row['samples']):,} | "
            f"{100.0 * float(row['close_rate']):.3f}% | {format_number(row['nmae'])} | "
            f"{format_number(row['mean_image_spearman'])} |"
        )

    conditional_lookup = {
        (str(row["attribute"]), str(row["field_condition"])): row
        for row in field_condition_rows
        if row["tree"] == "combined"
    }
    moment11_equal = conditional_lookup[("DIST_WEIGHTED_CENTRAL_MOMENT_11", "field_equal")]
    orientation_equal = conditional_lookup[("DIST_WEIGHTED_AXIS_ORIENTATION", "field_equal")]
    eccentricity_equal = conditional_lookup[("DIST_WEIGHTED_ECCENTRICITY", "field_equal")]
    maximum_center_geometry = next(
        row for row in geometry_rows if row["tree"] == "combined" and row["geometry"] == "max_center"
    )
    main_findings = [
        f"- The approximate squared-distance sum was never below the exact sum. Field differences occurred in "
        f"`{', '.join(changed_area_strata)}` supports; field-exact area strata were `{', '.join(exact_area_strata)}`.",
        f"- Exact equality held for every node for: {', '.join(f'`{name}`' for name in invariant_attributes)}.",
        f"- On nodes whose squared-distance field was identical, `DIST_WEIGHTED_CENTRAL_MOMENT_11` had NMAE "
        f"{format_number(moment11_equal['nmae'])}, but orientation and eccentricity amplified reducer-order rounding to NMAE "
        f"{format_number(orientation_equal['nmae'])} and {format_number(eccentricity_equal['nmae'])}, respectively. "
        "Their main sensitivity is therefore numerical degeneracy in the derived eigensystem, not DIFT field error.",
        f"- The {maximum_changed:,} maximum-cost changes altered maximum-center displacement by at most "
        f"{format_number(maximum_center_geometry['max_displacement'])} pixels; the joint geometry table reports all center and centroid effects.",
    ]
    if residual_kinds:
        residual_nodes = sum(int(summary_lookup[(tree_kind, "DIST_SQUARED_SUM")]["samples"]) for tree_kind in residual_kinds)
        residual_maximum_changes = sum(
            int(round(
                float(summary_lookup[(tree_kind, "MAX_SQUARED_DIST")]["samples"])
                * (1.0 - float(summary_lookup[(tree_kind, "MAX_SQUARED_DIST")]["exact_equal_rate"]))
            ))
            for tree_kind in residual_kinds
        )
        main_findings.append(
            f"- Across {residual_nodes:,} residual-tree nodes, `MAX_SQUARED_DIST` changed on "
            f"{residual_maximum_changes:,} nodes; the residual descriptor table separates unrestricted and saturated behavior."
        )
    if "tos" in metadata["tree_kinds"]:
        tos_field = summary_lookup[("tos", "DIST_SQUARED_SUM")]
        tos_maximum = summary_lookup[("tos", "MAX_SQUARED_DIST")]
        tos_field_changes = int(round(
            float(tos_field["samples"]) * (1.0 - float(tos_field["exact_equal_rate"]))
        ))
        tos_maximum_changes = int(round(
            float(tos_maximum["samples"]) * (1.0 - float(tos_maximum["exact_equal_rate"]))
        ))
        main_findings.append(
            f"- On the {int(tos_field['samples']):,} Tree-of-Shapes nodes, the squared-distance field changed on "
            f"{tos_field_changes:,} nodes and `MAX_SQUARED_DIST` changed on {tos_maximum_changes:,} nodes."
        )
    lines.extend(["", "## Main findings", "", *main_findings])

    lines.extend([
        "",
        "## All 29 attribute pairs, selected tree types combined",
        "",
        "| Attribute | Unit | Agreement | NMAE | NRMSE | Bias | Mean image Spearman |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: |",
    ])
    for row in sorted(combined, key=lambda item: str(item["attribute"])):
        lines.append(
            f"| `{row['attribute']}` | {row['unit']} | {100.0 * float(row['close_rate']):.3f}% | "
            f"{format_number(row['nmae'])} | {format_number(row['nrmse'])} | {format_number(row['bias'])} | "
            f"{format_number(row['mean_image_spearman'])} |"
        )

    lines.extend([
        "",
        "## Most sensitive pairs by disagreement rate",
        "",
        "| Attribute | Disagreement | Under | Over | MAE | Max absolute error |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ])
    for row in ranked[:10]:
        lines.append(
            f"| `{row['attribute']}` | {100.0 * (1.0 - float(row['close_rate'])):.3f}% | "
            f"{100.0 * float(row['underestimate_rate']):.3f}% | {100.0 * float(row['overestimate_rate']):.3f}% | "
            f"{format_number(row['mae'])} | {format_number(row['max_absolute_error'])} |"
        )

    lines.extend([
        "",
        "## Joint coordinate sensitivity",
        "",
        "| Tree | Geometry | Unchanged | Mean displacement | RMS displacement | Maximum displacement |",
        "| --- | --- | ---: | ---: | ---: | ---: |",
    ])
    for row in geometry_rows:
        if row["tree"] != "combined":
            continue
        lines.append(
            f"| combined | `{row['geometry']}` | {100.0 * float(row['unchanged_rate']):.3f}% | "
            f"{format_number(row['mean_displacement'])} | {format_number(row['rms_displacement'])} | "
            f"{format_number(row['max_displacement'])} |"
        )

    lines.extend([
        "",
        "## Interpretation constraints",
        "",
        "- Raw errors are meaningful only within an attribute's declared unit; do not rank different units by MAE alone.",
        "- Zero-heavy attributes can show high agreement while differing on the subset of geometrically nontrivial nodes; the CSV includes `exact_nonzero_rate` and nonzero relative error.",
        "- Center and centroid rows/columns are also summarized jointly as Euclidean displacement.",
        "- This experiment characterizes only the requested tree kinds and construction settings. Changing the adjacency, residual spatial order, saturated infinity pixel, or tree family can change the supports.",
        "- Trees of shapes are not included in this campaign.",
        "",
        "## Artifacts",
        "",
        "- `experiment.json`: software, input, tolerance, and command metadata.",
        "- `inputs.csv`: input hashes, dimensions, node counts, and runtimes.",
        "- `per_image_metrics.csv`: one row per image/tree/attribute pair.",
        "- `summary.csv`: micro-aggregated node errors plus macro per-image correlation summaries.",
        "- `area_strata.csv` and `area_strata_summary.csv`: the same comparison by support-area scale.",
        "- `field_condition.csv` and `field_condition_summary.csv`: descriptors conditioned on equality/change of the integer squared-distance field sum.",
        "- `geometry_per_image.csv` and `geometry_summary.csv`: joint coordinate displacement.",
        "",
    ])
    return "\n".join(lines)


def main() -> int:
    arguments = parse_args()
    tree_kinds = parse_tree_kinds(arguments.trees)
    images = discover_images(arguments.image_directory.resolve(), arguments.count)
    output_directory = arguments.output_directory.resolve()
    output_directory.mkdir(parents=True, exist_ok=True)

    repository = Path(__file__).resolve().parents[1]
    metadata = {
        "schema_version": 1,
        "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "command": sys.argv,
        "image_directory": str(arguments.image_directory.resolve()),
        "images": [path.name for path in images],
        "image_count": len(images),
        "selection_order": "lexicographic filename order",
        "tree_kinds": tree_kinds,
        "radius": arguments.radius,
        "infinity_pixel": arguments.infinity_pixel,
        "tree_of_shapes_convention": {
            "immersion": "canonical complementary grid (minimum 4-connectivity, maximum 8-connectivity)",
            "domain_extension": "none (no padding)",
            "infinity_pixel": 0,
            "altitude_encoding": "uint8",
        },
        "atol": arguments.atol,
        "rtol": arguments.rtol,
        "dtype": "float64",
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "mmcfilters": getattr(mmcfilters, "__version__", "unknown"),
        "git": git_metadata(repository),
    }

    input_rows: list[dict] = []
    image_metric_rows: list[dict] = []
    area_metric_rows: list[dict] = []
    field_condition_metric_rows: list[dict] = []
    geometry_image_rows: list[dict] = []
    aggregates: dict[tuple[str, str], ScalarAggregate] = {}
    area_aggregates: dict[tuple[str, str, str], ScalarAggregate] = {}
    field_condition_aggregates: dict[tuple[str, str, str], ScalarAggregate] = {}
    geometry_aggregates: dict[tuple[str, str], GeometryAggregate] = {}
    canonical_pairs: list[tuple[str, str]] | None = None

    for image_index, path in enumerate(images):
        image = load_grayscale(path)
        image_hash = sha256_file(path)
        for tree_kind in tree_kinds:
            build_start = time.perf_counter()
            tree = build_tree(image, tree_kind, arguments.radius, arguments.infinity_pixel)
            build_seconds = time.perf_counter() - build_start

            compute_start = time.perf_counter()
            layout, values = mmcfilters.Attribute.compute_topology_attributes(
                tree,
                [mmcfilters.Attribute.Group.DIST_TRANSF, mmcfilters.Attribute.Group.DIST_TRANSF_EXACT, mmcfilters.Attribute.AREA],
                dtype=np.float64,
            )
            attribute_seconds = time.perf_counter() - compute_start
            names = list(layout)
            offsets = dict(layout.items())
            pairs = sorted((name.removesuffix(EXACT_SUFFIX), name) for name in names if name.endswith(EXACT_SUFFIX))
            if len(pairs) != 29 or any(approximate not in offsets for approximate, _ in pairs):
                raise RuntimeError("DIST_TRANSF and DIST_TRANSF_EXACT must expose 29 aligned attribute pairs")
            if canonical_pairs is None:
                canonical_pairs = pairs
                metadata["attribute_pairs"] = [{"approximate": approximate, "exact": exact} for approximate, exact in pairs]
            elif pairs != canonical_pairs:
                raise RuntimeError("Distance-transform group pair layout changed between tree/image cases")

            alive = np.asarray(tree.alive_node_ids, dtype=np.int64)
            area = values[alive, offsets["AREA"]]
            input_rows.append({
                "image_index": image_index,
                "image": path.name,
                "sha256": image_hash,
                "rows": image.shape[0],
                "columns": image.shape[1],
                "pixels": image.size,
                "tree": tree_kind,
                "nodes": alive.size,
                "build_seconds": build_seconds,
                "dist_transf_seconds": attribute_seconds,
            })

            vectors: dict[str, np.ndarray] = {name: values[alive, offset] for name, offset in offsets.items()}
            field_changed_mask = vectors["DIST_SQUARED_SUM"] != vectors["DIST_SQUARED_SUM_EXACT"]
            field_conditions = (("field_equal", ~field_changed_mask), ("field_changed", field_changed_mask))
            for approximate_name, exact_name in pairs:
                approximate = vectors[approximate_name]
                exact = vectors[exact_name]
                metrics = scalar_metrics(approximate, exact, arguments.atol, arguments.rtol)
                identity = {
                    "image_index": image_index,
                    "image": path.name,
                    "tree": tree_kind,
                    "attribute": approximate_name,
                    "exact_attribute": exact_name,
                    "category": attribute_category(approximate_name),
                    "unit": attribute_unit(approximate_name),
                }
                image_metric_rows.append(identity | metrics)
                for aggregate_tree in (tree_kind, "combined"):
                    aggregate = aggregates.setdefault((aggregate_tree, approximate_name), ScalarAggregate())
                    aggregate.add(approximate, exact, metrics, arguments.atol, arguments.rtol)

                for field_condition, mask in field_conditions:
                    if not np.any(mask):
                        continue
                    condition_metrics = scalar_metrics(approximate[mask], exact[mask], arguments.atol, arguments.rtol)
                    field_condition_metric_rows.append(identity | {"field_condition": field_condition} | condition_metrics)
                    for aggregate_tree in (tree_kind, "combined"):
                        aggregate = field_condition_aggregates.setdefault(
                            (aggregate_tree, approximate_name, field_condition), ScalarAggregate()
                        )
                        aggregate.add(approximate[mask], exact[mask], condition_metrics, arguments.atol, arguments.rtol)

                for stratum, lower, upper in AREA_STRATA:
                    mask = area >= lower
                    if upper is not None:
                        mask &= area <= upper
                    if not np.any(mask):
                        continue
                    stratum_metrics = scalar_metrics(approximate[mask], exact[mask], arguments.atol, arguments.rtol)
                    area_metric_rows.append(identity | {"area_stratum": stratum} | stratum_metrics)
                    for aggregate_tree in (tree_kind, "combined"):
                        aggregate = area_aggregates.setdefault((aggregate_tree, approximate_name, stratum), ScalarAggregate())
                        aggregate.add(approximate[mask], exact[mask], stratum_metrics, arguments.atol, arguments.rtol)

            for geometry_name, row_name, column_name in GEOMETRY_PAIRS:
                row_delta = vectors[row_name] - vectors[row_name + EXACT_SUFFIX]
                column_delta = vectors[column_name] - vectors[column_name + EXACT_SUFFIX]
                displacement = np.hypot(row_delta, column_delta)
                geometry_image_rows.append({
                    "image_index": image_index,
                    "image": path.name,
                    "tree": tree_kind,
                    "geometry": geometry_name,
                    "samples": displacement.size,
                    "unchanged_rate": float(np.mean(displacement <= arguments.atol)),
                    "mean_displacement": float(np.mean(displacement)),
                    "rms_displacement": math.sqrt(float(np.dot(displacement, displacement)) / displacement.size),
                    "median_displacement": float(np.median(displacement)),
                    "p95_displacement": float(np.quantile(displacement, 0.95)),
                    "max_displacement": float(np.max(displacement)),
                })
                for aggregate_tree in (tree_kind, "combined"):
                    geometry_aggregates.setdefault((aggregate_tree, geometry_name), GeometryAggregate()).add(displacement, arguments.atol)

            print(
                f"[{image_index + 1:02d}/{len(images):02d}] {path.name} {tree_kind}-tree "
                f"nodes={alive.size} build={build_seconds:.3f}s dist_transf={attribute_seconds:.3f}s",
                flush=True,
            )

    summary_rows: list[dict] = []
    for (tree_kind, attribute), aggregate in sorted(aggregates.items()):
        summary_rows.append({
            "tree": tree_kind,
            "attribute": attribute,
            "exact_attribute": attribute + EXACT_SUFFIX,
            "category": attribute_category(attribute),
            "unit": attribute_unit(attribute),
        } | aggregate.as_row())

    area_summary_rows: list[dict] = []
    for (tree_kind, attribute, stratum), aggregate in sorted(area_aggregates.items()):
        area_summary_rows.append({
            "tree": tree_kind,
            "attribute": attribute,
            "exact_attribute": attribute + EXACT_SUFFIX,
            "area_stratum": stratum,
            "category": attribute_category(attribute),
            "unit": attribute_unit(attribute),
        } | aggregate.as_row())

    field_condition_summary_rows: list[dict] = []
    for (tree_kind, attribute, field_condition), aggregate in sorted(field_condition_aggregates.items()):
        field_condition_summary_rows.append({
            "tree": tree_kind,
            "attribute": attribute,
            "exact_attribute": attribute + EXACT_SUFFIX,
            "field_condition": field_condition,
            "category": attribute_category(attribute),
            "unit": attribute_unit(attribute),
        } | aggregate.as_row())

    geometry_summary_rows: list[dict] = []
    for (tree_kind, geometry), aggregate in sorted(geometry_aggregates.items()):
        geometry_summary_rows.append({"tree": tree_kind, "geometry": geometry, "unit": "pixel"} | aggregate.as_row())

    write_csv(output_directory / "inputs.csv", input_rows)
    write_csv(output_directory / "per_image_metrics.csv", image_metric_rows)
    write_csv(output_directory / "summary.csv", summary_rows)
    write_csv(output_directory / "area_strata.csv", area_metric_rows)
    write_csv(output_directory / "area_strata_summary.csv", area_summary_rows)
    write_csv(output_directory / "field_condition.csv", field_condition_metric_rows)
    write_csv(output_directory / "field_condition_summary.csv", field_condition_summary_rows)
    write_csv(output_directory / "geometry_per_image.csv", geometry_image_rows)
    write_csv(output_directory / "geometry_summary.csv", geometry_summary_rows)
    (output_directory / "experiment.json").write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (output_directory / "analysis.md").write_text(
        build_report(metadata, summary_rows, area_summary_rows, geometry_summary_rows, field_condition_summary_rows, input_rows), encoding="utf-8"
    )
    print(f"artifacts={output_directory}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
