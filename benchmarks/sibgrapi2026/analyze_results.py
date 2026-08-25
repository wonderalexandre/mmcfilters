#!/usr/bin/env python3
"""Validate raw SIBGRAPI 2026 measurements and generate the paper table."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path


RESOLUTIONS = ("480p", "720p", "1080p")
DIMENSIONS = {"480p": (853, 480), "720p": (1280, 720), "1080p": (1920, 1080)}
HIERARCHIES = ("max_tree", "min_tree", "tree_of_shapes")
METHODS = {
    "max_tree": ("proposed", "ref6"),
    "min_tree": ("proposed", "ref6"),
    "tree_of_shapes": ("proposed", "ref5_original"),
}
BASELINES = {"max_tree": "ref6", "min_tree": "ref6", "tree_of_shapes": "ref5_original"}
CONNECTIVITY = {"max_tree": "8", "min_tree": "8", "tree_of_shapes": "4/8"}
INTEGER_FIELDS = (
    "image_index",
    "width",
    "height",
    "pixels",
    "run",
    "order_position",
    "nodes",
    "baseline_nodes",
    "compared_nodes",
    "mismatch_nodes",
    "q1_mismatch_nodes",
    "q2_mismatch_nodes",
    "qd_mismatch_nodes",
    "q3_mismatch_nodes",
    "q4_mismatch_nodes",
    "max_abs_family_error",
    "checksum",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--count", type=int, default=100)
    args = parser.parse_args()
    if args.start < 0 or args.count < 1 or args.start + args.count > 100:
        parser.error("require START >= 0, COUNT >= 1, and START + COUNT <= 100")
    return args


def read_rows(filename: Path, resolution: str) -> list[dict[str, object]]:
    with filename.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        rows: list[dict[str, object]] = []
        for source in reader:
            row: dict[str, object] = dict(source)
            for field in INTEGER_FIELDS:
                row[field] = int(source[field])
            for field in ("time_ms", "tree_build_ms", "bridge_import_ms"):
                row[field] = float(source[field])
            rows.append(row)
    if any(row["resolution"] != resolution for row in rows):
        raise ValueError(f"resolution column disagrees with {filename.name}")
    return rows


def validate_rows(rows: list[dict[str, object]], resolution: str, start: int, count: int) -> None:
    expected_rows = count * len(HIERARCHIES) * 2 * 3
    if len(rows) != expected_rows:
        raise ValueError(f"{resolution}: expected {expected_rows} rows, found {len(rows)}")
    expected_indices = set(range(start, start + count))
    if {int(row["image_index"]) for row in rows} != expected_indices:
        raise ValueError(f"{resolution}: image-index set does not match the requested range")
    expected_width, expected_height = DIMENSIONS[resolution]
    if any((row["width"], row["height"]) != (expected_width, expected_height) for row in rows):
        raise ValueError(f"{resolution}: unexpected image dimensions")
    if {str(row["contract_mode"]) for row in rows} != {"UNCHECKED"}:
        raise ValueError(f"{resolution}: paper timings require a Release/UNCHECKED build")

    groups: dict[tuple[int, str, str], list[dict[str, object]]] = defaultdict(list)
    pairs: dict[tuple[int, str, int], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        image = int(row["image_index"])
        hierarchy = str(row["hierarchy"])
        method = str(row["method"])
        if hierarchy not in METHODS or method not in METHODS[hierarchy]:
            raise ValueError(f"{resolution}: unexpected hierarchy/method pair {hierarchy}/{method}")
        if row["image"] != f"val_{image:03d}.png":
            raise ValueError(f"{resolution}: unexpected filename for image index {image}")
        if row["connectivity"] != CONNECTIVITY[hierarchy]:
            raise ValueError(f"{resolution}: unexpected connectivity for {hierarchy}")
        if row["reference_method"] != "proposed":
            raise ValueError(f"{resolution}: unexpected comparison reference for {hierarchy}/{method}")
        if not math.isfinite(float(row["time_ms"])) or float(row["time_ms"]) <= 0.0:
            raise ValueError(f"{resolution}: invalid elapsed time for {hierarchy}/{method}")
        groups[(image, hierarchy, method)].append(row)
        pairs[(image, hierarchy, int(row["run"]))].append(row)

    expected_groups = count * len(HIERARCHIES) * 2
    if len(groups) != expected_groups:
        raise ValueError(f"{resolution}: expected {expected_groups} method groups, found {len(groups)}")
    for key, group in groups.items():
        if len(group) != 3 or {int(row["run"]) for row in group} != {0, 1, 2}:
            raise ValueError(f"{resolution}: {key} does not contain exactly runs 0, 1, and 2")
        for field in ("checksum", "nodes", "baseline_nodes", "compared_nodes", "mismatch_nodes", "q1_mismatch_nodes", "q2_mismatch_nodes",
                      "qd_mismatch_nodes", "q3_mismatch_nodes", "q4_mismatch_nodes", "max_abs_family_error"):
            if len({row[field] for row in group}) != 1:
                raise ValueError(f"{resolution}: {field} changes across repetitions for {key}")

    for key, pair in pairs.items():
        hierarchy = key[1]
        if len(pair) != 2 or {str(row["method"]) for row in pair} != set(METHODS[hierarchy]):
            raise ValueError(f"{resolution}: incomplete paired run {key}")
        if {int(row["order_position"]) for row in pair} != {0, 1}:
            raise ValueError(f"{resolution}: invalid alternating-order positions for {key}")
        image, hierarchy, run = key
        order_seed = image * len(HIERARCHIES) + HIERARCHIES.index(hierarchy)
        expected_proposed_position = 1 if (order_seed + run) % 2 == 0 else 0
        proposed_position = next(int(row["order_position"]) for row in pair if row["method"] == "proposed")
        if proposed_position != expected_proposed_position:
            raise ValueError(f"{resolution}: paired order does not alternate as specified for {key}")

    for image in expected_indices:
        for hierarchy in HIERARCHIES:
            proposed = groups[(image, hierarchy, "proposed")]
            baseline_name = BASELINES[hierarchy]
            baseline = groups[(image, hierarchy, baseline_name)]
            if hierarchy != "tree_of_shapes" and proposed[0]["checksum"] != baseline[0]["checksum"]:
                raise ValueError(f"{resolution}: component-tree checksums disagree for image {image}, {hierarchy}")
            if any(int(row["mismatch_nodes"]) != 0 for row in proposed):
                raise ValueError(f"{resolution}: proposed rows unexpectedly report mismatches")
            if hierarchy != "tree_of_shapes" and any(int(row["mismatch_nodes"]) != 0 for row in baseline):
                raise ValueError(f"{resolution}: specialized component-tree reference disagrees with the proposed method")


def image_means(rows: list[dict[str, object]]) -> dict[tuple[str, str, int], float]:
    samples: dict[tuple[str, str, int], list[float]] = defaultdict(list)
    for row in rows:
        samples[(str(row["hierarchy"]), str(row["method"]), int(row["image_index"]))].append(float(row["time_ms"]))
    return {key: statistics.fmean(values) for key, values in samples.items()}


def write_summary(output: Path, all_rows: dict[str, list[dict[str, object]]]) -> dict[tuple[str, str, str], float]:
    table_means: dict[tuple[str, str, str], float] = {}
    with (output / "summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ("resolution", "hierarchy", "method", "images", "repetitions_per_image", "mean_ms", "between_image_sd_ms", "min_image_mean_ms", "max_image_mean_ms")
        )
        for resolution in RESOLUTIONS:
            means = image_means(all_rows[resolution])
            for hierarchy in HIERARCHIES:
                for method in METHODS[hierarchy]:
                    values = [value for (h, m, _), value in means.items() if (h, m) == (hierarchy, method)]
                    overall = statistics.fmean(values)
                    table_means[(resolution, hierarchy, method)] = overall
                    writer.writerow(
                        (
                            resolution,
                            hierarchy,
                            method,
                            len(values),
                            3,
                            f"{overall:.6f}",
                            f"{statistics.stdev(values):.6f}" if len(values) > 1 else "0.000000",
                            f"{min(values):.6f}",
                            f"{max(values):.6f}",
                        )
                    )
    return table_means


def write_speedups(output: Path, all_rows: dict[str, list[dict[str, object]]], table_means: dict[tuple[str, str, str], float]) -> None:
    with (output / "speedups.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ("resolution", "hierarchy", "baseline", "ratio_of_means_baseline_over_proposed", "mean_per_image_ratio", "sd_per_image_ratio", "proposed_time_difference_percent")
        )
        for resolution in RESOLUTIONS:
            means = image_means(all_rows[resolution])
            for hierarchy in HIERARCHIES:
                baseline = BASELINES[hierarchy]
                indices = sorted(index for (h, m, index) in means if (h, m) == (hierarchy, "proposed"))
                ratios = [means[(hierarchy, baseline, index)] / means[(hierarchy, "proposed", index)] for index in indices]
                proposed = table_means[(resolution, hierarchy, "proposed")]
                reference = table_means[(resolution, hierarchy, baseline)]
                writer.writerow(
                    (
                        resolution,
                        hierarchy,
                        baseline,
                        f"{reference / proposed:.6f}",
                        f"{statistics.fmean(ratios):.6f}",
                        f"{statistics.stdev(ratios):.6f}" if len(ratios) > 1 else "0.000000",
                        f"{100.0 * (proposed - reference) / reference:.6f}",
                    )
                )


def write_correctness(output: Path, all_rows: dict[str, list[dict[str, object]]]) -> None:
    fields = ("mismatch_nodes", "q1_mismatch_nodes", "q2_mismatch_nodes", "qd_mismatch_nodes", "q3_mismatch_nodes", "q4_mismatch_nodes")
    with (output / "correctness.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("resolution", "hierarchy", "method", "images", *fields, "max_abs_family_error", "interpretation"))
        for resolution in RESOLUTIONS:
            unique: dict[tuple[int, str, str], dict[str, object]] = {}
            for row in all_rows[resolution]:
                key = (int(row["image_index"]), str(row["hierarchy"]), str(row["method"]))
                unique.setdefault(key, row)
            for hierarchy in HIERARCHIES:
                for method in METHODS[hierarchy]:
                    selected = [row for (_, h, m), row in unique.items() if (h, m) == (hierarchy, method)]
                    interpretation = (
                        "timing_only_original_patterns_differ_including_Q2_QD"
                        if hierarchy == "tree_of_shapes" and method == "ref5_original"
                        else "exact_agreement"
                    )
                    writer.writerow(
                        (
                            resolution,
                            hierarchy,
                            method,
                            len(selected),
                            *(sum(int(row[field]) for row in selected) for field in fields),
                            max(int(row["max_abs_family_error"]) for row in selected),
                            interpretation,
                        )
                    )


def write_latex_table(output: Path, means: dict[tuple[str, str, str], float], image_count: int) -> None:
    names = {
        "max_tree": ("Max-tree", "Silva et al.~\\cite{silva2020incremental}"),
        "min_tree": ("Min-tree", "Silva et al.~\\cite{silva2020incremental}"),
        "tree_of_shapes": ("ToS", "da Silva et al.~\\cite{dasilva2019incremental} (original)"),
    }
    lines = [
        "\\begin{table}[t]",
        f"\\caption{{Mean warm attribute-API time (ms) over {image_count} images. Tree construction, topology import, and decision-table loading are excluded.}}",
        "\\label{tab:bitquad-runtime}",
        "\\centering",
        "\\setlength{\\tabcolsep}{3.2pt}",
        "\\renewcommand{\\arraystretch}{1.05}",
        "\\resizebox{\\columnwidth}{!}{%",
        "\\begin{tabular}{@{}llrrr@{}}",
        "\\toprule",
        "Tree & Method & $853{\\times}480$ & $1280{\\times}720$ & $1920{\\times}1080$ \\\\",
        "\\midrule",
    ]
    for hierarchy_index, hierarchy in enumerate(HIERARCHIES):
        tree_name, baseline_name = names[hierarchy]
        baseline = BASELINES[hierarchy]
        proposed_values = [means[(resolution, hierarchy, "proposed")] for resolution in RESOLUTIONS]
        baseline_values = [means[(resolution, hierarchy, baseline)] for resolution in RESOLUTIONS]
        proposed_cells = [f"\\textbf{{{value:.2f}}}" if value <= other else f"{value:.2f}" for value, other in zip(proposed_values, baseline_values)]
        baseline_cells = [f"\\textbf{{{value:.2f}}}" if value < other else f"{value:.2f}" for value, other in zip(baseline_values, proposed_values)]
        lines.append(f"{tree_name} & Ours & " + " & ".join(proposed_cells) + " \\\\")
        lines.append(" & " + baseline_name + " & " + " & ".join(baseline_cells) + " \\\\")
        if hierarchy_index + 1 != len(HIERARCHIES):
            lines.append("\\addlinespace")
    lines.extend(("\\bottomrule", "\\end{tabular}", "}", "\\end{table}", ""))
    (output / "table.tex").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    all_rows: dict[str, list[dict[str, object]]] = {}
    for resolution in RESOLUTIONS:
        filename = args.output_dir / f"raw-{resolution}.csv"
        if not filename.is_file():
            raise ValueError(f"missing raw result: {filename}")
        rows = read_rows(filename, resolution)
        validate_rows(rows, resolution, args.start, args.count)
        all_rows[resolution] = rows
    table_means = write_summary(args.output_dir, all_rows)
    write_speedups(args.output_dir, all_rows, table_means)
    write_correctness(args.output_dir, all_rows)
    write_latex_table(args.output_dir, table_means, args.count)
    print(f"validated three repetitions per image; generated {args.output_dir / 'table.tex'}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, KeyError, ZeroDivisionError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
