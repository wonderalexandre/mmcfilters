#!/usr/bin/env python3
"""Create publication-oriented tables from a persisted API benchmark experiment."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import statistics
import sys
from collections import defaultdict
from pathlib import Path


OUTCOME_FIELDS = (
    "steps",
    "primary_nodes_removed",
    "secondary_nodes_removed",
    "complete_validation_commits",
    "incremental_validation_commits",
    "light_threshold",
    "medium_threshold",
    "heavy_threshold",
)
BOOTSTRAP_REPLICATES = 20_000


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyze a persisted CHECKED/UNCHECKED scientific API experiment.")
    parser.add_argument("experiment_dir", type=Path, help="directory created by compare_validation_modes.py --output-dir")
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_json_lines(path: Path) -> list[dict[str, object]]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"]
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(lines)


def percentile(sorted_values: list[float], probability: float) -> float:
    if not sorted_values:
        raise RuntimeError("cannot compute a percentile of an empty sample")
    position = probability * (len(sorted_values) - 1)
    lower = int(position)
    upper = min(lower + 1, len(sorted_values) - 1)
    fraction = position - lower
    return sorted_values[lower] * (1.0 - fraction) + sorted_values[upper] * fraction


def cliffs_delta(checked: list[float], unchecked: list[float]) -> float:
    greater = sum(unchecked_value > checked_value for unchecked_value in unchecked for checked_value in checked)
    lower = sum(unchecked_value < checked_value for unchecked_value in unchecked for checked_value in checked)
    return (greater - lower) / (len(checked) * len(unchecked))


def bootstrap_relative_interval(scenario: str, checked: list[float], unchecked: list[float]) -> tuple[float, float]:
    seed = int.from_bytes(hashlib.sha256(scenario.encode("utf-8")).digest()[:8], "big")
    generator = random.Random(seed)
    estimates: list[float] = []
    for _ in range(BOOTSTRAP_REPLICATES):
        checked_center = statistics.median(generator.choices(checked, k=len(checked)))
        unchecked_center = statistics.median(generator.choices(unchecked, k=len(unchecked)))
        estimates.append(100.0 * (unchecked_center / checked_center - 1.0))
    estimates.sort()
    return percentile(estimates, 0.025), percentile(estimates, 0.975)


def compute_sample_statistics(
    metadata: dict[str, object], summary: list[dict[str, str]], samples_path: Path
) -> list[dict[str, object]]:
    if not samples_path.is_file():
        return []
    sample_rows = read_csv(samples_path)
    benchmark = metadata.get("benchmark_metadata", {})
    expected_repetitions = int(benchmark["repetitions"])
    expected_processes = int(metadata["process_runs_per_mode"])
    summary_by_scenario = {row["scenario"]: row for row in summary}
    grouped: dict[tuple[str, str, int], list[tuple[int, float]]] = defaultdict(list)
    for row in sample_rows:
        try:
            scenario = row["scenario"]
            mode = row["contract_mode"]
            process_run = int(row["process_run"])
            sample_index = int(row["sample_index"])
            sample_ms = float(row["sample_ms"])
        except (KeyError, ValueError) as error:
            raise RuntimeError("timed-sample table contains a malformed row") from error
        if row.get("record") != "sample" or scenario not in summary_by_scenario:
            raise RuntimeError("timed-sample table contains an unexpected record or scenario")
        if mode not in ("CHECKED", "UNCHECKED") or not 1 <= process_run <= expected_processes:
            raise RuntimeError(f"{scenario} has an invalid contract mode or process-run index")
        if not math.isfinite(sample_ms) or sample_ms < 0.0:
            raise RuntimeError("timed-sample table contains an invalid duration")
        if row.get("checksum") != summary_by_scenario[scenario]["checksum"]:
            raise RuntimeError(f"{scenario} timed samples have an inconsistent checksum")
        grouped[(scenario, mode, process_run)].append((sample_index, sample_ms))

    statistics_rows: list[dict[str, object]] = []
    for scenario, summary_row in sorted(summary_by_scenario.items()):
        process_samples: dict[str, list[list[float]]] = {"CHECKED": [], "UNCHECKED": []}
        for mode in ("CHECKED", "UNCHECKED"):
            for process_run in range(1, expected_processes + 1):
                indexed = sorted(grouped.get((scenario, mode, process_run), []))
                if [index for index, _ in indexed] != list(range(1, expected_repetitions + 1)):
                    raise RuntimeError(f"{scenario} has an incomplete timed-sample sequence for {mode} process {process_run}")
                process_samples[mode].append([value for _, value in indexed])

        checked_medians = [statistics.median(values) for values in process_samples["CHECKED"]]
        unchecked_medians = [statistics.median(values) for values in process_samples["UNCHECKED"]]
        if any(value <= 0.0 for value in (*checked_medians, *unchecked_medians)):
            raise RuntimeError(f"{scenario} has a non-positive process-median duration")
        checked_center = statistics.median(checked_medians)
        unchecked_center = statistics.median(unchecked_medians)
        if not math.isclose(checked_center, float(summary_row["checked_median_ms"]), rel_tol=1e-6, abs_tol=2e-6):
            raise RuntimeError(f"{scenario} CHECKED summary is inconsistent with its raw samples")
        if not math.isclose(unchecked_center, float(summary_row["unchecked_median_ms"]), rel_tol=1e-6, abs_tol=2e-6):
            raise RuntimeError(f"{scenario} UNCHECKED summary is inconsistent with its raw samples")

        relative = 100.0 * (unchecked_center / checked_center - 1.0)
        interval_low, interval_high = bootstrap_relative_interval(scenario, checked_medians, unchecked_medians)
        if interval_high < 0.0:
            inference = "UNCHECKED_FASTER"
        elif interval_low > 0.0:
            inference = "UNCHECKED_SLOWER"
        else:
            inference = "INCONCLUSIVE"
        checked_flat = [value for process in process_samples["CHECKED"] for value in process]
        unchecked_flat = [value for process in process_samples["UNCHECKED"] for value in process]
        statistics_rows.append(
            {
                "scenario": scenario,
                "suite": summary_row["suite"],
                "scope": summary_row["scope"],
                "name": summary_row["name"],
                "process_runs_per_mode": expected_processes,
                "repetitions_per_process": expected_repetitions,
                "samples_per_mode": expected_processes * expected_repetitions,
                "checked_median_of_process_medians_ms": checked_center,
                "unchecked_median_of_process_medians_ms": unchecked_center,
                "unchecked_vs_checked_percent": relative,
                "bootstrap_ci95_low_percent": interval_low,
                "bootstrap_ci95_high_percent": interval_high,
                "bootstrap_replicates": BOOTSTRAP_REPLICATES,
                "cliffs_delta_process_medians": cliffs_delta(checked_medians, unchecked_medians),
                "inference": inference,
                "checked_process_medians_ms": json.dumps(checked_medians, separators=(",", ":")),
                "unchecked_process_medians_ms": json.dumps(unchecked_medians, separators=(",", ":")),
                "checked_pooled_sample_median_ms": statistics.median(checked_flat),
                "unchecked_pooled_sample_median_ms": statistics.median(unchecked_flat),
                "checksum": summary_row["checksum"],
            }
        )
    return statistics_rows


def main() -> int:
    arguments = parse_arguments()
    experiment_dir = arguments.experiment_dir.resolve()
    metadata_path = experiment_dir / "experiment.json"
    summary_path = experiment_dir / "summary.csv"
    raw_path = experiment_dir / "raw.jsonl"
    for path in (metadata_path, summary_path, raw_path):
        if not path.is_file():
            raise RuntimeError(f"missing experiment artifact: {path}")

    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    summary = read_csv(summary_path)
    raw = load_json_lines(raw_path)
    if not summary:
        raise RuntimeError("experiment summary is empty")

    samples_path = experiment_dir / "samples.csv"
    if metadata.get("samples_captured") and not samples_path.is_file():
        raise RuntimeError(f"experiment metadata declares timed samples, but the artifact is missing: {samples_path}")
    sample_statistics = compute_sample_statistics(metadata, summary, samples_path)
    if sample_statistics:
        write_csv(experiment_dir / "sample-statistics.csv", sample_statistics, list(sample_statistics[0]))

    summary_by_scenario = {row["scenario"]: row for row in summary}
    relative_differences = [float(row["unchecked_vs_checked_percent"]) for row in summary]
    checked_total = sum(float(row["checked_median_ms"]) for row in summary)
    unchecked_total = sum(float(row["unchecked_median_ms"]) for row in summary)

    scope_groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in summary:
        scope_groups[(row["suite"], row["scope"])].append(row)
    scope_rows: list[dict[str, object]] = []
    for (suite, scope), rows in sorted(scope_groups.items()):
        checked = [float(row["checked_median_ms"]) for row in rows]
        unchecked = [float(row["unchecked_median_ms"]) for row in rows]
        scope_rows.append(
            {
                "suite": suite,
                "scope": scope,
                "scenario_count": len(rows),
                "checked_sum_of_medians_ms": sum(checked),
                "unchecked_sum_of_medians_ms": sum(unchecked),
                "checked_median_scenario_ms": statistics.median(checked),
                "unchecked_median_scenario_ms": statistics.median(unchecked),
            }
        )
    write_csv(
        experiment_dir / "scope-summary.csv",
        scope_rows,
        [
            "suite",
            "scope",
            "scenario_count",
            "checked_sum_of_medians_ms",
            "unchecked_sum_of_medians_ms",
            "checked_median_scenario_ms",
            "unchecked_median_scenario_ms",
        ],
    )

    scope_pair_definitions = (
        (
            "casf_area_sequence",
            "casf.end_to_end.pipeline_area_sequence",
            "casf.established_input.incremental_area_sequence",
        ),
    )
    scope_pair_rows: list[dict[str, object]] = []
    for pair_name, end_to_end_name, established_name in scope_pair_definitions:
        end_to_end = summary_by_scenario.get(end_to_end_name)
        established = summary_by_scenario.get(established_name)
        if end_to_end is None or established is None:
            continue
        if end_to_end["checksum"] != established["checksum"]:
            raise RuntimeError(f"scope pair checksum mismatch: {pair_name}")
        row: dict[str, object] = {
            "pair": pair_name,
            "end_to_end_scenario": end_to_end_name,
            "established_input_scenario": established_name,
            "checksum": end_to_end["checksum"],
        }
        for mode in ("checked", "unchecked"):
            end_to_end_ms = float(end_to_end[f"{mode}_median_ms"])
            established_ms = float(established[f"{mode}_median_ms"])
            row[f"{mode}_end_to_end_ms"] = end_to_end_ms
            row[f"{mode}_established_input_ms"] = established_ms
            row[f"{mode}_setup_overhead_ms"] = end_to_end_ms - established_ms
            row[f"{mode}_setup_share_percent"] = 100.0 * (end_to_end_ms - established_ms) / end_to_end_ms
        scope_pair_rows.append(row)
    scope_pair_fields = [
        "pair",
        "end_to_end_scenario",
        "established_input_scenario",
        "checked_end_to_end_ms",
        "checked_established_input_ms",
        "checked_setup_overhead_ms",
        "checked_setup_share_percent",
        "unchecked_end_to_end_ms",
        "unchecked_established_input_ms",
        "unchecked_setup_overhead_ms",
        "unchecked_setup_share_percent",
        "checksum",
    ]
    write_csv(experiment_dir / "scope-pairs.csv", scope_pair_rows, scope_pair_fields)

    attribute_rows: list[dict[str, object]] = []
    for sequential_name, sequential in sorted(summary_by_scenario.items()):
        suffix = "_sequential_scalars"
        if not sequential_name.endswith(suffix):
            continue
        grouped_name = sequential_name[: -len(suffix)]
        grouped = summary_by_scenario.get(grouped_name)
        if grouped is None:
            continue
        if grouped["checksum"] != sequential["checksum"]:
            raise RuntimeError(f"attribute pair checksum mismatch: {grouped_name}")
        for mode in ("checked", "unchecked"):
            grouped_ms = float(grouped[f"{mode}_median_ms"])
            sequential_ms = float(sequential[f"{mode}_median_ms"])
            attribute_rows.append(
                {
                    "contract_mode": mode.upper(),
                    "suite": grouped["suite"],
                    "name": grouped["name"],
                    "grouped_ms": grouped_ms,
                    "sequential_scalars_ms": sequential_ms,
                    "sequential_over_grouped_ratio": sequential_ms / grouped_ms,
                    "grouped_time_saved_percent": 100.0 * (1.0 - grouped_ms / sequential_ms),
                    "checksum": grouped["checksum"],
                }
            )
    write_csv(
        experiment_dir / "attribute-pairs.csv",
        attribute_rows,
        [
            "contract_mode",
            "suite",
            "name",
            "grouped_ms",
            "sequential_scalars_ms",
            "sequential_over_grouped_ratio",
            "grouped_time_saved_percent",
            "checksum",
        ],
    )

    scenario_records = [record for record in raw if record.get("record") == "scenario"]
    invariant_outcomes: dict[str, dict[str, object]] = {}
    for record in scenario_records:
        name = str(record["name"])
        if not name.startswith("casf."):
            continue
        outcome = {field: record.get(field, 0) for field in OUTCOME_FIELDS}
        previous = invariant_outcomes.setdefault(name, outcome)
        if outcome != previous:
            raise RuntimeError(f"CASF outcome changed across process runs: {name}")
    casf_rows: list[dict[str, object]] = []
    for name, outcome in sorted(invariant_outcomes.items()):
        summary_row = summary_by_scenario[name]
        casf_rows.append(
            {
                "scenario": name,
                "checked_median_ms": float(summary_row["checked_median_ms"]),
                "unchecked_median_ms": float(summary_row["unchecked_median_ms"]),
                "unchecked_vs_checked_percent": float(summary_row["unchecked_vs_checked_percent"]),
                **outcome,
                "checksum": summary_row["checksum"],
            }
        )
    write_csv(
        experiment_dir / "casf.csv",
        casf_rows,
        [
            "scenario",
            "checked_median_ms",
            "unchecked_median_ms",
            "unchecked_vs_checked_percent",
            *OUTCOME_FIELDS,
            "checksum",
        ],
    )

    benchmark = metadata.get("benchmark_metadata", {})
    git = metadata.get("git", {})
    host = metadata.get("host", {})
    faster_count = sum(value < 0 for value in relative_differences)
    report_lines = [
        "# Scientific API validation-contract experiment",
        "",
        "## Experimental context",
        "",
        f"- Workload: `{benchmark.get('workload', 'unknown')}` ({benchmark.get('input_source', 'unknown')})",
        f"- Domain: {benchmark.get('rows', '?')} x {benchmark.get('cols', '?')}",
        f"- Repetitions per process: {benchmark.get('repetitions', '?')}",
        f"- Process runs per mode: {metadata.get('process_runs_per_mode', '?')}",
        f"- Git commit: `{git.get('commit', 'unknown')}`; dirty worktree: `{git.get('dirty', 'unknown')}`",
        f"- CPU: {host.get('cpu', 'unknown')}",
        f"- Operating system: {host.get('operating_system', 'unknown')}",
        f"- Compiler: {benchmark.get('compiler', 'unknown')}",
        "",
        "All CHECKED and UNCHECKED scenario checksums matched; the runner would have aborted on the first mismatch.",
        "",
        "## Validation-policy comparison",
        "",
        f"The experiment contains {len(summary)} scenarios. UNCHECKED was faster in {faster_count} scenarios. "
        f"The median per-scenario time difference (UNCHECKED relative to CHECKED) was {statistics.median(relative_differences):+.2f}%.",
        "",
        f"The descriptive sum of scenario medians was {checked_total:.3f} ms in CHECKED and {unchecked_total:.3f} ms in UNCHECKED "
        f"({100.0 * (unchecked_total / checked_total - 1.0):+.2f}%). These sums are not an end-to-end application runtime because scenarios are independent.",
        "",
    ]
    if sample_statistics:
        supported_faster = [row for row in sample_statistics if row["inference"] == "UNCHECKED_FASTER"]
        supported_slower = [row for row in sample_statistics if row["inference"] == "UNCHECKED_SLOWER"]
        inconclusive = len(sample_statistics) - len(supported_faster) - len(supported_slower)
        strongest_supported = sorted(supported_faster, key=lambda row: float(row["unchecked_vs_checked_percent"]))[:12]
        report_lines.extend(
            [
                "## Raw-sample inference",
                "",
                f"Every scenario retained {sample_statistics[0]['samples_per_mode']} timed samples per contract mode, nested in "
                f"{sample_statistics[0]['process_runs_per_mode']} independent process runs. Confidence intervals use a deterministic "
                f"{BOOTSTRAP_REPLICATES}-replicate cluster bootstrap over process medians; repetitions inside one process are not treated as independent experiments.",
                "",
                f"The 95% interval supports UNCHECKED as faster in {len(supported_faster)} scenarios and slower in {len(supported_slower)}; "
                f"{inconclusive} scenarios include zero and remain inconclusive.",
                "",
            ]
        )
        if strongest_supported:
            report_lines.extend(
                [
                    markdown_table(
                        ["Scenario", "Difference", "Bootstrap 95% CI", "Cliff's delta"],
                        [
                            [
                                str(row["scenario"]),
                                f"{float(row['unchecked_vs_checked_percent']):+.2f}%",
                                f"[{float(row['bootstrap_ci95_low_percent']):+.2f}%, {float(row['bootstrap_ci95_high_percent']):+.2f}%]",
                                f"{float(row['cliffs_delta_process_medians']):+.2f}",
                            ]
                            for row in strongest_supported
                        ],
                    ),
                    "",
                ]
            )
    else:
        report_lines.extend(
            [
                "## Raw-sample inference",
                "",
                "This experiment predates timed-sample capture. Re-run it with `compare_validation_modes.py --capture-samples` before making inferential claims.",
                "",
            ]
        )
    report_lines.extend(
        [
            "## Measurement scopes",
            "",
            markdown_table(
                ["Suite", "Scope", "Scenarios", "CHECKED sum (ms)", "UNCHECKED sum (ms)"],
                [
                    [
                        str(row["suite"]),
                        str(row["scope"]),
                        str(row["scenario_count"]),
                        f"{float(row['checked_sum_of_medians_ms']):.3f}",
                        f"{float(row['unchecked_sum_of_medians_ms']):.3f}",
                    ]
                    for row in scope_rows
                ],
            ),
            "",
            "Scope totals are descriptive only: `end_to_end` and `established_input` contain different scenario sets and must not be treated as paired observations.",
            "",
            "The following scope comparison is paired by scientific output: both CASF scenarios apply the same fixed threshold sequence and have the same checksum; only construction/setup is excluded from `established_input`.",
            "",
            markdown_table(
                ["Pair", "Mode", "End-to-end (ms)", "Established (ms)", "Setup (ms)", "Setup share"],
                [
                    [
                        str(row["pair"]),
                        mode.upper(),
                        f"{float(row[f'{mode}_end_to_end_ms']):.3f}",
                        f"{float(row[f'{mode}_established_input_ms']):.3f}",
                        f"{float(row[f'{mode}_setup_overhead_ms']):.3f}",
                        f"{float(row[f'{mode}_setup_share_percent']):.2f}%",
                    ]
                    for row in scope_pair_rows
                    for mode in ("checked", "unchecked")
                ],
            ),
            "",
            "## Attribute groups and bundles",
            "",
        ]
    )
    if attribute_rows:
        report_lines.append(
            markdown_table(
                ["Mode", "Bundle/group", "Grouped (ms)", "Scalars (ms)", "Scalars/grouped"],
                [
                    [
                        str(row["contract_mode"]),
                        f"{row['suite']}.{row['name']}",
                        f"{float(row['grouped_ms']):.3f}",
                        f"{float(row['sequential_scalars_ms']):.3f}",
                        f"{float(row['sequential_over_grouped_ratio']):.2f}x",
                    ]
                    for row in attribute_rows
                ],
            )
        )
    else:
        report_lines.append("No grouped/sequential attribute pairs were selected by this workload.")
    report_lines.extend(
        [
            "",
            "## CASF incremental work",
            "",
            markdown_table(
                ["Scenario", "CHECKED (ms)", "UNCHECKED (ms)", "Primary removed", "Secondary removed", "Incremental commits"],
                [
                    [
                        str(row["scenario"]),
                        f"{float(row['checked_median_ms']):.3f}",
                        f"{float(row['unchecked_median_ms']):.3f}",
                        str(row["primary_nodes_removed"]),
                        str(row["secondary_nodes_removed"]),
                        str(row["incremental_validation_commits"]),
                    ]
                    for row in casf_rows
                ],
            ),
            "",
            "The full numeric tables are stored in `summary.csv`, `scope-summary.csv`, `scope-pairs.csv`, `attribute-pairs.csv`, and `casf.csv`. When samples are captured, `samples.csv`, `samples.jsonl`, and `sample-statistics.csv` retain repetition-level data and inference; process-level observations remain in `raw.jsonl` and `raw/`.",
            "",
        ]
    )
    (experiment_dir / "analysis.md").write_text("\n".join(report_lines), encoding="utf-8")
    print(f"analysis={experiment_dir / 'analysis.md'}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as error:
        print(f"analysis failed: {error}", file=sys.stderr)
        raise SystemExit(1)
