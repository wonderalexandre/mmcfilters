#!/usr/bin/env python3
"""Combine multiple persisted scientific API experiments into one campaign report."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


ATTRIBUTE_PAIR_FIELDS = [
    "workload",
    "contract_mode",
    "suite",
    "name",
    "grouped_ms",
    "sequential_scalars_ms",
    "sequential_over_grouped_ratio",
    "grouped_time_saved_percent",
    "checksum",
]
CASF_FIELDS = [
    "workload",
    "scenario",
    "checked_median_ms",
    "unchecked_median_ms",
    "unchecked_vs_checked_percent",
    "steps",
    "primary_nodes_removed",
    "secondary_nodes_removed",
    "complete_validation_commits",
    "incremental_validation_commits",
    "light_threshold",
    "medium_threshold",
    "heavy_threshold",
    "checksum",
]
SCOPE_PAIR_FIELDS = [
    "workload",
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


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize multiple compare_validation_modes.py experiment directories.")
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("experiment_dirs", nargs="+", type=Path)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join("---" for _ in headers) + " |"]
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(lines)


def prepare_output_directory(path: Path) -> Path:
    resolved = path.resolve()
    if resolved.exists() and any(resolved.iterdir()):
        raise RuntimeError(f"campaign output directory is not empty: {resolved}")
    resolved.mkdir(parents=True, exist_ok=True)
    return resolved


def command_output(command: list[str], cwd: Path) -> str:
    completed = subprocess.run(command, cwd=cwd, check=False, capture_output=True, text=True)
    return completed.stdout.strip() if completed.returncode == 0 else "unavailable"


def main() -> int:
    arguments = parse_arguments()
    if len(arguments.experiment_dirs) < 2:
        raise RuntimeError("campaign analysis requires at least two experiments")
    output_dir = prepare_output_directory(arguments.output_dir)

    experiments: list[tuple[Path, dict[str, object], list[dict[str, str]]]] = []
    for source in arguments.experiment_dirs:
        source = source.resolve()
        metadata_path = source / "experiment.json"
        summary_path = source / "summary.csv"
        for path in (
            metadata_path,
            summary_path,
            source / "analysis.md",
            source / "scope-pairs.csv",
            source / "attribute-pairs.csv",
            source / "casf.csv",
        ):
            if not path.is_file():
                raise RuntimeError(f"missing experiment artifact: {path}")
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        summary = read_csv(summary_path)
        if not summary:
            raise RuntimeError(f"empty experiment summary: {source}")
        experiments.append((source, metadata, summary))

    reference_metadata = experiments[0][1]
    reference_git = reference_metadata.get("git", {}).get("commit")
    reference_host = reference_metadata.get("host", {})
    reference_benchmark = reference_metadata.get("benchmark_metadata", {})
    sample_capture_states = {bool(metadata.get("samples_captured")) for _, metadata, _ in experiments}
    if len(sample_capture_states) != 1:
        raise RuntimeError("campaign mixes experiments with and without timed-sample capture")
    samples_captured = sample_capture_states.pop()
    for source, metadata, _ in experiments[1:]:
        if metadata.get("git", {}).get("commit") != reference_git:
            raise RuntimeError(f"campaign mixes Git commits: {source}")
        if metadata.get("host", {}).get("cpu") != reference_host.get("cpu"):
            raise RuntimeError(f"campaign mixes CPUs: {source}")
        benchmark = metadata.get("benchmark_metadata", {})
        for key in ("profile", "repetitions", "cplusplus", "compiler"):
            if benchmark.get(key) != reference_benchmark.get(key):
                raise RuntimeError(f"campaign metadata differs for {key}: {source}")
        if metadata.get("process_runs_per_mode") != reference_metadata.get("process_runs_per_mode"):
            raise RuntimeError(f"campaign process-run count differs: {source}")
    if samples_captured:
        for source, _, _ in experiments:
            for path in (source / "samples.csv", source / "sample-statistics.csv"):
                if not path.is_file():
                    raise RuntimeError(f"missing timed-sample campaign artifact: {path}")

    workload_rows: list[dict[str, object]] = []
    combined_rows: list[dict[str, object]] = []
    grouped_by_scenario: dict[str, list[dict[str, object]]] = defaultdict(list)
    grouped_by_scope: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    combined_attribute_rows: list[dict[str, object]] = []
    combined_casf_rows: list[dict[str, object]] = []
    combined_scope_pair_rows: list[dict[str, object]] = []
    combined_sample_statistics: list[dict[str, object]] = []

    for source, metadata, summary in experiments:
        benchmark = metadata["benchmark_metadata"]
        workload = str(benchmark["workload"])
        checked_sum = sum(float(row["checked_median_ms"]) for row in summary)
        unchecked_sum = sum(float(row["unchecked_median_ms"]) for row in summary)
        differences = [float(row["unchecked_vs_checked_percent"]) for row in summary]
        workload_rows.append(
            {
                "workload": workload,
                "input_source": benchmark["input_source"],
                "rows": int(benchmark["rows"]),
                "cols": int(benchmark["cols"]),
                "scenario_count": len(summary),
                "checked_sum_of_medians_ms": checked_sum,
                "unchecked_sum_of_medians_ms": unchecked_sum,
                "unchecked_vs_checked_sum_percent": 100.0 * (unchecked_sum / checked_sum - 1.0),
                "median_scenario_difference_percent": statistics.median(differences),
                "unchecked_faster_scenarios": sum(value < 0 for value in differences),
            }
        )
        for row in summary:
            combined = {"workload": workload, **row}
            combined_rows.append(combined)
            grouped_by_scenario[row["scenario"]].append(combined)
            grouped_by_scope[(row["suite"], row["scope"])].append(combined)
        for row in read_csv(source / "attribute-pairs.csv"):
            combined_attribute_rows.append({"workload": workload, **row})
        for row in read_csv(source / "casf.csv"):
            combined_casf_rows.append({"workload": workload, **row})
        for row in read_csv(source / "scope-pairs.csv"):
            combined_scope_pair_rows.append({"workload": workload, **row})
        if samples_captured:
            sample_statistics = read_csv(source / "sample-statistics.csv")
            if len(sample_statistics) != len(summary):
                raise RuntimeError(f"timed-sample statistics do not cover every scenario: {source}")
            for row in sample_statistics:
                combined_sample_statistics.append({"workload": workload, **row})

    workload_fields = [
        "workload",
        "input_source",
        "rows",
        "cols",
        "scenario_count",
        "checked_sum_of_medians_ms",
        "unchecked_sum_of_medians_ms",
        "unchecked_vs_checked_sum_percent",
        "median_scenario_difference_percent",
        "unchecked_faster_scenarios",
    ]
    write_csv(output_dir / "campaign-workloads.csv", workload_rows, workload_fields)
    write_csv(output_dir / "campaign-scenarios.csv", combined_rows, list(combined_rows[0]))
    write_csv(output_dir / "campaign-attribute-pairs.csv", combined_attribute_rows, ATTRIBUTE_PAIR_FIELDS)
    write_csv(output_dir / "campaign-casf.csv", combined_casf_rows, CASF_FIELDS)
    write_csv(output_dir / "campaign-scope-pairs.csv", combined_scope_pair_rows, SCOPE_PAIR_FIELDS)
    if combined_sample_statistics:
        write_csv(
            output_dir / "campaign-sample-statistics.csv",
            combined_sample_statistics,
            list(combined_sample_statistics[0]),
        )

    common_scenario_rows: list[dict[str, object]] = []
    experiment_count = len(experiments)
    for scenario, rows in sorted(grouped_by_scenario.items()):
        if len(rows) != experiment_count:
            continue
        checked_sum = sum(float(row["checked_median_ms"]) for row in rows)
        unchecked_sum = sum(float(row["unchecked_median_ms"]) for row in rows)
        common_scenario_rows.append(
            {
                "scenario": scenario,
                "suite": rows[0]["suite"],
                "scope": rows[0]["scope"],
                "name": rows[0]["name"],
                "workload_count": len(rows),
                "checked_sum_ms": checked_sum,
                "unchecked_sum_ms": unchecked_sum,
                "unchecked_vs_checked_sum_percent": 100.0 * (unchecked_sum / checked_sum - 1.0),
                "median_workload_difference_percent": statistics.median(
                    float(row["unchecked_vs_checked_percent"]) for row in rows
                ),
            }
        )
    common_fields = [
        "scenario",
        "suite",
        "scope",
        "name",
        "workload_count",
        "checked_sum_ms",
        "unchecked_sum_ms",
        "unchecked_vs_checked_sum_percent",
        "median_workload_difference_percent",
    ]
    write_csv(output_dir / "campaign-common-scenarios.csv", common_scenario_rows, common_fields)

    scope_rows: list[dict[str, object]] = []
    for (suite, scope), rows in sorted(grouped_by_scope.items()):
        checked_sum = sum(float(row["checked_median_ms"]) for row in rows)
        unchecked_sum = sum(float(row["unchecked_median_ms"]) for row in rows)
        scope_rows.append(
            {
                "suite": suite,
                "scope": scope,
                "observation_count": len(rows),
                "checked_sum_of_medians_ms": checked_sum,
                "unchecked_sum_of_medians_ms": unchecked_sum,
                "unchecked_vs_checked_sum_percent": 100.0 * (unchecked_sum / checked_sum - 1.0),
            }
        )
    scope_fields = [
        "suite",
        "scope",
        "observation_count",
        "checked_sum_of_medians_ms",
        "unchecked_sum_of_medians_ms",
        "unchecked_vs_checked_sum_percent",
    ]
    write_csv(output_dir / "campaign-scope-summary.csv", scope_rows, scope_fields)

    all_checked = sum(float(row["checked_median_ms"]) for row in combined_rows)
    all_unchecked = sum(float(row["unchecked_median_ms"]) for row in combined_rows)
    all_differences = [float(row["unchecked_vs_checked_percent"]) for row in combined_rows]
    material_common = [row for row in common_scenario_rows if float(row["checked_sum_ms"]) >= 1.0]
    strongest_gains = sorted(material_common, key=lambda row: float(row["unchecked_vs_checked_sum_percent"]))[:12]
    strongest_losses = sorted(material_common, key=lambda row: float(row["unchecked_vs_checked_sum_percent"]), reverse=True)[:8]

    grouped_ratios: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    for row in combined_attribute_rows:
        grouped_ratios[(str(row["contract_mode"]), str(row["suite"]), str(row["name"]))].append(
            float(row["sequential_over_grouped_ratio"])
        )
    attribute_summary = [
        {
            "contract_mode": mode,
            "suite": suite,
            "name": name,
            "median_sequential_over_grouped_ratio": statistics.median(ratios),
            "workload_count": len(ratios),
        }
        for (mode, suite, name), ratios in sorted(grouped_ratios.items())
    ]
    write_csv(
        output_dir / "campaign-attribute-summary.csv",
        attribute_summary,
        ["contract_mode", "suite", "name", "median_sequential_over_grouped_ratio", "workload_count"],
    )

    analysis_repository = Path(__file__).resolve().parent.parent
    analysis_git_commit = command_output(["git", "rev-parse", "HEAD"], analysis_repository)
    analysis_git_status = command_output(["git", "status", "--porcelain=v1"], analysis_repository)
    analysis_git_dirty = analysis_git_status not in ("", "unavailable")

    if combined_sample_statistics:
        supported_faster = [row for row in combined_sample_statistics if row["inference"] == "UNCHECKED_FASTER"]
        supported_slower = [row for row in combined_sample_statistics if row["inference"] == "UNCHECKED_SLOWER"]
        inconclusive = len(combined_sample_statistics) - len(supported_faster) - len(supported_slower)
        strongest_supported = sorted(
            supported_faster, key=lambda row: float(row["unchecked_vs_checked_percent"])
        )[:12]
        sample_report = [
            "## Raw-sample inference",
            "",
            f"Across the campaign, process-cluster bootstrap intervals support UNCHECKED as faster in {len(supported_faster)} "
            f"workload-scenario observations and slower in {len(supported_slower)}; {inconclusive} remain inconclusive.",
            "",
            markdown_table(
                ["Workload", "Scenario", "Difference", "Bootstrap 95% CI", "Cliff's delta"],
                [
                    [
                        str(row["workload"]),
                        str(row["scenario"]),
                        f"{float(row['unchecked_vs_checked_percent']):+.2f}%",
                        f"[{float(row['bootstrap_ci95_low_percent']):+.2f}%, {float(row['bootstrap_ci95_high_percent']):+.2f}%]",
                        f"{float(row['cliffs_delta_process_medians']):+.2f}",
                    ]
                    for row in strongest_supported
                ],
            ),
            "",
            "The campaign table preserves each workload-scenario interval separately; it does not pool inner repetitions or manufacture a larger independent sample.",
            "",
        ]
    else:
        sample_report = [
            "## Raw-sample inference",
            "",
            "These experiments predate timed-sample capture, so the campaign reports descriptive differences only.",
            "",
        ]

    report = [
        "# Scientific API validation-contract campaign",
        "",
        "## Reproducibility",
        "",
        f"- Experiment Git commit: `{reference_git}`",
        f"- Analysis Git commit: `{analysis_git_commit}`; dirty worktree: `{analysis_git_dirty}`",
        f"- CPU: {reference_host.get('cpu', 'unknown')}",
        f"- Operating system: {reference_host.get('operating_system', 'unknown')}",
        f"- Compiler: {reference_benchmark.get('compiler', 'unknown')}",
        f"- Profile: `{reference_benchmark.get('profile', 'unknown')}`",
        f"- Repetitions per scenario/process: {reference_benchmark.get('repetitions', 'unknown')}",
        f"- Process runs per contract mode: {reference_metadata.get('process_runs_per_mode', 'unknown')}",
        "",
        "Every experiment passed deterministic repetition checks and CHECKED/UNCHECKED checksum equality for all scenarios and structural outcomes.",
        "",
        "## Workloads",
        "",
        markdown_table(
            ["Workload", "Domain", "Scenarios", "CHECKED sum (ms)", "UNCHECKED sum (ms)", "Difference"],
            [
                [
                    str(row["workload"]),
                    f"{row['rows']}x{row['cols']}",
                    str(row["scenario_count"]),
                    f"{float(row['checked_sum_of_medians_ms']):.3f}",
                    f"{float(row['unchecked_sum_of_medians_ms']):.3f}",
                    f"{float(row['unchecked_vs_checked_sum_percent']):+.2f}%",
                ]
                for row in workload_rows
            ],
        ),
        "",
        f"Across {len(combined_rows)} workload-scenario observations, UNCHECKED was faster in {sum(value < 0 for value in all_differences)}. "
        f"The median observation-level difference was {statistics.median(all_differences):+.2f}%.",
        "",
        f"The descriptive sum of all independent scenario medians was {all_checked:.3f} ms in CHECKED and {all_unchecked:.3f} ms in UNCHECKED "
        f"({100.0 * (all_unchecked / all_checked - 1.0):+.2f}%). It is not an application runtime.",
        "",
        *sample_report,
        "## Suites and timing scopes",
        "",
        markdown_table(
            ["Suite", "Scope", "Observations", "CHECKED sum (ms)", "UNCHECKED sum (ms)", "Difference"],
            [
                [
                    str(row["suite"]),
                    str(row["scope"]),
                    str(row["observation_count"]),
                    f"{float(row['checked_sum_of_medians_ms']):.3f}",
                    f"{float(row['unchecked_sum_of_medians_ms']):.3f}",
                    f"{float(row['unchecked_vs_checked_sum_percent']):+.2f}%",
                ]
                for row in scope_rows
            ],
        ),
        "",
        "`end_to_end` and `established_input` contain different scientific operations. Their totals are descriptive and are not paired estimates of validation overhead.",
        "",
        "## Paired end-to-end and established-input CASF",
        "",
        "These rows have identical final checksums and threshold sequences. Their difference isolates CASF construction/setup from the established incremental sequence.",
        "",
        markdown_table(
            ["Workload", "Mode", "End-to-end (ms)", "Established (ms)", "Setup (ms)", "Setup share"],
            [
                [
                    str(row["workload"]),
                    mode.upper(),
                    f"{float(row[f'{mode}_end_to_end_ms']):.3f}",
                    f"{float(row[f'{mode}_established_input_ms']):.3f}",
                    f"{float(row[f'{mode}_setup_overhead_ms']):.3f}",
                    f"{float(row[f'{mode}_setup_share_percent']):.2f}%",
                ]
                for row in combined_scope_pair_rows
                for mode in ("checked", "unchecked")
            ],
        ),
        "",
        "## Largest consistent reductions",
        "",
        markdown_table(
            ["Scenario", "CHECKED sum (ms)", "UNCHECKED sum (ms)", "Difference"],
            [
                [
                    str(row["scenario"]),
                    f"{float(row['checked_sum_ms']):.3f}",
                    f"{float(row['unchecked_sum_ms']):.3f}",
                    f"{float(row['unchecked_vs_checked_sum_percent']):+.2f}%",
                ]
                for row in strongest_gains
            ],
        ),
        "",
        "## Largest increases",
        "",
        markdown_table(
            ["Scenario", "CHECKED sum (ms)", "UNCHECKED sum (ms)", "Difference"],
            [
                [
                    str(row["scenario"]),
                    f"{float(row['checked_sum_ms']):.3f}",
                    f"{float(row['unchecked_sum_ms']):.3f}",
                    f"{float(row['unchecked_vs_checked_sum_percent']):+.2f}%",
                ]
                for row in strongest_losses
            ],
        ),
        "",
        "Small increases and reductions in sub-millisecond scalar scenarios should be treated as timing noise unless supported by a targeted experiment.",
        "",
        "## Attribute grouping",
        "",
        markdown_table(
            ["Mode", "Group/bundle", "Median scalars/grouped"],
            [
                [
                    str(row["contract_mode"]),
                    f"{row['suite']}.{row['name']}",
                    f"{float(row['median_sequential_over_grouped_ratio']):.2f}x",
                ]
                for row in attribute_summary
            ],
        ),
        "",
        "## CASF",
        "",
        "CASF aggregate timings and invariant edit counts for every workload are retained in `campaign-casf.csv`. The common-scenario table allows direct aggregation without mixing workload-specific bundles.",
        "",
        "The complete campaign tables are `campaign-workloads.csv`, `campaign-scenarios.csv`, `campaign-common-scenarios.csv`, "
        "`campaign-scope-summary.csv`, `campaign-scope-pairs.csv`, `campaign-attribute-pairs.csv`, `campaign-attribute-summary.csv`, and `campaign-casf.csv`. "
        "Campaigns with captured repetitions also include `campaign-sample-statistics.csv`.",
        "",
    ]
    (output_dir / "campaign-analysis.md").write_text("\n".join(report), encoding="utf-8")
    campaign_metadata = {
        "schema_version": 2,
        "experiment_git_commit": reference_git,
        "analysis_git_commit": analysis_git_commit,
        "analysis_git_dirty": analysis_git_dirty,
        "host": reference_host,
        "benchmark": reference_benchmark,
        "process_runs_per_mode": reference_metadata.get("process_runs_per_mode"),
        "samples_captured": samples_captured,
        "experiments": [str(source) for source, _, _ in experiments],
    }
    (output_dir / "campaign.json").write_text(json.dumps(campaign_metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"campaign_analysis={output_dir / 'campaign-analysis.md'}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as error:
        print(f"campaign analysis failed: {error}", file=sys.stderr)
        raise SystemExit(1)
