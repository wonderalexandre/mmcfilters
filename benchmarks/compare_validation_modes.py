#!/usr/bin/env python3
"""Compare CHECKED and UNCHECKED scientific benchmark executables."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import shlex
import statistics
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


MEDIAN_SUFFIX = "_median_ms"
CHECKSUM_SUFFIX = "_checksum"
SAMPLES_SUFFIX = "_samples_ms"
INVARIANT_SUFFIXES = (
    "_checksum",
    "_pixels",
    "_proper_parts",
    "_primary_node_slots",
    "_primary_live_nodes",
    "_secondary_node_slots",
    "_secondary_live_nodes",
    "_edges",
    "_steps",
    "_primary_nodes_removed",
    "_secondary_nodes_removed",
    "_complete_validation_commits",
    "_incremental_validation_commits",
    "_light_threshold",
    "_medium_threshold",
    "_heavy_threshold",
)


SCENARIO_JSON_FIELDS = (
    "median_ms",
    "mad_ms",
    "minimum_ms",
    "maximum_ms",
    "checksum",
    "pixels",
    "proper_parts",
    "primary_node_slots",
    "primary_live_nodes",
    "secondary_node_slots",
    "secondary_live_nodes",
    "edges",
    "steps",
    "primary_nodes_removed",
    "secondary_nodes_removed",
    "complete_validation_commits",
    "incremental_validation_commits",
    "light_threshold",
    "medium_threshold",
    "heavy_threshold",
)


def scalar_text(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def parse_output(output: str) -> dict[str, str]:
    stripped_lines = [line.strip() for line in output.splitlines() if line.strip()]
    if stripped_lines and stripped_lines[0].startswith("{"):
        parsed: dict[str, str] = {}
        for line in stripped_lines:
            record = json.loads(line)
            if record.get("record") == "metadata":
                for key, value in record.items():
                    if key != "record":
                        parsed[key] = scalar_text(value)
            elif record.get("record") == "scenario":
                name = record.get("name")
                if not isinstance(name, str) or not name:
                    raise RuntimeError("benchmark JSONL scenario omitted its name")
                for field in SCENARIO_JSON_FIELDS:
                    if field in record:
                        parsed[f"{name}_{field}"] = scalar_text(record[field])
                if "samples_ms" in record:
                    parsed[f"{name}{SAMPLES_SUFFIX}"] = json.dumps(record["samples_ms"], separators=(",", ":"))
        return parsed

    parsed: dict[str, str] = {}
    for line in output.splitlines():
        key, separator, value = line.partition("=")
        if separator:
            parsed[key.strip()] = value.strip()
    return parsed


def run_benchmark(
    executable: Path,
    mode: str,
    benchmark_arguments: list[str],
    expected_legacy_dimensions: tuple[int, int, int] | None,
) -> tuple[dict[str, str], str]:
    completed = subprocess.run(
        [str(executable), *benchmark_arguments],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise RuntimeError(f"{mode} benchmark failed: {detail}")

    parsed = parse_output(completed.stdout)
    if parsed.get("contract_mode") != mode:
        actual = parsed.get("contract_mode", "missing")
        raise RuntimeError(f"expected {mode} executable, received contract_mode={actual}")
    if expected_legacy_dimensions is not None:
        rows, cols, repetitions = expected_legacy_dimensions
        for key, expected in (("rows", rows), ("cols", cols), ("repetitions", repetitions)):
            if parsed.get(key) != str(expected):
                raise RuntimeError(f"{mode} benchmark reported inconsistent {key}")
    return parsed, completed.stdout


def scenario_names(result: dict[str, str]) -> list[str]:
    return sorted(key[: -len(MEDIAN_SUFFIX)] for key in result if key.endswith(MEDIAN_SUFFIX))


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Alternate CHECKED/UNCHECKED runs, verify result hashes, and report median-of-medians timings.",
        epilog="Place a standalone -- after the two executables to forward arbitrary arguments to both benchmarks.",
    )
    parser.add_argument("checked", type=Path, help="path to the CHECKED benchmark executable")
    parser.add_argument("unchecked", type=Path, help="path to the UNCHECKED benchmark executable")
    parser.add_argument("--rows", type=int, default=256)
    parser.add_argument("--cols", type=int, default=256)
    parser.add_argument("--repetitions", type=int, default=9, help="timed samples within each process")
    parser.add_argument("--process-runs", type=int, default=5, help="number of process invocations per contract mode")
    parser.add_argument("--output-dir", type=Path, help="persist raw runs, reproducibility metadata, JSONL, and CSV summaries")
    parser.add_argument("--experiment-name", default="validation-contract", help="stable label stored with persisted results")
    parser.add_argument("--capture-samples", action="store_true", help="request and persist every timed repetition")
    raw_arguments = sys.argv[1:]
    benchmark_arguments: list[str] = []
    if "--" in raw_arguments:
        separator = raw_arguments.index("--")
        benchmark_arguments = raw_arguments[separator + 1 :]
        raw_arguments = raw_arguments[:separator]
    arguments = parser.parse_args(raw_arguments)
    if arguments.rows <= 0 or arguments.cols <= 0 or arguments.repetitions <= 0 or arguments.process_runs <= 0:
        parser.error("rows, cols, repetitions, and process-runs must be positive")
    arguments.benchmark_arguments = benchmark_arguments
    if arguments.capture_samples and not benchmark_arguments:
        parser.error("--capture-samples requires forwarded benchmark arguments after --")
    return arguments


def command_output(command: list[str], cwd: Path | None = None) -> str:
    completed = subprocess.run(command, cwd=cwd, check=False, capture_output=True, text=True)
    return completed.stdout.strip() if completed.returncode == 0 else "unavailable"


def cpu_name() -> str:
    if sys.platform == "darwin":
        value = command_output(["sysctl", "-n", "machdep.cpu.brand_string"])
        if value != "unavailable":
            return value
    value = platform.processor().strip()
    if value:
        return value
    if sys.platform.startswith("linux"):
        try:
            for line in Path("/proc/cpuinfo").read_text(encoding="utf-8").splitlines():
                if line.lower().startswith("model name"):
                    return line.partition(":")[2].strip()
        except OSError:
            pass
    return platform.machine() or "unavailable"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def find_cmake_cache(executable: Path) -> Path | None:
    for directory in (executable.parent, *executable.parents):
        candidate = directory / "CMakeCache.txt"
        if candidate.is_file():
            return candidate
    return None


def cmake_metadata(executable: Path) -> dict[str, str]:
    cache_path = find_cmake_cache(executable)
    if cache_path is None:
        return {"cache": "unavailable"}
    selected_prefixes = (
        "CMAKE_BUILD_TYPE",
        "CMAKE_CXX_COMPILER:",
        "CMAKE_CXX_COMPILER_VERSION",
        "CMAKE_CXX_FLAGS:",
        "CMAKE_CXX_FLAGS_DEBUG",
        "CMAKE_CXX_FLAGS_RELEASE",
        "CMAKE_CXX_FLAGS_RELWITHDEBINFO",
        "CMAKE_GENERATOR:",
        "MMCFILTERS_CONTRACT_MODE",
    )
    metadata = {"cache": str(cache_path.resolve())}
    for line in cache_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("//") or line.startswith("#") or "=" not in line:
            continue
        declaration, value = line.split("=", 1)
        key = declaration.split(":", 1)[0]
        if declaration.startswith(selected_prefixes):
            metadata[key] = value
    build_type = metadata.get("CMAKE_BUILD_TYPE", "")
    flags = [metadata.get("CMAKE_CXX_FLAGS", "")]
    if build_type:
        flags.append(metadata.get(f"CMAKE_CXX_FLAGS_{build_type.upper()}", ""))
    metadata["effective_cxx_flags"] = " ".join(item for item in flags if item).strip()
    return metadata


def git_metadata(repository: Path) -> dict[str, object]:
    commit = command_output(["git", "rev-parse", "HEAD"], repository)
    branch = command_output(["git", "branch", "--show-current"], repository)
    status = command_output(["git", "status", "--porcelain=v1"], repository)
    return {
        "repository": str(repository),
        "commit": commit,
        "branch": branch,
        "dirty": status not in ("", "unavailable"),
        "status": status,
    }


def prepare_output_directory(path: Path | None) -> Path | None:
    if path is None:
        return None
    resolved = path.resolve()
    if resolved.exists() and any(resolved.iterdir()):
        raise RuntimeError(f"output directory is not empty: {resolved}")
    resolved.mkdir(parents=True, exist_ok=True)
    (resolved / "raw").mkdir()
    return resolved


def parse_samples(result: dict[str, str], scenario: str) -> list[float] | None:
    raw = result.get(scenario + SAMPLES_SUFFIX)
    if raw is None:
        return None
    try:
        values = json.loads(raw) if raw.startswith("[") else [float(item) for item in raw.split(",") if item]
    except (ValueError, json.JSONDecodeError) as error:
        raise RuntimeError(f"{scenario} reported malformed timed samples") from error
    if not isinstance(values, list) or not values:
        raise RuntimeError(f"{scenario} reported no timed samples")
    samples = [float(value) for value in values]
    if any(not math.isfinite(value) or value < 0.0 for value in samples):
        raise RuntimeError(f"{scenario} reported an invalid timed sample")
    return samples


def validate_sample_summary(result: dict[str, str], scenario: str, samples: list[float], repetitions: int) -> None:
    if len(samples) != repetitions:
        raise RuntimeError(f"{scenario} reported {len(samples)} samples for {repetitions} repetitions")
    center = statistics.median(samples)
    mad = statistics.median(abs(value - center) for value in samples)
    expected = {
        "median_ms": center,
        "mad_ms": mad,
        "minimum_ms": min(samples),
        "maximum_ms": max(samples),
    }
    for field, value in expected.items():
        reported = float(result[scenario + "_" + field])
        if not math.isclose(reported, value, rel_tol=1e-6, abs_tol=2e-6):
            raise RuntimeError(f"{scenario} {field} is inconsistent with its timed samples")


def scenario_record(result: dict[str, str], scenario: str, samples: list[float] | None = None) -> dict[str, object]:
    record: dict[str, object] = {"record": "scenario", "name": scenario}
    for field in SCENARIO_JSON_FIELDS:
        key = scenario + "_" + field
        if key in result:
            value = result[key]
            if field == "checksum":
                record[field] = value
            elif field.endswith("_ms") or field.endswith("threshold"):
                record[field] = float(value)
            else:
                record[field] = int(value)
    parts = scenario.split(".", 2)
    if len(parts) == 3:
        record["suite"], record["scope"], record["scenario"] = parts
    if samples is not None:
        record["samples_ms"] = samples
    return record


def write_json_lines(path: Path, records: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        for record in records:
            stream.write(json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n")


def write_summary_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError("cannot persist an empty benchmark summary")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_samples_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError("cannot persist an empty timed-sample table")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    arguments = parse_arguments()
    executables = {"CHECKED": arguments.checked.resolve(), "UNCHECKED": arguments.unchecked.resolve()}
    for mode, executable in executables.items():
        if not executable.is_file():
            raise RuntimeError(f"{mode} executable does not exist: {executable}")
    output_directory = prepare_output_directory(arguments.output_dir)
    repository = Path(__file__).resolve().parent.parent
    experiment_metadata: dict[str, object] = {
        "record": "experiment",
        "schema_version": 2,
        "name": arguments.experiment_name,
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "working_directory": os.getcwd(),
        "command": [sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]],
        "host": {
            "cpu": cpu_name(),
            "logical_cpu_count": os.cpu_count(),
            "machine": platform.machine(),
            "operating_system": platform.platform(),
            "python": platform.python_version(),
        },
        "git": git_metadata(repository),
        "executables": {
            mode: {
                "path": str(executable),
                "sha256": sha256_file(executable),
                "cmake": cmake_metadata(executable),
            }
            for mode, executable in executables.items()
        },
    }

    timings: dict[str, dict[str, list[float]]] = {
        "CHECKED": defaultdict(list),
        "UNCHECKED": defaultdict(list),
    }
    normalized_records: list[dict[str, object]] = [experiment_metadata]
    sample_records: list[dict[str, object]] = []
    samples_by_mode: dict[str, dict[str, list[list[float]]]] = {
        "CHECKED": defaultdict(list),
        "UNCHECKED": defaultdict(list),
    }
    checksums: dict[str, str] = {}
    invariant_values: dict[str, str] = {}
    expected_scenarios: list[str] | None = None
    expected_metadata: dict[str, str] | None = None
    if arguments.benchmark_arguments:
        benchmark_arguments = list(arguments.benchmark_arguments)
        if arguments.capture_samples and "--emit-samples" not in benchmark_arguments:
            benchmark_arguments.append("--emit-samples")
        expected_legacy_dimensions = None
    else:
        benchmark_arguments = [str(arguments.rows), str(arguments.cols), str(arguments.repetitions)]
        expected_legacy_dimensions = (arguments.rows, arguments.cols, arguments.repetitions)

    for process_run in range(arguments.process_runs):
        order = ("CHECKED", "UNCHECKED") if process_run % 2 == 0 else ("UNCHECKED", "CHECKED")
        for order_index, mode in enumerate(order):
            result, raw_output = run_benchmark(executables[mode], mode, benchmark_arguments, expected_legacy_dimensions)
            scenarios = scenario_names(result)
            if not scenarios:
                raise RuntimeError(f"{mode} benchmark reported no timed scenarios")
            if expected_scenarios is None:
                expected_scenarios = scenarios
            elif scenarios != expected_scenarios:
                raise RuntimeError(f"{mode} benchmark reported a different scenario set")

            metadata_keys = (
                "benchmark",
                "profile",
                "input",
                "input_source",
                "input_checksum",
                "workload",
                "manifest",
                "casf_quantiles",
                "attribute_bundles",
                "suite",
                "rows",
                "cols",
                "repetitions",
                "scenario_count",
                "samples_emitted",
                "cplusplus",
                "compiler",
            )
            metadata = {key: result[key] for key in metadata_keys if key in result}
            if expected_metadata is None:
                expected_metadata = metadata
            elif metadata != expected_metadata:
                raise RuntimeError(f"{mode} benchmark reported inconsistent run metadata")

            raw_name = f"process-{process_run + 1:03d}-order-{order_index + 1}-{mode.lower()}.txt"
            if output_directory is not None:
                (output_directory / "raw" / raw_name).write_text(raw_output, encoding="utf-8")
            normalized_records.append(
                {
                    "record": "process",
                    "process_run": process_run + 1,
                    "order": order_index + 1,
                    "contract_mode": mode,
                    "raw_output": f"raw/{raw_name}",
                    "benchmark_metadata": metadata,
                }
            )

            for scenario in scenarios:
                checksum_key = scenario + CHECKSUM_SUFFIX
                if checksum_key not in result:
                    raise RuntimeError(f"{mode} benchmark omitted {checksum_key}")
                timings[mode][scenario].append(float(result[scenario + MEDIAN_SUFFIX]))
                checksum = result[checksum_key]
                previous = checksums.setdefault(scenario, checksum)
                if checksum != previous:
                    raise RuntimeError(
                        f"{scenario} checksum mismatch: expected {previous}, received {checksum} from {mode}"
                    )
                samples = parse_samples(result, scenario)
                samples_expected = arguments.capture_samples or result.get("samples_emitted") == "true"
                if samples_expected and samples is None:
                    raise RuntimeError(f"{mode} benchmark omitted timed samples for {scenario}")
                if samples is not None:
                    repetitions = int(result["repetitions"])
                    validate_sample_summary(result, scenario, samples, repetitions)
                    samples_by_mode[mode][scenario].append(samples)
                    parts = scenario.split(".", 2)
                    for sample_index, sample in enumerate(samples, start=1):
                        sample_records.append(
                            {
                                "record": "sample",
                                "process_run": process_run + 1,
                                "order": order_index + 1,
                                "contract_mode": mode,
                                "scenario": scenario,
                                "suite": parts[0] if len(parts) == 3 else "",
                                "scope": parts[1] if len(parts) == 3 else "",
                                "name": parts[2] if len(parts) == 3 else scenario,
                                "sample_index": sample_index,
                                "sample_ms": sample,
                                "checksum": checksum,
                            }
                        )
                record = scenario_record(result, scenario, samples)
                record.update({"process_run": process_run + 1, "order": order_index + 1, "contract_mode": mode})
                normalized_records.append(record)
            for key, value in result.items():
                if key.endswith(INVARIANT_SUFFIXES):
                    previous = invariant_values.setdefault(key, value)
                    if value != previous:
                        raise RuntimeError(f"{key} mismatch: expected {previous}, received {value} from {mode}")

    assert expected_scenarios is not None
    if arguments.benchmark_arguments:
        print(f"benchmark_arguments={shlex.join(benchmark_arguments)}")
        assert expected_metadata is not None
        for key, value in expected_metadata.items():
            print(f"{key}={value}")
    else:
        print(f"rows={arguments.rows}")
        print(f"cols={arguments.cols}")
        print(f"repetitions_per_process={arguments.repetitions}")
    print(f"process_runs_per_mode={arguments.process_runs}")
    print("scenario\tchecked_median_ms\tunchecked_median_ms\tunchecked_vs_checked_percent\tchecksum")
    summary_rows: list[dict[str, object]] = []
    for scenario in expected_scenarios:
        checked_median = statistics.median(timings["CHECKED"][scenario])
        unchecked_median = statistics.median(timings["UNCHECKED"][scenario])
        relative_percent = 100.0 * (unchecked_median / checked_median - 1.0)
        parts = scenario.split(".", 2)
        summary_rows.append(
            {
                "scenario": scenario,
                "suite": parts[0] if len(parts) == 3 else "",
                "scope": parts[1] if len(parts) == 3 else "",
                "name": parts[2] if len(parts) == 3 else scenario,
                "checked_median_ms": checked_median,
                "unchecked_median_ms": unchecked_median,
                "unchecked_minus_checked_ms": unchecked_median - checked_median,
                "unchecked_vs_checked_percent": relative_percent,
                "checked_process_medians_ms": json.dumps(timings["CHECKED"][scenario], separators=(",", ":")),
                "unchecked_process_medians_ms": json.dumps(timings["UNCHECKED"][scenario], separators=(",", ":")),
                "checked_sample_count": sum(len(values) for values in samples_by_mode["CHECKED"].get(scenario, [])),
                "unchecked_sample_count": sum(len(values) for values in samples_by_mode["UNCHECKED"].get(scenario, [])),
                "checksum": checksums[scenario],
            }
        )
        print(
            f"{scenario}\t{checked_median:.3f}\t{unchecked_median:.3f}\t"
            f"{relative_percent:+.2f}\t{checksums[scenario]}"
        )

    if output_directory is not None:
        experiment_metadata["completed_at_utc"] = datetime.now(timezone.utc).isoformat()
        experiment_metadata["process_runs_per_mode"] = arguments.process_runs
        experiment_metadata["benchmark_arguments"] = benchmark_arguments
        experiment_metadata["benchmark_metadata"] = expected_metadata or {}
        experiment_metadata["samples_captured"] = bool(sample_records)
        (output_directory / "experiment.json").write_text(
            json.dumps(experiment_metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        write_json_lines(output_directory / "raw.jsonl", normalized_records)
        write_json_lines(
            output_directory / "summary.jsonl",
            [{**experiment_metadata, "record": "experiment_summary"}, *({"record": "summary", **row} for row in summary_rows)],
        )
        write_summary_csv(output_directory / "summary.csv", summary_rows)
        if sample_records:
            write_json_lines(output_directory / "samples.jsonl", sample_records)
            write_samples_csv(output_directory / "samples.csv", sample_records)
        print(f"results_directory={output_directory}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except RuntimeError as error:
        print(f"comparison failed: {error}", file=sys.stderr)
        sys.exit(1)
