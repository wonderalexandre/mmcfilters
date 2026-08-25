#!/usr/bin/env python3
"""Run the complete SIBGRAPI 2026 bitquad timing campaign."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import platform
import struct
import subprocess
import sys
from pathlib import Path
from typing import Optional


REFERENCE5_COMMIT = "da32cf7666a774c25d11dc0200a63ebb3f1fe574"
BASE_RELEASE_TAG = "v4.3.0"
RESOLUTIONS = (
    ("480p", 853, 480),
    ("720p", 1280, 720),
    ("1080p", 1920, 1080),
)
REQUIRED_REFERENCE_FILES = (
    "include/morphotree/adjacency/adjacency8c.hpp",
    "include/morphotree/attributes/bitquads/quadCountComputer.hpp",
    "include/morphotree/attributes/bitquads/quadCountTreeOfShapesComputer.hpp",
    "include/morphotree/tree/treeOfShapes/kgrid.hpp",
    "include/morphotree/tree/treeOfShapes/order_image.hpp",
    "include/morphotree/tree/treeOfShapes/tos.hpp",
    "src/adjacency/adjacency8c.cpp",
    "src/adjacency/adjacencyuc.cpp",
    "src/attributes/bitquads/quads.cpp",
    "src/attributes/bitquads/quadCountTreeOfShapesComputer.cpp",
    "src/core/box.cpp",
    "resource/quads/dt-max-tree-8c.dat",
    "resource/quads/dt-min-tree-8c.dat",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the fixed SIBGRAPI 2026 protocol and generate the paper table."
    )
    parser.add_argument("--runner", required=True, type=Path, help="compiled bitquad benchmark")
    parser.add_argument(
        "--morphotree-source", required=True, type=Path, help="unmodified MorphoTree snapshot at the documented commit"
    )
    for label, _, _ in RESOLUTIONS:
        parser.add_argument(f"--dataset-{label}", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--count", type=int, default=100)
    parser.add_argument(
        "--allow-dirty",
        action="store_true",
        help="allow development runs from a modified mmcfilters checkout",
    )
    args = parser.parse_args()
    if args.start < 0 or args.count < 1 or args.start + args.count > 100:
        parser.error("require START >= 0, COUNT >= 1, and START + COUNT <= 100")
    return args


def sha256(filename: Path) -> str:
    digest = hashlib.sha256()
    with filename.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def png_dimensions(filename: Path) -> tuple[int, int]:
    with filename.open("rb") as stream:
        header = stream.read(24)
    if len(header) != 24 or header[:8] != b"\x89PNG\r\n\x1a\n" or header[12:16] != b"IHDR":
        raise ValueError(f"not a valid PNG file: {filename}")
    return struct.unpack(">II", header[16:24])


def git_value(repository: Path, *arguments: str) -> Optional[str]:
    result = subprocess.run(
        ["git", "-C", str(repository), *arguments],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        text=True,
    )
    return result.stdout.strip() if result.returncode == 0 else None


def git_succeeds(repository: Path, *arguments: str) -> bool:
    return subprocess.run(
        ["git", "-C", str(repository), *arguments],
        check=False,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    ).returncode == 0


def validate_mmcfilters_checkout(repository: Path, allow_dirty: bool) -> dict[str, object]:
    if git_value(repository, "rev-parse", "--show-toplevel") != str(repository.resolve()):
        raise ValueError(f"cannot resolve the mmcfilters repository at {repository}")
    if not git_succeeds(repository, "merge-base", "--is-ancestor", BASE_RELEASE_TAG, "HEAD"):
        raise ValueError(f"the experiment commit must descend from {BASE_RELEASE_TAG}")
    commit = git_value(repository, "rev-parse", "HEAD")
    if commit is None:
        raise ValueError(f"cannot resolve the mmcfilters commit at {repository}")
    status = git_value(repository, "status", "--porcelain")
    if status is None:
        raise ValueError(f"cannot inspect the mmcfilters checkout at {repository}")
    if status and not allow_dirty:
        raise ValueError(
            "mmcfilters has local changes; commit the protocol first or pass --allow-dirty for a development run"
        )
    return {
        "commit": commit,
        "dirty": bool(status),
        "base_release": BASE_RELEASE_TAG,
        "exact_tag": git_value(repository, "describe", "--tags", "--exact-match", "HEAD"),
    }


def ensure_empty_output(directory: Path) -> None:
    if directory.exists():
        if not directory.is_dir():
            raise ValueError(f"output path is not a directory: {directory}")
        if any(directory.iterdir()):
            raise ValueError(f"refusing to reuse nonempty output directory: {directory}")
    else:
        directory.mkdir(parents=True)


def validate_reference(source: Path) -> None:
    missing = [relative for relative in REQUIRED_REFERENCE_FILES if not (source / relative).is_file()]
    if missing:
        raise ValueError("MorphoTree source is missing required original files: " + ", ".join(missing))
    # A source archive has no .git entry. Avoid asking Git to walk upward into
    # an unrelated parent repository in that case.
    if (source / ".git").exists():
        head = git_value(source, "rev-parse", "HEAD")
        # A copied worktree may retain a stale .git pointer; treat it like a
        # source archive when Git cannot resolve the pointer.
        if head is not None and head != REFERENCE5_COMMIT:
            raise ValueError(f"MorphoTree checkout is at {head}, expected {REFERENCE5_COMMIT}")
        dirty = git_value(source, "status", "--porcelain") if head is not None else None
        if dirty:
            raise ValueError("MorphoTree checkout has local changes; reference [5] must remain unmodified")


def write_dataset_manifest(args: argparse.Namespace, output: Path) -> None:
    with (output / "dataset-manifest.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("resolution", "image_index", "filename", "width", "height", "sha256"))
        for label, expected_width, expected_height in RESOLUTIONS:
            dataset = getattr(args, f"dataset_{label}")
            if not dataset.is_dir():
                raise ValueError(f"dataset directory does not exist: {dataset}")
            for image_index in range(args.start, args.start + args.count):
                filename = f"val_{image_index:03d}.png"
                image = dataset / filename
                if not image.is_file():
                    raise ValueError(f"missing dataset image: {image}")
                width, height = png_dimensions(image)
                if (width, height) != (expected_width, expected_height):
                    raise ValueError(
                        f"unexpected dimensions for {label}/{filename}: {width}x{height}; "
                        f"expected {expected_width}x{expected_height}"
                    )
                writer.writerow((label, image_index, filename, width, height, sha256(image)))


def symbolic_commands(args: argparse.Namespace) -> list[list[str]]:
    commands: list[list[str]] = []
    for label, _, _ in RESOLUTIONS:
        commands.append(
            [
                "<RUNNER>",
                "--image-dir",
                f"<DATASET_{label.upper()}>",
                "--resolution",
                label,
                "--output",
                f"raw-{label}.csv",
                "--dt-max-8c",
                "<MORPHOTREE>/resource/quads/dt-max-tree-8c.dat",
                "--dt-min-8c",
                "<MORPHOTREE>/resource/quads/dt-min-tree-8c.dat",
                "--start",
                str(args.start),
                "--count",
                str(args.count),
            ]
        )
    return commands


def write_experiment_manifest(args: argparse.Namespace, output: Path, mmcfilters_provenance: dict[str, object]) -> None:
    manifest = {
        "schema": "mmcfilters.sibgrapi2026.experiment.v2",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "paper": "Unifying Local Attribute Computation on Component Trees and the Tree of Shapes",
        "protocol": {
            "images": args.count,
            "first_image_index": args.start,
            "untimed_warmups_per_method_and_image": 1,
            "timed_repetitions_per_method_and_image": 3,
            "paired_order": "deterministic alternating order",
            "timed_scope": "attribute API only",
            "resampling": "none",
        },
        "mmcfilters": mmcfilters_provenance,
        "morphotree_reference": {
            "repository": "https://github.com/dennisjosesilva/morphotree",
            "commit": REFERENCE5_COMMIT,
            "paper_references": ["5", "6"],
            "source_policy": "compiled from the unmodified source files",
            "source_files_sha256": {
                relative: sha256(args.morphotree_source / relative) for relative in REQUIRED_REFERENCE_FILES
            },
        },
        "runner_sha256": sha256(args.runner),
        "environment": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "processor": platform.processor(),
            "python": platform.python_version(),
        },
        "commands": symbolic_commands(args),
    }
    (output / "experiment.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def run_resolution(args: argparse.Namespace, output: Path, label: str) -> None:
    reference = args.morphotree_source.resolve()
    command = [
        str(args.runner.resolve()),
        "--image-dir",
        str(getattr(args, f"dataset_{label}").resolve()),
        "--resolution",
        label,
        "--output",
        str((output / f"raw-{label}.csv").resolve()),
        "--dt-max-8c",
        str(reference / "resource/quads/dt-max-tree-8c.dat"),
        "--dt-min-8c",
        str(reference / "resource/quads/dt-min-tree-8c.dat"),
        "--start",
        str(args.start),
        "--count",
        str(args.count),
    ]
    with (output / f"run-{label}.log").open("w", encoding="utf-8") as log:
        log.write("command=" + " ".join(symbolic_commands(args)[[item[0] for item in RESOLUTIONS].index(label)]) + "\n")
        log.flush()
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)


def main() -> int:
    args = parse_args()
    repository = Path(__file__).resolve().parents[2]
    if not args.runner.is_file():
        raise ValueError(f"benchmark runner does not exist: {args.runner}")
    mmcfilters_provenance = validate_mmcfilters_checkout(repository, args.allow_dirty)
    validate_reference(args.morphotree_source)
    ensure_empty_output(args.output_dir)
    output = args.output_dir.resolve()
    write_dataset_manifest(args, output)
    write_experiment_manifest(args, output, mmcfilters_provenance)
    for label, _, _ in RESOLUTIONS:
        print(f"running {label}: {args.count} images, three timed repetitions", flush=True)
        run_resolution(args, output, label)
    analyzer = Path(__file__).with_name("analyze_results.py")
    subprocess.run(
        [sys.executable, str(analyzer), str(output), "--start", str(args.start), "--count", str(args.count)],
        check=True,
    )
    print(f"complete: {output}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
