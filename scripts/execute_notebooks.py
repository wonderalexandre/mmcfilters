#!/usr/bin/env python3
"""Execute maintained notebooks into isolated output copies."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from nbclient import NotebookClient
import nbformat

from validate_notebooks import EXPECTED_NOTEBOOKS, NOTEBOOK_DIRECTORY


SMOKE_NOTEBOOKS = (
    "Attribute_Filters.ipynb",
    "Higra_Attribute_Interoperability.ipynb",
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", choices=("smoke", "full"), default="smoke")
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=NOTEBOOK_DIRECTORY.parent / "build" / "notebook-runs",
    )
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    names = SMOKE_NOTEBOOKS if arguments.profile == "smoke" else EXPECTED_NOTEBOOKS
    arguments.output_dir.mkdir(parents=True, exist_ok=True)

    for name in names:
        source_path = NOTEBOOK_DIRECTORY / name
        output_path = arguments.output_dir / name
        print(f"Executing {name}...", flush=True)
        notebook = nbformat.read(source_path, as_version=4)
        client = NotebookClient(
            notebook,
            timeout=arguments.timeout,
            kernel_name="python3",
            allow_errors=False,
            record_timing=False,
            resources={"metadata": {"path": str(NOTEBOOK_DIRECTORY)}},
        )
        try:
            client.execute()
        except Exception as error:  # noqa: BLE001 - show failing notebook in CI
            print(f"{name} failed: {error}", file=sys.stderr)
            return 1
        nbformat.write(notebook, output_path)

    print(f"Executed {len(names)} notebooks into {arguments.output_dir}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
