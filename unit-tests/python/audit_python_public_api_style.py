#!/usr/bin/env python3
"""Prevent new nonconforming names from entering the public Python API."""

from __future__ import annotations

import argparse
import importlib.util
import inspect
import json
import pathlib
import re
import sys
import types


SNAKE_CASE = re.compile(r"^[a-z][a-z0-9]*(?:_[a-z0-9]+)*$")
PASCAL_CASE = re.compile(r"^[A-Z][A-Za-z0-9]*$")
UPPER_SNAKE_CASE = re.compile(r"^[A-Z][A-Z0-9]*(?:_[A-Z0-9]+)*$")


def is_public(name: str) -> bool:
    return bool(name) and not name.startswith("_")


def module_name_conforms(name: str, value: object) -> bool:
    if isinstance(value, type) or isinstance(value, types.ModuleType):
        return PASCAL_CASE.fullmatch(name) is not None or SNAKE_CASE.fullmatch(name) is not None
    return SNAKE_CASE.fullmatch(name) is not None or UPPER_SNAKE_CASE.fullmatch(name) is not None


def member_name_conforms(name: str, value: object) -> bool:
    if isinstance(value, type):
        return PASCAL_CASE.fullmatch(name) is not None
    return SNAKE_CASE.fullmatch(name) is not None or UPPER_SNAKE_CASE.fullmatch(name) is not None


def safe_static_attribute(owner: object, name: str) -> object:
    try:
        return inspect.getattr_static(owner, name)
    except AttributeError:
        return object()


def collect_nonconforming(module: types.ModuleType) -> set[str]:
    failures: set[str] = set()
    exported = getattr(module, "__all__", [name for name in dir(module) if is_public(name)])
    for name in sorted(set(exported)):
        if name == "__version__" or not hasattr(module, name):
            continue
        value = getattr(module, name)
        if not module_name_conforms(name, value):
            failures.add(name)
        if not isinstance(value, type):
            continue
        for member in dir(value):
            if not is_public(member):
                continue
            member_value = safe_static_attribute(value, member)
            if not member_name_conforms(member, member_value):
                failures.add(f"{name}.{member}")
    return failures


def load_build_package(build_dir: pathlib.Path) -> types.ModuleType:
    package_dir = build_dir.resolve() / "python" / "mmcfilters"
    package_init = package_dir / "__init__.py"
    if not package_init.is_file():
        raise RuntimeError(f"package init not found: {package_init}")

    for name in list(sys.modules):
        if name == "mmcfilters" or name.startswith("mmcfilters."):
            sys.modules.pop(name, None)

    def handles_mmcfilters(finder: object) -> bool:
        known_modules: dict[str, object] = {}
        known_modules.update(getattr(finder, "known_source_files", {}))
        known_modules.update(getattr(finder, "known_wheel_files", {}))
        return any(name == "mmcfilters" or name.startswith("mmcfilters.") for name in known_modules)

    sys.meta_path = [finder for finder in sys.meta_path if not handles_mmcfilters(finder)]
    spec = importlib.util.spec_from_file_location("mmcfilters", package_init, submodule_search_locations=[str(package_dir)])
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot create module specification for {package_init}")
    module = importlib.util.module_from_spec(spec)
    sys.modules["mmcfilters"] = module
    spec.loader.exec_module(module)
    return module


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", required=True, type=pathlib.Path)
    parser.add_argument("--baseline", type=pathlib.Path)
    parser.add_argument("--report", action="store_true")
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    mmcfilters = load_build_package(arguments.build_dir)

    observed = collect_nonconforming(mmcfilters)
    if arguments.report:
        print(json.dumps(sorted(observed), indent=2))
        return 0

    allowed: set[str] = set()
    if arguments.baseline is not None:
        baseline = json.loads(arguments.baseline.read_text(encoding="utf-8"))
        if baseline.get("schema_version") != 1:
            print(f"{arguments.baseline}: schema_version must be 1", file=sys.stderr)
            return 1
        allowed = set(baseline.get("allowed_nonconforming", []))
    new_failures = sorted(observed - allowed)
    if new_failures:
        print("Python public API style audit found new nonconforming names:", file=sys.stderr)
        print("\n".join(new_failures), file=sys.stderr)
        return 1
    if arguments.baseline is None:
        print("Python public API style audit passed: all public names conform")
    else:
        print(
            f"Python public API style audit passed: {len(observed)} baseline nonconforming names remain, "
            f"{len(allowed - observed)} have been removed"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
