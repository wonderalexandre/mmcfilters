#!/usr/bin/env python3
"""Validate refactoring metadata and progressively reject legacy vocabulary."""

from __future__ import annotations

import argparse
import bisect
import json
import pathlib
import re
import sys
import tempfile


TEXT_SUFFIXES = {
    ".c",
    ".cc",
    ".cmake",
    ".cpp",
    ".cxx",
    ".h",
    ".hpp",
    ".html",
    ".in",
    ".ini",
    ".ipynb",
    ".js",
    ".json",
    ".md",
    ".py",
    ".rst",
    ".toml",
    ".txt",
    ".xml",
    ".yaml",
    ".yml",
}
ALLOWED_STATUSES = {"not_started", "in_progress", "implemented", "verified", "not_applicable"}
MATRIX_GROUPS = {
    "G": 5,
    "T": 14,
    "V": 14,
    "L": 12,
    "B": 14,
    "R": 18,
    "F": 12,
    "P": 9,
    "C": 6,
}
CANONICAL_UPPERCASE_CPP_ENUM_MEMBERS = {"Q1", "Q2", "QD", "Q3", "Q4"}
PYTHON_SNAKE_CASE = re.compile(r"^[a-z][a-z0-9]*(?:_[a-z0-9]+)*$")
PROJECT_TEXT_PATHS = (
    ".github",
    "CMakeLists.txt",
    "MANIFEST.in",
    "README.md",
    "benchmarks",
    "cmake",
    "docs",
    "examples",
    "mmcfilters",
    "notebooks",
    "pybind11.cmake",
    "pybinds",
    "pyproject.toml",
    "python",
    "scripts",
    "unit-tests",
)
COORDINATE_VOCABULARY_EXCLUSIONS = {
    pathlib.PurePosixPath("benchmarks/results"),
    pathlib.PurePosixPath("docs/refactoring"),
    pathlib.PurePosixPath("notebooks/.ipynb_checkpoints"),
    pathlib.PurePosixPath("unit-tests/cmake/audit_api_refactoring.py"),
}


def expected_matrix_ids() -> set[str]:
    return {f"{prefix}-{index:02d}" for prefix, count in MATRIX_GROUPS.items() for index in range(1, count + 1)}


def load_json(path: pathlib.Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"cannot read {path}: {error}") from error


def validate_migration_status(path: pathlib.Path) -> list[str]:
    failures: list[str] = []
    document = load_json(path)
    if document.get("schema_version") != 1:
        failures.append(f"{path}: schema_version must be 1")

    rows = document.get("rows")
    if not isinstance(rows, dict):
        return failures + [f"{path}: rows must be an object keyed by migration ID"]

    expected = expected_matrix_ids()
    actual = set(rows)
    for missing in sorted(expected - actual):
        failures.append(f"{path}: missing migration row {missing}")
    for unexpected in sorted(actual - expected):
        failures.append(f"{path}: unexpected migration row {unexpected}")

    for row_id, row in sorted(rows.items()):
        if not isinstance(row, dict):
            failures.append(f"{path}: {row_id} must be an object")
            continue
        status = row.get("status")
        if status not in ALLOWED_STATUSES:
            failures.append(f"{path}: {row_id} has invalid status {status!r}")
        change_set = row.get("change_set")
        if not isinstance(change_set, str) or re.fullmatch(r"RF-(?:0\d|1\d|2[01])", change_set) is None:
            failures.append(f"{path}: {row_id} has invalid change_set {change_set!r}")
        core_required = row.get("core_required")
        if not isinstance(core_required, bool):
            failures.append(f"{path}: {row_id} core_required must be Boolean")
        elif row_id.startswith("C-") == core_required:
            failures.append(f"{path}: {row_id} core_required must be false only for CFP rows")
        if status == "not_applicable" and not row.get("justification"):
            failures.append(f"{path}: {row_id} marked not_applicable requires justification")

    if len(rows) != 104:
        failures.append(f"{path}: expected 104 migration rows, found {len(rows)}")

    program_state = document.get("program_state")
    if not isinstance(program_state, dict):
        failures.append(f"{path}: program_state must be an object")
        return failures
    completed = program_state.get("completed_change_sets")
    next_change_set = program_state.get("next_change_set")
    if not isinstance(completed, list) or any(not isinstance(item, str) for item in completed):
        failures.append(f"{path}: completed_change_sets must be an array of strings")
        completed = []
    elif len(completed) != len(set(completed)):
        failures.append(f"{path}: completed_change_sets must not contain duplicates")
    if not isinstance(next_change_set, str) or re.fullmatch(r"RF-(?:0\d|1\d|2[01])", next_change_set) is None:
        failures.append(f"{path}: next_change_set must identify RF-00 through RF-21")
    elif next_change_set in completed:
        failures.append(f"{path}: next_change_set {next_change_set} is already completed")

    if "RF-17" in completed:
        non_verified_core = sorted(
            row_id for row_id, row in rows.items()
            if isinstance(row, dict) and row.get("core_required") is True and row.get("status") != "verified"
        )
        if non_verified_core:
            failures.append(
                f"{path}: completed RF-17 requires every core migration row to be verified; remaining {non_verified_core}"
            )
    return failures


def cpp_without_comments_and_literals(text: str) -> str:
    """Replace comments and quoted literals while retaining source line numbers."""
    output: list[str] = []
    index = 0
    while index < len(text):
        if text.startswith("//", index):
            end = text.find("\n", index)
            if end < 0:
                output.extend(" " for _ in text[index:])
                break
            output.extend(" " for _ in text[index:end])
            index = end
            continue
        if text.startswith("/*", index):
            end = text.find("*/", index + 2)
            if end < 0:
                end = len(text) - 2
            end += 2
            output.extend("\n" if character == "\n" else " " for character in text[index:end])
            index = end
            continue
        if text[index] in {'"', "'"}:
            quote = text[index]
            end = index + 1
            while end < len(text):
                if text[end] == "\\":
                    end += 2
                    continue
                end += 1
                if text[end - 1] == quote:
                    break
            output.extend("\n" if character == "\n" else " " for character in text[index:end])
            index = end
            continue
        output.append(text[index])
        index += 1
    return "".join(output)


def is_cpp_pascal_case_enum_member(member: str) -> bool:
    if member in CANONICAL_UPPERCASE_CPP_ENUM_MEMBERS:
        return True
    if re.fullmatch(r"[A-Z][A-Za-z0-9]*", member) is None:
        return False
    return len(member) == 1 or re.fullmatch(r"[A-Z][A-Z0-9]*", member) is None


def validate_cpp_enum_members(root: pathlib.Path) -> list[str]:
    """Require PascalCase members in project headers, with canonical bitquad-code exceptions."""
    failures: list[str] = []
    header_root = root / "mmcfilters"
    for path in sorted(header_root.rglob("*.hpp")):
        cleaned = cpp_without_comments_and_literals(path.read_text(encoding="utf-8"))
        for enum_match in re.finditer(r"\benum(?:\s+class)?\s+([A-Za-z_]\w*)[^;{]*\{([^{}]*)\}", cleaned, re.DOTALL):
            enum_name, body = enum_match.groups()
            body_start = enum_match.start(2)
            for member_match in re.finditer(r"(?:^|,)\s*([A-Za-z_]\w*)", body):
                member = member_match.group(1)
                if is_cpp_pascal_case_enum_member(member):
                    continue
                line = cleaned.count("\n", 0, body_start + member_match.start(1)) + 1
                relative = path.relative_to(root)
                failures.append(
                    f"{relative}:{line}: C++ enum member {enum_name}::{member} must use PascalCase"
                )
    return failures


def validate_identifier_domains(root: pathlib.Path) -> list[str]:
    """Reject declarations and sentinels that confuse pixel and node domains."""
    failures: list[str] = []
    source_roots = (
        root / "mmcfilters",
        root / "pybinds",
        root / "unit-tests",
        root / "benchmarks",
        root / "examples",
    )
    pixel_names = (
        "anchorPixel",
        "boundaryPixel",
        "flatZoneRepresentative",
        "infinityPixel",
        "pixel",
        "pixelId",
        "properPart",
        "samplePixel",
        "spatialMinimum",
    )
    pixel_name_expression = "|".join(pixel_names)
    patterns = (
        (
            re.compile(rf"\b(?:const\s+)?(?:int|NodeId)\s+(?:const\s+)?(?:{pixel_name_expression})\b"),
            "pixel-domain identifier must use PixelId, not int or NodeId",
        ),
        (
            re.compile(
                r"\bstd::(?:vector|span)<\s*(?:const\s+)?NodeId\s*>\s*"
                r"[A-Za-z0-9_]*(?:BoundaryPixels|FlatZoneMergeRepresentatives|ProperParts|SupportPixels|"
                r"boundaryPixels|flatZoneMergeRepresentatives|properParts|supportPixels|unitProperParts|unitPixels)\b"
            ),
            "pixel-domain collection must contain PixelId, not NodeId",
        ),
        (
            re.compile(r"\bPixelId\s+(?:const\s+)?properPart\b"),
            "a pixel-domain element must be named pixel, not properPart",
        ),
        (
            re.compile(
                r"(?:\bPixelId\s+[A-Za-z_]\w*|\bstd::vector<\s*PixelId\s*>\s+[A-Za-z_]\w*)"
                r"[^;\n]*\bInvalidNode\b"
            ),
            "PixelId storage must use InvalidPixel, not InvalidNode",
        ),
        (
            re.compile(
                rf"\b(?:{pixel_name_expression})\s*(?:==|!=)\s*InvalidNode\b|"
                rf"\bInvalidNode\s*(?:==|!=)\s*(?:{pixel_name_expression})\b"
            ),
            "pixel-domain comparison must use InvalidPixel, not InvalidNode",
        ),
        (
            re.compile(
                r"\b(?:const\s+)?(?:NodeId|PixelId)\s+(?:const\s+)?"
                r"(?:count|num|number|size)[A-Za-z0-9_]*\b"
            ),
            "identifier types represent domain elements, not counts",
        ),
    )

    for source_root in source_roots:
        if not source_root.is_dir():
            continue
        for path in sorted(source_root.rglob("*")):
            if path.suffix not in {".cpp", ".hpp"}:
                continue
            cleaned = cpp_without_comments_and_literals(path.read_text(encoding="utf-8"))
            relative = path.relative_to(root)
            for pattern, message in patterns:
                for match in pattern.finditer(cleaned):
                    line = cleaned.count("\n", 0, match.start()) + 1
                    failures.append(f"{relative}:{line}: {message}: {match.group(0).strip()}")
    return failures


def uses_column_abbreviation(identifier: str) -> bool:
    """Recognize column abbreviations without confusing color, collection, or protocol."""
    if re.search(r"(?:^|_)(?:col|cols)(?:_|$)", identifier, re.IGNORECASE):
        return True
    if re.search(r"(?:^|[a-z0-9])Cols?(?=[A-Z0-9]|$)", identifier):
        return True
    if re.search(r"^Cols?(?=[A-Z0-9]|$)", identifier):
        return True
    return re.fullmatch(r"ncols", identifier, re.IGNORECASE) is not None


def validate_coordinate_vocabulary(root: pathlib.Path) -> list[str]:
    """Require the complete row/column words throughout project source and documentation."""
    failures: list[str] = []
    for configured in PROJECT_TEXT_PATHS:
        base = root / configured
        paths = [base] if base.is_file() else base.rglob("*") if base.is_dir() else []
        for path in sorted(paths):
            if not path.is_file():
                continue
            relative = pathlib.PurePosixPath(path.relative_to(root).as_posix())
            if any(
                relative == exclusion or exclusion in relative.parents
                for exclusion in COORDINATE_VOCABULARY_EXCLUSIONS
            ):
                continue
            if path.name != "CMakeLists.txt" and path.suffix not in TEXT_SUFFIXES:
                continue
            try:
                text = audit_text(path)
            except (UnicodeDecodeError, json.JSONDecodeError):
                continue
            for match in re.finditer(r"\b[A-Za-z_][A-Za-z0-9_]*\b", text):
                identifier = match.group(0)
                if not uses_column_abbreviation(identifier):
                    continue
                line = text.count("\n", 0, match.start()) + 1
                failures.append(
                    f"{relative}:{line}: column vocabulary must use column/columns, not {identifier!r}"
                )
    return failures


def validate_python_binding_names(root: pathlib.Path) -> list[str]:
    """Require snake_case for bound functions, methods, properties, and named arguments."""
    failures: list[str] = []
    binding_root = root / "pybinds"
    public_name_pattern = re.compile(
        r'(?:\.(?:def|def_static|def_property|def_property_readonly|def_readwrite|def_readonly)|\bm\.def)\s*\(\s*"([^"]+)"'
    )
    named_argument_patterns = (
        re.compile(r'"([^"]+)"_a'),
        re.compile(r'py::arg\(\s*"([^"]+)"\s*\)'),
    )

    for path in sorted(binding_root.rglob("*")):
        if path.suffix not in {".cpp", ".hpp"}:
            continue
        text = path.read_text(encoding="utf-8")
        relative = path.relative_to(root)
        for match in public_name_pattern.finditer(text):
            name = match.group(1)
            if name.startswith("__") and name.endswith("__"):
                continue
            if PYTHON_SNAKE_CASE.fullmatch(name) is None:
                line = text.count("\n", 0, match.start(1)) + 1
                failures.append(f"{relative}:{line}: Python binding name {name!r} must use snake_case")
        for pattern in named_argument_patterns:
            for match in pattern.finditer(text):
                name = match.group(1)
                if PYTHON_SNAKE_CASE.fullmatch(name) is None:
                    line = text.count("\n", 0, match.start(1)) + 1
                    failures.append(f"{relative}:{line}: Python named argument {name!r} must use snake_case")
    return failures


def validate_doxygen_semantics(root: pathlib.Path) -> list[str]:
    """Reject tautological Doxygen templates that do not document semantic roles."""
    failures: list[str] = []
    patterns = (
        (re.compile(r"@brief\s+Stores the\b"), "field brief must describe meaning, not merely storage"),
        (
            re.compile(r"@param\s+([A-Za-z_]\w*)\s+[^.\n]+ represented by `\1`\."),
            "parameter documentation must not repeat the parameter name as its meaning",
        ),
        (re.compile(r"@param\s+[^\n]+ used by the operation\."), "parameter documentation must state the parameter's role"),
        (re.compile(r"@return\s+Reference to the resulting object\."), "return documentation must identify the returned state"),
    )
    for source_root in (root / "mmcfilters", root / "pybinds"):
        if not source_root.is_dir():
            continue
        for path in sorted(source_root.rglob("*")):
            if path.suffix not in {".cpp", ".hpp"}:
                continue
            text = path.read_text(encoding="utf-8")
            for pattern, message in patterns:
                for match in pattern.finditer(text):
                    line = text.count("\n", 0, match.start()) + 1
                    failures.append(f"{path.relative_to(root)}:{line}: {message}")
    return failures


def pattern_for(entry: dict) -> re.Pattern[str]:
    term = entry["term"]
    mode = entry.get("mode", "identifier")
    escaped = re.escape(term)
    if mode == "identifier":
        expression = rf"(?<![A-Za-z0-9_]){escaped}(?![A-Za-z0-9_])"
    elif mode == "prefix":
        expression = rf"(?<![A-Za-z0-9_]){escaped}[A-Za-z0-9_]*"
    elif mode == "phrase":
        expression = escaped.replace(r"\ ", r"\s+")
    elif mode == "substring":
        expression = escaped
    else:
        raise ValueError(f"unsupported vocabulary match mode {mode!r} for {term!r}")
    flags = re.IGNORECASE if entry.get("ignore_case", False) else 0
    return re.compile(expression, flags)


def candidate_paths(root: pathlib.Path, config: dict):
    excluded = {pathlib.PurePosixPath(item) for item in config.get("exclude", [])}
    for configured in config.get("scan", []):
        base = root / configured
        paths = [base] if base.is_file() else base.rglob("*") if base.is_dir() else []
        for path in paths:
            if not path.is_file():
                continue
            relative = pathlib.PurePosixPath(path.relative_to(root).as_posix())
            if any(relative == item or item in relative.parents for item in excluded):
                continue
            if path.name == "CMakeLists.txt" or path.suffix in TEXT_SUFFIXES:
                yield path


def audit_text(path: pathlib.Path) -> str:
    text = path.read_text(encoding="utf-8")
    if path.suffix != ".ipynb":
        return text
    notebook = json.loads(text)
    sources: list[str] = []
    for cell in notebook.get("cells", []):
        source = cell.get("source", [])
        sources.append("".join(source) if isinstance(source, list) else str(source))
    return "\n".join(sources)


def vocabulary_observations(root: pathlib.Path, config: dict) -> dict[str, list[str]]:
    entries = config.get("entries", [])
    compiled = [(entry, pattern_for(entry)) for entry in entries]
    observations = {entry["term"]: [] for entry in entries}
    for path in sorted(set(candidate_paths(root, config))):
        try:
            text = audit_text(path)
        except (UnicodeDecodeError, json.JSONDecodeError):
            continue
        relative = path.relative_to(root)
        line_starts = [0]
        line_starts.extend(match.end() for match in re.finditer("\n", text))
        lowercase_text: str | None = None
        for entry, pattern in compiled:
            search_hint = entry["term"].split()[0]
            if entry.get("ignore_case", False):
                if lowercase_text is None:
                    lowercase_text = text.lower()
                if search_hint.lower() not in lowercase_text:
                    continue
            elif search_hint not in text:
                continue
            for match in pattern.finditer(text):
                line_number = bisect.bisect_right(line_starts, match.start())
                column = match.start() - line_starts[line_number - 1] + 1
                observations[entry["term"]].append(f"{relative}:{line_number}:{column}:{match.group(0)}")
    return observations


def validate_vocabulary(root: pathlib.Path, path: pathlib.Path, report_only: bool) -> tuple[list[str], dict[str, int]]:
    failures: list[str] = []
    config = load_json(path)
    if config.get("schema_version") != 1:
        failures.append(f"{path}: schema_version must be 1")
    complete_ledger = config.get("complete_ledger", False)
    if not isinstance(complete_ledger, bool):
        failures.append(f"{path}: complete_ledger must be Boolean")
    entries = config.get("entries")
    if not isinstance(entries, list):
        return failures + [f"{path}: entries must be an array"], {}

    terms = [entry.get("term") for entry in entries if isinstance(entry, dict)]
    if len(terms) != len(set(terms)) or any(not isinstance(term, str) or not term for term in terms):
        failures.append(f"{path}: vocabulary terms must be unique nonempty strings")

    observations = vocabulary_observations(root, config)
    counts = {term: len(matches) for term, matches in observations.items()}
    if report_only:
        return failures, counts

    known_ids = expected_matrix_ids()
    for entry in entries:
        term = entry["term"]
        state = entry.get("state")
        if state not in {"pending", "active"}:
            failures.append(f"{path}: {term!r} state must be pending or active")
            continue
        matrix_rows = entry.get("matrix_rows", [])
        if not matrix_rows or any(row not in known_ids for row in matrix_rows):
            failures.append(f"{path}: {term!r} must reference valid matrix_rows")
        count = counts[term]
        if state == "active" and count:
            failures.append(f"active removed term {term!r} has {count} public matches: {observations[term][:8]}")
        if state == "pending":
            baseline_maximum = entry.get("baseline_max_matches")
            if not isinstance(baseline_maximum, int) or baseline_maximum < 0:
                failures.append(f"{path}: pending term {term!r} needs nonnegative baseline_max_matches")
            elif count > baseline_maximum:
                failures.append(
                    f"pending legacy term {term!r} increased from at most {baseline_maximum} to {count}: {observations[term][:8]}"
                )
    if complete_ledger:
        pending_terms = sorted(entry["term"] for entry in entries if entry.get("state") != "active")
        if pending_terms:
            failures.append(f"{path}: complete_ledger requires every removal-ledger term to be active; pending {pending_terms}")
    return failures, counts


def artifact_candidate_paths(root: pathlib.Path):
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        if path.name == "CMakeLists.txt" or path.suffix.lower() in TEXT_SUFFIXES:
            yield path


def validate_artifact_vocabulary(artifact_roots: list[pathlib.Path], config: dict) -> list[str]:
    """Reject active legacy vocabulary in installed or generated publication artifacts."""
    failures: list[str] = []
    active_entries = [entry for entry in config.get("entries", []) if entry.get("state") == "active"]
    compiled = [(entry, pattern_for(entry)) for entry in active_entries]
    for artifact_root in artifact_roots:
        resolved = artifact_root.resolve()
        if not resolved.is_dir():
            failures.append(f"artifact root does not exist or is not a directory: {resolved}")
            continue
        for path in sorted(artifact_candidate_paths(resolved)):
            try:
                text = audit_text(path)
            except (UnicodeDecodeError, json.JSONDecodeError):
                continue
            lowercase_text: str | None = None
            for entry, pattern in compiled:
                search_hint = entry["term"].split()[0]
                if entry.get("ignore_case", False):
                    if lowercase_text is None:
                        lowercase_text = text.lower()
                    if search_hint.lower() not in lowercase_text:
                        continue
                elif search_hint not in text:
                    continue
                match = pattern.search(text)
                if match is None:
                    continue
                line = text.count("\n", 0, match.start()) + 1
                failures.append(
                    f"{resolved.name}/{path.relative_to(resolved)}:{line}: published artifact contains removed term "
                    f"{entry['term']!r} as {match.group(0)!r}"
                )
    return failures


def validate_contextual_vocabulary(root: pathlib.Path, path: pathlib.Path) -> tuple[list[str], dict]:
    """Report generic terms by context while rejecting only scientifically invalid roles."""
    failures: list[str] = []
    config = load_json(path)
    if config.get("schema_version") != 1:
        failures.append(f"{path}: schema_version must be 1")
    if not isinstance(config.get("reviewed_on"), str):
        failures.append(f"{path}: reviewed_on must be a date string")
    terms = config.get("terms")
    if not isinstance(terms, list):
        return failures + [f"{path}: terms must be an array"], {}

    report = {
        "schema_version": 1,
        "reviewed_on": config.get("reviewed_on"),
        "terms": {},
    }
    paths = sorted(set(candidate_paths(root, config)))
    seen_terms: set[str] = set()
    for entry in terms:
        if not isinstance(entry, dict):
            failures.append(f"{path}: contextual term entries must be objects")
            continue
        term = entry.get("term")
        valid_roles = entry.get("valid_roles")
        forbidden_patterns = entry.get("forbidden_patterns")
        if not isinstance(term, str) or not term or term in seen_terms:
            failures.append(f"{path}: contextual terms must be unique nonempty strings")
            continue
        seen_terms.add(term)
        if not isinstance(valid_roles, list) or not valid_roles or any(not isinstance(role, str) or not role for role in valid_roles):
            failures.append(f"{path}: contextual term {term!r} requires nonempty valid_roles")
            valid_roles = []
        if not isinstance(forbidden_patterns, list) or any(not isinstance(pattern, str) or not pattern for pattern in forbidden_patterns):
            failures.append(f"{path}: contextual term {term!r} forbidden_patterns must be an array of nonempty regex strings")
            forbidden_patterns = []

        term_pattern = re.compile(re.escape(term), re.IGNORECASE)
        compiled_forbidden: list[re.Pattern[str]] = []
        for expression in forbidden_patterns:
            try:
                compiled_forbidden.append(re.compile(expression, re.IGNORECASE))
            except re.error as error:
                failures.append(f"{path}: invalid forbidden pattern {expression!r} for {term!r}: {error}")

        occurrence_count = 0
        representative_contexts: list[str] = []
        forbidden_matches: list[str] = []
        for source_path in paths:
            try:
                text = audit_text(source_path)
            except (UnicodeDecodeError, json.JSONDecodeError):
                continue
            relative = source_path.relative_to(root)
            lines = text.splitlines()
            for line_number, line_text in enumerate(lines, 1):
                matches = list(term_pattern.finditer(line_text))
                occurrence_count += len(matches)
                if matches and len(representative_contexts) < 25:
                    representative_contexts.append(f"{relative}:{line_number}:{line_text.strip()[:240]}")
                for forbidden in compiled_forbidden:
                    match = forbidden.search(line_text)
                    if match is None:
                        continue
                    location = f"{relative}:{line_number}:{match.group(0)}"
                    forbidden_matches.append(location)
                    failures.append(f"contextual term {term!r} has invalid role at {location}")

        report["terms"][term] = {
            "occurrences": occurrence_count,
            "valid_roles": valid_roles,
            "representative_contexts": representative_contexts,
            "forbidden_matches": forbidden_matches,
        }
    return failures, report


def run_self_test() -> int:
    with tempfile.TemporaryDirectory() as directory:
        root = pathlib.Path(directory)
        (root / "public").mkdir()
        (root / "history").mkdir()
        (root / "mmcfilters").mkdir()
        (root / "public" / "api.hpp").write_text("OldName pendingName pendingName\n", encoding="utf-8")
        (root / "history" / "migration.md").write_text("OldName\n", encoding="utf-8")
        config = {
            "schema_version": 1,
            "scan": ["public", "history"],
            "exclude": ["history"],
            "entries": [
                {"term": "OldName", "state": "active", "matrix_rows": ["T-01"], "baseline_max_matches": 0},
                {"term": "pendingName", "state": "pending", "matrix_rows": ["T-01"], "baseline_max_matches": 2},
            ],
        }
        config_path = root / "config.json"
        config_path.write_text(json.dumps(config), encoding="utf-8")
        failures, counts = validate_vocabulary(root, config_path, report_only=False)
        if counts != {"OldName": 1, "pendingName": 2}:
            print(f"API refactoring audit self-test observed unexpected counts: {counts}", file=sys.stderr)
            return 1
        if len(failures) != 1 or "active removed term" not in failures[0]:
            print(f"API refactoring audit self-test failed to reject active term: {failures}", file=sys.stderr)
            return 1
        (root / "public" / "api.hpp").write_text("pendingName pendingName pendingName\n", encoding="utf-8")
        failures, _ = validate_vocabulary(root, config_path, report_only=False)
        if len(failures) != 1 or "increased" not in failures[0]:
            print(f"API refactoring audit self-test failed to enforce pending ratchet: {failures}", file=sys.stderr)
            return 1
        artifact_root = root / "artifact"
        artifact_root.mkdir()
        artifact_file = artifact_root / "published.hpp"
        artifact_file.write_text("OldName\n", encoding="utf-8")
        artifact_failures = validate_artifact_vocabulary([artifact_root], config)
        if len(artifact_failures) != 1 or "published artifact contains removed term" not in artifact_failures[0]:
            print(f"API refactoring audit self-test failed to reject a published legacy symbol: {artifact_failures}", file=sys.stderr)
            return 1
        artifact_file.write_text("CurrentName\n", encoding="utf-8")
        artifact_failures = validate_artifact_vocabulary([artifact_root], config)
        if artifact_failures:
            print(f"API refactoring audit self-test rejected a clean artifact: {artifact_failures}", file=sys.stderr)
            return 1

        contextual_config = {
            "schema_version": 1,
            "reviewed_on": "2026-08-13",
            "scan": ["public"],
            "exclude": [],
            "terms": [{
                "term": "owner",
                "valid_roles": ["resource ownership and lifetime"],
                "forbidden_patterns": [r"proper[ _-]*part[ _-]*owner"],
            }],
        }
        contextual_path = root / "contextual.json"
        contextual_path.write_text(json.dumps(contextual_config), encoding="utf-8")
        (root / "public" / "api.hpp").write_text("resource owner remains alive\n", encoding="utf-8")
        contextual_failures, contextual_report = validate_contextual_vocabulary(root, contextual_path)
        if contextual_failures or contextual_report["terms"]["owner"]["occurrences"] != 1:
            print(f"API refactoring audit self-test rejected a valid contextual use: {contextual_failures}", file=sys.stderr)
            return 1
        (root / "public" / "api.hpp").write_text("proper part owner\n", encoding="utf-8")
        contextual_failures, _ = validate_contextual_vocabulary(root, contextual_path)
        if len(contextual_failures) != 1 or "invalid role" not in contextual_failures[0]:
            print(f"API refactoring audit self-test failed to reject an invalid contextual use: {contextual_failures}", file=sys.stderr)
            return 1

        identifier_path = root / "mmcfilters" / "semantic_ids.hpp"
        identifier_path.write_text(
            "NodeId infinityPixel = InvalidNode;\n"
            "std::vector<NodeId> properParts;\n"
            "PixelId properPart = InvalidPixel;\n"
            "PixelId anchorPixel = InvalidNode;\n"
            "if (samplePixel == InvalidNode) {}\n"
            "NodeId numPixels = 4;\n"
            "int boundaryPixel = 0;\n",
            encoding="utf-8",
        )
        identifier_failures = validate_identifier_domains(root)
        if len(identifier_failures) != 7:
            print(f"API refactoring audit self-test missed identifier-domain errors: {identifier_failures}", file=sys.stderr)
            return 1
        identifier_path.write_text(
            "PixelId infinityPixel = InvalidPixel;\n"
            "std::vector<PixelId> properParts;\n"
            "PixelId pixel = InvalidPixel;\n"
            "PixelId anchorPixel = InvalidPixel;\n"
            "if (samplePixel == InvalidPixel) {}\n"
            "int numPixels = 4;\n"
            "PixelId boundaryPixel = 0;\n",
            encoding="utf-8",
        )
        identifier_failures = validate_identifier_domains(root)
        if identifier_failures:
            print(f"API refactoring audit self-test rejected valid identifier domains: {identifier_failures}", file=sys.stderr)
            return 1

        doxygen_path = root / "mmcfilters" / "doxygen.hpp"
        doxygen_path.write_text(
            "/** @brief Stores the node. */\n"
            "/** @param nodeId Identifier represented by `nodeId`. */\n"
            "/** @return Reference to the resulting object. */\n",
            encoding="utf-8",
        )
        doxygen_failures = validate_doxygen_semantics(root)
        if len(doxygen_failures) != 3:
            print(f"API refactoring audit self-test missed Doxygen tautologies: {doxygen_failures}", file=sys.stderr)
            return 1
        doxygen_path.write_text(
            "/** @brief Dense node identifier held by this record. */\n"
            "/** @param nodeId Dense internal node identifier. */\n"
            "/** @return Mutable reference to the updated object. */\n",
            encoding="utf-8",
        )
        doxygen_failures = validate_doxygen_semantics(root)
        if doxygen_failures:
            print(f"API refactoring audit self-test rejected semantic Doxygen: {doxygen_failures}", file=sys.stderr)
            return 1

        coordinate_path = root / "mmcfilters" / "coordinates.hpp"
        coordinate_path.write_text(
            "int numCols; int box_col_min; int col; int ncols;\n",
            encoding="utf-8",
        )
        coordinate_failures = validate_coordinate_vocabulary(root)
        if len(coordinate_failures) != 4:
            print(f"API refactoring audit self-test missed column abbreviations: {coordinate_failures}", file=sys.stderr)
            return 1
        coordinate_path.write_text(
            "int numColumns; int box_column_min; int column; int color; class Collection {}; const char* protocol;\n",
            encoding="utf-8",
        )
        coordinate_failures = validate_coordinate_vocabulary(root)
        if coordinate_failures:
            print(f"API refactoring audit self-test rejected valid column vocabulary: {coordinate_failures}", file=sys.stderr)
            return 1
    print("API refactoring audit self-test passed")
    return 0


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=pathlib.Path)
    parser.add_argument("--status", type=pathlib.Path)
    parser.add_argument("--vocabulary", type=pathlib.Path)
    parser.add_argument("--contextual", type=pathlib.Path)
    parser.add_argument("--context-report", type=pathlib.Path)
    parser.add_argument("--artifact-root", action="append", default=[], type=pathlib.Path)
    parser.add_argument("--report", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    if arguments.self_test:
        return run_self_test()
    if arguments.root is None:
        print("--root is required outside --self-test", file=sys.stderr)
        return 2
    if (arguments.status is None) != (arguments.vocabulary is None):
        print("--status and --vocabulary must be provided together", file=sys.stderr)
        return 2

    root = arguments.root.resolve()
    failures: list[str] = []
    if arguments.status is not None:
        failures.extend(validate_migration_status(arguments.status.resolve()))
    failures.extend(validate_cpp_enum_members(root))
    failures.extend(validate_identifier_domains(root))
    failures.extend(validate_coordinate_vocabulary(root))
    failures.extend(validate_python_binding_names(root))
    failures.extend(validate_doxygen_semantics(root))
    counts: dict[str, int] = {}
    if arguments.vocabulary is not None:
        vocabulary_failures, counts = validate_vocabulary(
            root, arguments.vocabulary.resolve(), arguments.report
        )
        failures.extend(vocabulary_failures)
    contextual_report: dict | None = None
    if arguments.contextual is not None:
        contextual_failures, contextual_report = validate_contextual_vocabulary(root, arguments.contextual.resolve())
        failures.extend(contextual_failures)
    elif arguments.context_report is not None:
        failures.append("--context-report requires --contextual")

    if arguments.artifact_root and arguments.vocabulary is None:
        failures.append("--artifact-root requires --vocabulary")
    elif arguments.artifact_root:
        vocabulary_config = load_json(arguments.vocabulary.resolve())
        failures.extend(validate_artifact_vocabulary(arguments.artifact_root, vocabulary_config))

    if arguments.context_report is not None and contextual_report is not None:
        report_path = arguments.context_report.resolve()
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(contextual_report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if arguments.report:
        print(json.dumps(counts, indent=2, sort_keys=True))
    if failures:
        print("API refactoring audit failed:", file=sys.stderr)
        print("\n".join(failures), file=sys.stderr)
        return 1
    metadata_summary = (
        f"104 migration rows, {len(counts)} removed vocabulary terms"
        if arguments.status is not None
        else "no internal refactoring metadata"
    )
    print(
        f"API refactoring audit passed: {metadata_summary}, PascalCase C++ enums, "
        f"separate pixel/node identifier domains, semantic Doxygen roles, complete column vocabulary, "
        f"snake_case Python bindings, and "
        f"{0 if contextual_report is None else len(contextual_report['terms'])} contextual terms"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
