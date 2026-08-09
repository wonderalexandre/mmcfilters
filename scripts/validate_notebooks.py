#!/usr/bin/env python3
"""Validate repository notebook structure, language, and import policy."""

from __future__ import annotations

import json
from pathlib import Path
import re
import sys

import nbformat


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK_DIRECTORY = REPOSITORY_ROOT / "notebooks"

EXPECTED_NOTEBOOKS = (
    "Attribute_Filters.ipynb",
    "Comparative_Morphological_Tree_Filtering.ipynb",
    "Higra_Attribute_Interoperability.ipynb",
    "MaxDistExample.ipynb",
    "SimpleExamples.ipynb",
    "ToS_Contour_Example.ipynb",
    "UAO_Examples.ipynb",
)

LOCAL_ONLY_NOTEBOOKS = {
    "Saliency_Maps_Tutorial.ipynb",
    "Shape_Space_Saliency_Xu.ipynb",
}

FORBIDDEN_CONTENT = re.compile(
    r"_notebook_mmcfilters|load_local_mmcfilters|load_mmcfilters|"
    r"importlib\.util|runpy|sys\.path|%pip|!pip|pip\s+install|"
    r"build/python/mmcfilters|/(?:Users|Volumes)/",
    re.IGNORECASE,
)
PORTUGUESE_PROSE = re.compile(
    r"\b(?:objetivo|prepara[cç][aã]o|etapas|resultados|conclus[oõ]es|"
    r"pr[oó]ximos|filtragem|imagem|[aá]rvore|sali[eê]ncia|extin[cç][aã]o|"
    r"contornos?|valora[cç][aã]o|verifica[cç][oõ]es|diferen[cç]a|m[aá]ximo|"
    r"m[ií]nimo|n[oó]s|arestas|área|regi[oõ]es|n[ií]veis|limiar|"
    r"parti[cç][oõ]es|persist[eê]ncia|circularidade|raiz|folha)\b",
    re.IGNORECASE,
)
DIRECT_IMPORT = re.compile(r"(?:^|\n)\s*import\s+mmcfilters(?:\s|$)")
FROM_IMPORT = re.compile(r"(?:^|\n)\s*from\s+mmcfilters(?:\.|\s)")


def output_text(cell: nbformat.NotebookNode) -> str:
    return json.dumps(cell.get("outputs", []), ensure_ascii=False, default=str)


def validate_notebook(path: Path) -> list[str]:
    errors: list[str] = []
    try:
        notebook = nbformat.read(path, as_version=4)
        nbformat.validate(notebook)
    except Exception as error:  # noqa: BLE001 - aggregate all validation failures
        return [f"invalid nbformat: {error}"]

    if not notebook.cells or notebook.cells[0].cell_type != "markdown":
        errors.append("the first cell must be a Markdown title/overview")
    elif not notebook.cells[0].source.lstrip().startswith("# "):
        errors.append("the first Markdown cell must start with a level-one title")

    code_cells = [
        cell
        for cell in notebook.cells
        if cell.cell_type == "code" and cell.source.strip()
    ]
    if not code_cells:
        errors.append("the notebook has no non-empty code cell")
        return errors

    if not DIRECT_IMPORT.search(code_cells[0].source):
        errors.append("the first code cell must contain direct `import mmcfilters`")

    counts = [cell.execution_count for cell in code_cells]
    if any(count is not None for count in counts):
        if any(count is None for count in counts):
            errors.append("stored execution counts must cover every non-empty code cell")
        elif counts != list(range(1, len(code_cells) + 1)):
            errors.append("stored execution counts must be sequential 1..N")

    for index, cell in enumerate(notebook.cells):
        content = cell.source
        if cell.cell_type == "code":
            content += "\n" + output_text(cell)
            if FROM_IMPORT.search(cell.source):
                errors.append(f"cell {index}: use `import mmcfilters`, not a from-import")
            for output in cell.get("outputs", []):
                if output.get("output_type") == "error":
                    errors.append(f"cell {index}: committed error output")
        if FORBIDDEN_CONTENT.search(content):
            errors.append(f"cell {index}: loader, install, import-hook, or absolute-path content")
        if PORTUGUESE_PROSE.search(content):
            errors.append(f"cell {index}: Portuguese reader-facing prose or labels")

    return errors


def main() -> int:
    present = tuple(
        sorted(
            path.name
            for path in NOTEBOOK_DIRECTORY.glob("*.ipynb")
            if path.name not in LOCAL_ONLY_NOTEBOOKS
        )
    )
    expected = tuple(sorted(EXPECTED_NOTEBOOKS))
    failures: list[str] = []

    if present != expected:
        missing = sorted(set(expected) - set(present))
        unexpected = sorted(set(present) - set(expected))
        if missing:
            failures.append(f"missing notebooks: {', '.join(missing)}")
        if unexpected:
            failures.append(f"unexpected notebooks: {', '.join(unexpected)}")

    if (NOTEBOOK_DIRECTORY / "_notebook_mmcfilters.py").exists():
        failures.append("the removed local-extension notebook loader is present")

    for name in EXPECTED_NOTEBOOKS:
        path = NOTEBOOK_DIRECTORY / name
        if not path.exists():
            continue
        for error in validate_notebook(path):
            failures.append(f"{name}: {error}")

    if failures:
        print("Notebook validation failed:", file=sys.stderr)
        for failure in failures:
            print(f"- {failure}", file=sys.stderr)
        return 1

    print(f"Validated {len(EXPECTED_NOTEBOOKS)} notebooks.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
