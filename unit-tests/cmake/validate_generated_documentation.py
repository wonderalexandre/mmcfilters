#!/usr/bin/env python3
"""Validate structural properties of generated public documentation."""

from __future__ import annotations

import argparse
from html.parser import HTMLParser
from pathlib import Path
import re
import sys


ATTRIBUTE_ROW = re.compile(
    r"^\|\s*`(Attribute::[^`]+)`\s*\|"
    r"[^|]*\|[^|]*\|[^|]*\|"
    r"\s*`(Increasing|Non-increasing)`\s*\|"
)


class MarkdownTableParser(HTMLParser):
    """Collect text from tables emitted by Doxygen for Markdown tables."""

    def __init__(self) -> None:
        super().__init__()
        self.tables: list[list[list[str]]] = []
        self._table_depth = 0
        self._current_table: list[list[str]] | None = None
        self._current_row: list[str] | None = None
        self._current_cell: list[str] | None = None

    def handle_starttag(
        self, tag: str, attrs: list[tuple[str, str | None]]
    ) -> None:
        attributes = dict(attrs)
        if tag == "table":
            if self._table_depth:
                self._table_depth += 1
            elif "markdownTable" in (attributes.get("class") or "").split():
                self._table_depth = 1
                self._current_table = []
            return

        if self._table_depth != 1:
            return
        if tag == "tr":
            self._current_row = []
        elif tag in {"td", "th"} and self._current_row is not None:
            self._current_cell = []

    def handle_data(self, data: str) -> None:
        if self._current_cell is not None:
            self._current_cell.append(data)

    def handle_endtag(self, tag: str) -> None:
        if not self._table_depth:
            return

        if self._table_depth == 1:
            if tag in {"td", "th"} and self._current_cell is not None:
                assert self._current_row is not None
                self._current_row.append(" ".join("".join(self._current_cell).split()))
                self._current_cell = None
            elif tag == "tr" and self._current_row is not None:
                assert self._current_table is not None
                self._current_table.append(self._current_row)
                self._current_row = None

        if tag == "table":
            self._table_depth -= 1
            if not self._table_depth and self._current_table is not None:
                self.tables.append(self._current_table)
                self._current_table = None


def attribute_rows(source: Path) -> list[tuple[str, str]]:
    lines = source.read_text(encoding="utf-8").splitlines()
    candidates = [line for line in lines if line.startswith("| `Attribute::")]
    rows = [match.groups() for line in candidates if (match := ATTRIBUTE_ROW.match(line))]
    if len(rows) != len(candidates):
        raise RuntimeError(
            "every attribute catalog row must contain an Attribute class of "
            "`Increasing` or `Non-increasing`"
        )
    return rows


def find_catalog_html(artifact_root: Path) -> Path:
    candidates = sorted(artifact_root.glob("*attribute-catalog*.html"))
    if len(candidates) != 1:
        names = ", ".join(path.name for path in candidates) or "none"
        raise RuntimeError(
            "expected one generated attribute-catalog HTML file, "
            f"found {len(candidates)}: {names}"
        )
    return candidates[0]


def validate_attribute_catalog(root: Path, artifact_root: Path) -> None:
    source = root / "docs" / "attribute-catalog.md"
    expected_rows = attribute_rows(source)
    if not expected_rows:
        raise RuntimeError(f"no attribute rows found in {source}")
    expected_enumerators = [enumerator for enumerator, _ in expected_rows]
    expected_classes = [attribute_class for _, attribute_class in expected_rows]

    parser = MarkdownTableParser()
    html_path = find_catalog_html(artifact_root)
    parser.feed(html_path.read_text(encoding="utf-8"))

    catalog_tables = [
        table
        for table in parser.tables
        if table
        and "C++ enumerator" in table[0]
        and "Stable/Python name" in table[0]
        and "Attribute class" in table[0]
    ]
    if len(catalog_tables) != 1:
        raise RuntimeError(
            "expected one rendered attribute catalog table, "
            f"found {len(catalog_tables)}"
        )

    catalog = catalog_tables[0]
    class_column = catalog[0].index("Attribute class")
    rendered_rows = [row for row in catalog[1:] if row]
    if any(len(row) <= class_column for row in rendered_rows):
        raise RuntimeError("rendered attribute catalog contains an incomplete row")
    rendered_enumerators = [row[0] for row in rendered_rows]
    rendered_classes = [row[class_column] for row in rendered_rows]
    if rendered_enumerators != expected_enumerators:
        missing = [name for name in expected_enumerators if name not in rendered_enumerators]
        unexpected = [name for name in rendered_enumerators if name not in expected_enumerators]
        raise RuntimeError(
            "rendered attribute catalog does not match its Markdown source: "
            f"expected {len(expected_enumerators)} rows, found {len(rendered_enumerators)}; "
            f"missing={missing[:5]}, unexpected={unexpected[:5]}"
        )
    if rendered_classes != expected_classes:
        raise RuntimeError(
            "rendered attribute classes do not match their Markdown source"
        )

    print(
        f"Validated {len(rendered_enumerators)} rendered attribute catalog rows "
        "and inclusion classes."
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--artifact-root", type=Path, required=True)
    args = parser.parse_args()

    try:
        validate_attribute_catalog(args.root.resolve(), args.artifact_root.resolve())
    except RuntimeError as error:
        print(f"Documentation validation failed: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
