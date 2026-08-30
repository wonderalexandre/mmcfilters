#!/usr/bin/env python3
"""Validate structural properties of generated public documentation."""

from __future__ import annotations

import argparse
from html.parser import HTMLParser
from pathlib import Path
import re
import sys


ATTRIBUTE_ROW = re.compile(r"^\|\s*`(Attribute::[^`]+)`\s*\|")


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


def attribute_enumerators(source: Path) -> list[str]:
    return [
        match.group(1)
        for line in source.read_text(encoding="utf-8").splitlines()
        if (match := ATTRIBUTE_ROW.match(line))
    ]


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
    expected = attribute_enumerators(source)
    if not expected:
        raise RuntimeError(f"no attribute rows found in {source}")

    parser = MarkdownTableParser()
    html_path = find_catalog_html(artifact_root)
    parser.feed(html_path.read_text(encoding="utf-8"))

    catalog_tables = [
        table
        for table in parser.tables
        if table
        and "C++ enumerator" in table[0]
        and "Stable/Python name" in table[0]
    ]
    if len(catalog_tables) != 1:
        raise RuntimeError(
            "expected one rendered attribute catalog table, "
            f"found {len(catalog_tables)}"
        )

    rendered = [row[0] for row in catalog_tables[0][1:] if row]
    if rendered != expected:
        missing = [name for name in expected if name not in rendered]
        unexpected = [name for name in rendered if name not in expected]
        raise RuntimeError(
            "rendered attribute catalog does not match its Markdown source: "
            f"expected {len(expected)} rows, found {len(rendered)}; "
            f"missing={missing[:5]}, unexpected={unexpected[:5]}"
        )

    print(f"Validated {len(rendered)} rendered attribute catalog rows.")


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
