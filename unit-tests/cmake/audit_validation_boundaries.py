#!/usr/bin/env python3
"""Reject defensive operations from namespaces named ``detail::kernel``."""

from __future__ import annotations

import pathlib
import re
import sys


FORBIDDEN = (
    "MMCFILTERS_CONTRACT_",
    "throw ",
    ".getChildren(",
    ".getProperParts(",
    ".getNodeParent(",
    ".getProperPartOwner(",
    ".getNodeSubtree(",
    ".getAliveNodeIds(",
    ".getNumProperParts(",
    ".getNumChildren(",
    ".isRoot(",
    ".isLeaf(",
    ".isStrictAncestor(",
    ".isAlive(",
    ".getNeighborIndices(",
    ".getAltitude(",
    ".getNodeResidue(",
    ".getNumRowsOfGridDomain2D(",
    ".getNumColsOfGridDomain2D(",
    ".offset(",
    ".requireTopologyUnchanged(",
    ".requireMutationVersion(",
    ".dependencies.require(",
    ".dependencies.requireAll(",
    ".linearIndex(",
    ".getIndex(",
    ".at(",
    "DependencyResolver<",
    "findDependencySource(",
    "Computer::compute(",
    "EventEngine::",
    "BitquadLocalEventComputation::",
    "ContourSideLocalEventComputation::",
    "BitquadAttributeMaterialization::",
    "ContourSideAttributeMaterialization::",
    "ContoursComputedIncrementally::extract",
    "requireAttributeBufferShape(",
    "validateAltitudeBufferShape(",
)

KERNEL_NAMESPACE = re.compile(r"\bnamespace\s+(?:(?:[A-Za-z_][A-Za-z0-9_]*)::)*kernel\s*\{")


def kernel_namespace_bodies(text: str):
    for match in KERNEL_NAMESPACE.finditer(text):
        open_brace = text.find("{", match.start())
        depth = 0
        for index in range(open_brace, len(text)):
            if text[index] == "{":
                depth += 1
            elif text[index] == "}":
                depth -= 1
                if depth == 0:
                    yield text[open_brace : index + 1]
                    break


def main() -> int:
    root = pathlib.Path(sys.argv[1]).resolve()
    failures: list[str] = []
    paths = (path for path in (root / "mmcfilters").rglob("*") if path.suffix in {".h", ".hpp", ".cc", ".cpp"})
    for path in sorted(paths):
        text = path.read_text(encoding="utf-8")
        for body in kernel_namespace_bodies(text):
            for token in FORBIDDEN:
                if token in body:
                    failures.append(f"{path.relative_to(root)}: detail::kernel contains forbidden '{token}'")
    if failures:
        print("detail::kernel boundary audit failed:", file=sys.stderr)
        print("\n".join(failures), file=sys.stderr)
        return 1
    print("detail::kernel boundary audit passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
