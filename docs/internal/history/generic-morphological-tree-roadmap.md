# Generic Morphological-Tree Roadmap

## Objective

Consolidate `MorphologicalTree` as the topology of a rooted hierarchy of
partial partitions over a finite proper-part set. Altitude, regular 2D geometry,
and adjacency are explicit capabilities layered on top of that topology.
Concrete families such as max-tree, min-tree, and Tree of Shapes are
producers of the generic representation rather than special cases embedded in
the core API.

This roadmap intentionally does not introduce a 3D domain.

## Invariants

- Internal nodes and proper parts have independent dense id domains.
- Every proper part has exactly one live direct owner.
- Every committed live node has non-empty full subtree support, although it may
  own no direct proper part.
- The committed topology is one connected rooted tree.
- A declared ordered `AltitudeOrder` is validated with strict parent-child
  inequalities at construction and commit.
- Geometry-dependent operations require an attached `GridDomain2D`.
- Algorithms dispatch on capabilities, never on `MorphologicalTreeKind`.
- `MorphologicalTreeKind` remains descriptive metadata.

## Phase 1 — Capability Contracts

Deliverables:

- declare altitude, regular-grid, adjacency, and monotone-order requirements for
  every public scalar attribute;
- centralize validation before scheduler execution;
- expose the declarations through C++ and Python metadata APIs;
- produce diagnostics naming the incompatible attribute and missing capability.

Completion criteria:

- abstract-domain `AREA`, tree-topology attributes, and `GRAY_HEIGHT` work;
- geometry attributes fail before their computer is entered;
- `MAX_DIST` accepts/rejects by capabilities independently of descriptive kind.

## Phase 2 — Safe Editing

Deliverables:

- retain linear checked commit as the only publication boundary;
- make every public editor recoverable with explicit rollback;
- restore topology and altitude when a public editor is abandoned;
- use a lazy copy-on-first-write journal proportional to the mutation delta;
- cover cycles, detached nodes, empty support, reused slots, failed altitude
  validation, repair, rollback, and successful commit.

Completion criteria:

- failed validation preserves staged-repair behavior while the editor lives;
- all public editors provide the strong rollback guarantee;
- `edit()` uses a recoverable delta journal without a full copy;
- no public method can publish an unchecked invalid hierarchy.

## Phase 3 — Producer Boundary

Deliverables:

- define one typed native hierarchy view containing parent, ownership, altitude,
  root, optional 2D domain, and semantics;
- route native factory overloads and concrete ToS builders through that
  representation;
- keep Higra import as an adapter because it uses a distinct leaf-first id
  domain.

Completion criteria:

- builder-specific results do not alter the `MorphologicalTree` data model;
- abstract and regular-2D producers use the same validated materialization path.

## Phase 4 — Attribute Classification

Deliverables:

- make attribute requirements the canonical classification mechanism;
- preserve the topology-versus-altitude execution distinction;
- represent conditional altitude requirements such as directional BitQuads;
- document a complete attribute/capability matrix.

Completion criteria:

- all registered attributes have tested metadata;
- group requests validate the expanded scalar attributes consistently.

## Phase 5 — API Stabilization

Deliverables:

- use `descriptiveKind`, `gridDomain2D`, and `HierarchySemantics` as canonical
  names;
- remove deprecated aliases after migrating internal callers;
- add C++ and Python examples for abstract and regular-2D native hierarchies;
- verify the installed-header consumer.

Completion criteria:

- public C++ and Python APIs expose the same capability model;
- the public surface has one spelling per capability.

## Phase 6 — Scientific Validation And Performance

Deliverables:

- add deterministic randomized/property tests for generic native trees;
- compare `GRAY_HEIGHT` with a brute-force subtree oracle;
- retain contrast-inversion and backend conformance checks;
- add a benchmark for native materialization and linear commit validation on
  deep and branching hierarchies.

Completion criteria:

- focused invariant/property tests pass;
- full C++/Python test suite passes;
- benchmark reports construction/validation time without embedding timing
  thresholds in correctness tests.

## Validation Commands

```bash
cmake --build /tmp/mmcfilters-tos-python-build -j4
ctest --test-dir /tmp/mmcfilters-tos-python-build --output-on-failure
python -m compileall -q python/mmcfilters unit-tests/python
git diff --check
```

## Implementation Status

Completed on 2026-07-23:

- all six phases above are implemented;
- the installed C++ consumer exercises the capability registry and the common
  native-producer boundary;
- deterministic randomized tests compare generic-tree `GRAY_HEIGHT` against a
  brute-force subtree oracle and verify contrast inversion;
- the complete configured test suite passes (`54/54`);
- Python bytecode compilation succeeds;
- `mmcfilters_tree_validation_benchmark` runs on deep and branching native
  hierarchies without correctness timing thresholds.

One intentionally deferred direction remains outside this roadmap: adding a 3D
proper-part domain. The generic topology is dimension-independent, while the
only concrete geometry capability exposed here is `GridDomain2D`.
