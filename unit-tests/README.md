# unit-tests

This directory contains the authoritative regression suite for the active `MorphologicalTree` API.

## Goals

- validate the current dense `NodeId` interface that should remain stable;
- cover proper-part ownership and structural tree operations;
- catch regressions in Python bindings, filters, attributes, and tree mutations.

## Current coverage

The suite covers:

- tree construction, traversal, topology, and mutations;
- proper parts, reconstruction, and tree-of-shapes behaviour;
- incremental attributes and secondary attribute consumers;
- filters, extinction values, contours, and Python bindings;
- structural invariants and parent-array round trips.

## Running the suite

Configure with tests enabled and then run:

```bash
cmake -S . -B build -DMMCFILTERS_BUILD_TESTS=ON
cmake --build build
ctest --test-dir build --output-on-failure
```
