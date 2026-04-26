# Incremental Contours API

`ContoursComputedIncrementally` stores compact per-node contour deltas during
extraction and materializes final contours lazily when callers iterate them.

For implementation details, invariants, complexity, and benchmark
interpretation, see [contours-internals.md](contours-internals.md).

## C++

```cpp
auto contours = ContoursComputedIncrementally::extractCompactContours(tree);

for (int pixel : contours.getContour(nodeId)) {
    // use one contour pixel
}

for (auto [nodeId, contour] : contours.contoursByNode()) {
    for (int pixel : contour) {
        // use every node contour incrementally
    }
}
```

`getContour(nodeId)` returns a cache-aware range. The first iteration
materializes the requested subtree as needed; later iterations over already
materialized nodes only scan cached contiguous values.

`materializeAll()` is an explicit prefetch for workloads that will revisit many
contours repeatedly. It is not required for ordinary iteration.

## Python

```python
contours = mmcfilters.ContourComputation.extraction(tree)

for pixel in contours.getContour(node_id):
    pass

for node_id, contour in contours.contoursByNode():
    for pixel in contour:
        pass
```

The Python API uses the same contour access names: `getContour(node_id)` for
one node and `contoursByNode()` for all live nodes.

## Benchmark

Build examples and run:

```bash
cmake -S . -B build -DMMCFILTERS_BUILD_EXAMPLES=ON
cmake --build build --target mmcfilters_contour_benchmark

./build/examples/mmcfilters_contour_benchmark 1024 1024 3
./build/examples/mmcfilters_contour_benchmark path/to/image.png 3
```

The benchmark reports Component Tree and Tree of Shapes timings for extraction,
single-root access, full ordered iteration, full random-order iteration, and
explicit `materializeAll()` prefetch plus iteration.
