# Contour Traces

`ContourTraceComputation` is the geometric contour API. It represents
boundaries as oriented pixel sides and traces those sides into ordered loops.
This preserves information that pixel-contour sets cannot represent, including
which loops are external and which loops bound holes.

For the older pixel-contour API, see [Incremental Contours](contours.md).

## C++ API

```cpp
auto traces = ContourTraceComputation::extract(tree);

for (ContourTraceEdge edge : traces.getEdges(nodeId)) {
    // edge.pixel is the row-major support pixel id.
    // edge.side is North, West, East, or South.
}

for (const ContourTraceLoop& loop : traces.getLoops(nodeId)) {
    if (loop.kind == ContourLoopKind::External) {
        // outer boundary
    } else {
        // inner boundary around a hole
    }

    for (ContourTraceEdge edge : traces.getLoopEdges(loop)) {
        // ordered loop edge
    }
}
```

`getEdges(node)` materializes only the unordered side-level boundary of the
requested node. `getLoops(node)` first materializes edges if needed, then traces
ordered loops for that node. `materializeAll()` traces every live node.

## Python API

```python
traces = mmcfilters.ContourTraceComputation.extraction(tree)

edges = traces.getEdges(node_id)
loops = traces.getLoops(node_id)

for loop in loops:
    loop_edges = traces.getLoopEdges(loop)
```

The Python API returns lists for node-local queries. The C++ owner keeps the
lazy caches internally.

## Geometry Convention

Each boundary primitive is a side of a support pixel. Edges are oriented in
image coordinates, with rows growing downward and columns growing rightward.
The support pixel stays on the right side of each directed edge.

With this convention:

- external loops have positive doubled signed area;
- internal loops have negative doubled signed area;
- `ContourLoopKind` is derived from the signed-area sign.

At grid vertices with several outgoing boundary edges, tracing uses a
deterministic wall-following rule: right turn, then straight, then left, then
back. This keeps diagonal-touching supports split into stable loop components.

## Compatibility With Scalar Attributes

For a node `v`, let `B_v` be the number of traced boundary edges.

- `B_v == CONTOUR_PERIMETER(v)`.
- Directional edge counts match `CONTOUR_SIDE_NORTH`, `CONTOUR_SIDE_WEST`,
  `CONTOUR_SIDE_EAST`, and `CONTOUR_SIDE_SOUTH`.
- Projecting traced edges to unique pixels matches the older
  `ContoursComputedIncrementally::getContour(v)` pixel set.

`BITQUADS_NUMBER_HOLES` can be used as a sanity check on simple connected
supports, but it has its own digital-connectivity convention and is not the
source of loop classification.

## Complexity

Let:

- `P` be the number of image pixels;
- `N` be the number of internal node slots;
- `B_v` be the number of boundary edges of node `v`;
- `M(S)` be the number of not-yet-materialized nodes in a requested subtree;
- `B(S)` be the number of boundary edges committed while materializing that
  subtree.

Extraction with valid tree-query caches is:

```text
O(N + P)
```

First edge materialization for a subtree is:

```text
O(M(S) + B(S) + D(S))
```

where `D(S)` is the number of compact local edge additions/removals read for
the newly materialized part of `S`. The global compact delta count is bounded
by a constant factor of `P` because the stencil has four sides per pixel.

Tracing one materialized node is:

```text
O(B_v)
```

Materializing and tracing all live nodes is output-sensitive:

```text
O(N + P + sum_v B_v)
```

Any implementation that flood-fills or reconstructs the full image domain
independently for every node would regress to `O(NP)` and should not be used
for this API.

## Benchmark

Build examples and run:

```bash
cmake -S . -B build -DMMCFILTERS_BUILD_EXAMPLES=ON -DMMCFILTERS_BUILD_PYTHON=OFF
cmake --build build --target mmcfilters_contour_trace_benchmark

./build/examples/mmcfilters_contour_trace_benchmark 512 512 5
./build/examples/mmcfilters_contour_trace_benchmark path/to/image.png 5
```

The benchmark reports extraction, root edge access, root loop tracing, all-node
edge materialization, all-node loop tracing, and explicit `materializeAll()`.
For internal optimization work, `mmcfilters_contour_trace_profile` profiles the
loop-only phase after edge caches are ready and separates adjacency
construction, loop walking, cache commit, and scratch release time. The
profiler also reports outgoing-vertex degree distribution and
successor-selection scan counts, which helps audit whether future optimization
work should target ambiguous vertices, degree-1 walking, or a different
directed-edge layout.

### Reference run on a real image

Local `RelWithDebInfo` run on 2026-06-01 with `dat/mesa.png` (`2455 x
1305`, `3,203,775` pixels), one repeat:

```bash
./build-relwithdebinfo/examples/mmcfilters_contour_benchmark dat/mesa.png 1
./build-relwithdebinfo/examples/mmcfilters_contour_trace_benchmark dat/mesa.png 1
```

Common max-tree workload:

| API | Extract | Materialization/tracing after extract | Approx. cache memory | Cached output |
| --- | ---: | ---: | ---: | --- |
| Legacy compact pixel contours | `0.31 s` | `0.27 s` for `materializeAll()` | `241.8 MiB` | `37.1 M` contour pixel ids |
| Trace edges only | `0.36 s` | `0.41 s` for all-node `getEdges()` | `403.9 MiB` | `51.4 M` side edges |
| Trace root loop | `0.37 s` | `0.38 s` for `getLoops(root)` | `404.4 MiB` | reuses `51.4 M` side edges, `1` loop |
| Trace all loops | `0.36 s` | `1.50 s` for all-node `getLoops()` random order | `442.0 MiB` | reuses `51.4 M` side edges, `2.31 M` loops |

Interpretation:

- compact contours remain cheaper when the workload only needs unordered
  boundary pixels;
- trace edges are the cheaper path when side-level geometry is enough;
- `getLoops(node)` is lazy at loop level: it materializes the edges needed for
  `node`, but traces loops only for that requested node;
- traced loops are substantially more expensive on large real images because
  they order every side edge and cache loop metadata for every live node;
- explicit `materializeAll()` traces all max-tree loops in tree order faster
  than random node-local loop access in this run (`1.15 s` after extraction);
- isolated node-local loop tracing uses a sparse outgoing-vertex table when it
  is smaller than the dense image-vertex table; repeated node-local tracing
  switches back to the dense table after a small threshold to avoid slowing
  down all-node workloads;
- loop tracing reuses materialized side-edge storage when safe, so the loop
  cache no longer persists a second full copy of ordered side edges;
- trace-only vertex scratch is allocated lazily, so extraction and edge-only
  workloads no longer retain the full image vertex-head table;
- once all edge caches are ready, local edge deltas and edge-mark scratch are
  released, reducing memory for point-loop and all-loop workloads;
- the trace API is the definitive geometric API when callers need
  external/internal loop separation, ordered traversal, or signed loop area.
