# Examples

Build the example targets with:

```bash
cmake -S . -B build -DMMCFILTERS_BUILD_EXAMPLES=ON -DMMCFILTERS_BUILD_PYTHON=OFF
cmake --build build
```

Public examples use stable headers and are suitable as API references. Internal
examples may include `detail/` headers and are intended for validation,
benchmarking, and implementation debugging.

The Python filtering demo uses an installed package and does not require a
CMake example target:

```bash
python examples/python_attribute_filter_demo.py
```

## Public Targets

- `mmcfilters_example_editing_api`
  Demonstrates safe topology mutators, staged `TreeEditor` commits, weighted
  editor commits, and failed weighted validation. Run:
  `./build/examples/mmcfilters_example_editing_api`

- `mmcfilters_contour_benchmark`
  Measures compact contour extraction, lazy `getContour(node)` materialization,
  `contoursByNode()` iteration, random-order access, and `materializeAll()`.
  Run with a synthetic square image:
  `./build/examples/mmcfilters_contour_benchmark 512 512 5`
  or with a grayscale image file:
  `./build/examples/mmcfilters_contour_benchmark path/to/image.png 5`

- `mmcfilters_contour_trace_benchmark`
  Measures definitive geometric contour traces: extraction, side-level edge
  materialization, ordered loop tracing, random-order access, and
  `materializeAll()`. Run with a synthetic square image:
  `./build/examples/mmcfilters_contour_trace_benchmark 512 512 5`
  or with a grayscale image file:
  `./build/examples/mmcfilters_contour_trace_benchmark path/to/image.png 5`

- `mmcfilters_maxdist_benchmark`
  Measures `MAX_DIST` attribute computation on max-trees and min-trees with
  adjacency metadata. Run:
  `./build/examples/mmcfilters_maxdist_benchmark --sizes 128,256,512 --repeats 3 --radius 1.5`

- `mmcfilters_altitude_span_attribute_benchmark`
  Compares canonical `WeightedMorphologicalTree<uint8_t>` computation against
  `WeightedTreeView<T>` altitude-span computation for equivalent altitude
  buffers. Run:
  `./build/examples/mmcfilters_altitude_span_attribute_benchmark --sizes 128,256,512 --repeats 3 --suite core`

- `mmcfilters_typed_altitude_benchmark`
  Measures typed component-tree construction and attribute computation for
  `uint8_t`, `int32_t`, and `float` altitude owners. Output is CSV. Run:
  `./build/examples/mmcfilters_typed_altitude_benchmark --sizes 256,512 --repeats 3 --suite both`

## Internal Targets

- `mmcfilters_internal_edt_dift`
  Exercises the internal exact distance-transform support used by `MAX_DIST`
  and writes visualization PNGs. Run:
  `./build/examples/mmcfilters_internal_edt_dift tmp/edt_dift_example`

- `mmcfilters_local_events_bitquads_benchmark`
  Benchmarks local-event bitquad and contour-side computations against baseline
  attribute paths. Run:
  `./build/examples/mmcfilters_local_events_bitquads_benchmark 256 256 3 --tos`

- `mmcfilters_contour_trace_profile`
  Profiles contour trace loop materialization after edge caches are ready,
  separating adjacency construction, loop walking, cache commit, and scratch
  release time. It also reports outgoing-vertex degree distribution and
  successor-selection scans so loop-walking bottlenecks can be classified
  before changing the trace data layout. Run:
  `./build/examples/mmcfilters_contour_trace_profile path/to/image.png 3`

- `mmcfilters_tos_bitquad_projection_export`
  Exports Tree-of-Shapes bitquad projection CSV files for inspection. Run:
  `./build/examples/mmcfilters_tos_bitquad_projection_export --synthetic fixture --out-dir tmp/tos-bitquad-projections`

- `mmcfilters_tos_bitquad_naive_validation`
  Compares Tree-of-Shapes local-event bitquad data with a naive reference on a
  grayscale image. Run:
  `./build/examples/mmcfilters_tos_bitquad_naive_validation path/to/image.pgm`

## Input Contracts

The public examples follow the same contracts as the C++ API:

- component-tree benchmarks use positive image dimensions and positive repeat
  counts;
- `radius=1.0` selects 4-connectivity and `radius=1.5` selects 8-connectivity
  on the 2D square grid;
- `MAX_DIST` requires max-tree or min-tree component trees with adjacency
  metadata;
- typed-altitude examples use C++ `Image<T>` factories. Python remains limited
  to C-contiguous 2D `np.uint8` images.
