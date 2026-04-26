# Examples

This directory contains small standalone programs that exercise optional workflows without adding those dependencies to `mmcfilters::core`.

## Current examples

- `edt_dift_example.cpp`
  Runs `maxdist::EdtDIFT` on a simple synthetic shape and writes:
  - `binary.png`
  - `distance_transform.png`
- `editing_api_example.cpp`
  Demonstrates the public editing contract:
  - safe `pruneNode` and `mergeNodeIntoParent` calls;
  - staged topology edits through `TreeEditor`;
  - staged weighted edits through `WeightedTreeEditor`;
  - `commit()` rejecting invalid weighted altitude order.
- `contour_benchmark.cpp`
  Benchmarks incremental contour extraction/materialization for Component Tree
  and Tree of Shapes. It accepts either synthetic dimensions or an image path.

## Why this directory exists

`EdtDIFT` used to carry PNG export helpers directly in the core class. Those helpers required `stb_image_write`, which meant the core library needed an extra compiled support unit.

That dependency has now been moved out of the core:

- `mmcfilters::core` remains effectively header-only;
- `stb` is compiled only for example programs that need image export.

## Building examples

Configure with:

```bash
cmake -S . -B build -DMMCFILTERS_BUILD_EXAMPLES=ON
cmake --build build
```

The current executables are:

```bash
./build/examples/mmcfilters_example_edt_dift
./build/examples/mmcfilters_example_editing_api
./build/examples/mmcfilters_contour_benchmark 1024 1024 3
./build/examples/mmcfilters_contour_benchmark path/to/image.png 3
```
