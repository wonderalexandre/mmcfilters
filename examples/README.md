# Examples

This directory contains small standalone programs that exercise optional workflows without adding those dependencies to `mmcfilters::core`.

## Current example

- `edt_dift_example.cpp`
  Runs `maxdist::EdtDIFT` on a simple synthetic shape and writes:
  - `binary.png`
  - `distance_transform.png`

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

The current executable is:

```bash
./build/examples/mmcfilters_example_edt_dift
```
