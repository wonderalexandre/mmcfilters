# SIBGRAPI 2026 reproducibility

This directory documents the experiments for *Unifying Local Attribute Computation on Component Trees and the Tree of Shapes*. All SIBGRAPI 2026 reproducibility documentation is intentionally kept in this single file.

## Camera-ready artifact

- implementation base: `v4.3.0` (`cd531b8121057395d6255232f4d50dc3af7fbd23`);
- experiment-code commit: `c587d870802355766abb825aa4a9948e0f24e110`;
- MorphoTree reference commit: `da32cf7666a774c25d11dc0200a63ebb3f1fe574`.

[`reference-results/`](../../../benchmarks/sibgrapi2026/reference-results/) contains the complete campaign used by the camera-ready paper: 5,400 timed measurements from 100 images, three resolutions, three hierarchies, two methods, and exactly three repetitions. Its manifest records `dirty: false`, base release `v4.3.0`, no exact tag at measurement time, all 13 compiled MorphoTree source hashes, and the runner hash.

## Implementation map

The paper implementation is not confined to the benchmark directory. Its generic algorithm, bitquad specialization, public-library integration, validation, and experiment harness are separated as follows.

| Paper role | Implementation |
|---|---|
| Localization of a translated sample on an anchor branch through the inclusion join/LCA | [`ConnectedSubsetTreeLocalizer.hpp`](../../../mmcfilters/localAttributes/ConnectedSubsetTreeLocalizer.hpp) |
| Finite observation windows, visibility states, event deltas, node-local increments, final node attributes, and the decision/algebra concepts | [`LocalEventModel.hpp`](../../../mmcfilters/localAttributes/LocalEventModel.hpp) |
| Compilation of local state changes into hierarchy-attached event deltas and dense increments | [`FiniteWindowLocalEventCompiler.hpp`](../../../mmcfilters/localAttributes/FiniteWindowLocalEventCompiler.hpp) |
| Generic bottom-up aggregation and the complete finite-window API shared by component trees and trees of shapes | [`FiniteWindowLocalAttributeComputer.hpp`](../../../mmcfilters/localAttributes/FiniteWindowLocalAttributeComputer.hpp) |
| Canonical `2x2` bitquad windows, ownership rule, local decisions, event algebra, signed increments, and aggregated `Q1/Q2/QD/Q3/Q4` counts | [`BitquadFiniteWindowComputation.hpp`](../../../mmcfilters/attributes/computers/detail/BitquadFiniteWindowComputation.hpp) |
| Canonical bitquad state/family types and lookup convention | [`BitquadAttributeData.hpp`](../../../mmcfilters/attributes/computers/detail/BitquadAttributeData.hpp) |
| Connectivity and lower/upper shape policy used when converting family counts to scalar descriptors | [`BitquadConnectivityPolicy.hpp`](../../../mmcfilters/attributes/computers/detail/BitquadConnectivityPolicy.hpp) and [`BitquadAttributeProjection.hpp`](../../../mmcfilters/attributes/computers/detail/BitquadAttributeProjection.hpp) |
| Adapter that exposes the implementation through the built-in `BITQUAD_*` attributes | [`BitquadAttributeComputer.hpp`](../../../mmcfilters/attributes/computers/BitquadAttributeComputer.hpp) |
| Public attribute metadata and backend dispatch | [`AttributeRegistry.hpp`](../../../mmcfilters/attributes/AttributeRegistry.hpp) and [`TopologyAttributeBackend.hpp`](../../../mmcfilters/attributes/detail/TopologyAttributeBackend.hpp) |

The principal implementation path is therefore:

```text
BitquadAttributeComputer
  -> BitquadFiniteWindowComputation
  -> FiniteWindowLocalEventCompiler
  -> ConnectedSubsetTreeLocalizer / inclusion join
  -> bottom-up aggregation
  -> BitquadAttributeProjection
```

The detailed mathematical-to-C++ correspondence, data invariants, complexity, and extension example are documented in the [finite-window local-attribute guide](../../finite-window-local-attributes.md).

## Validation and experiment code

The permanent library tests are in [`test_finite_window_local_attribute_computations.cpp`](../../../unit-tests/attributes/test_finite_window_local_attribute_computations.cpp). They cover the generic event pipeline, randomized comparison with a direct implementation, canonical bitquad states, component trees, and both complementary-connectivity conventions for trees of shapes.

The paper-specific code is under [`benchmarks/sibgrapi2026/`](../../../benchmarks/sibgrapi2026/):

- [`bitquad_exactness_validation.cpp`](../../../benchmarks/sibgrapi2026/bitquad_exactness_validation.cpp) validates the proposed computation against a direct exhaustive scan;
- [`bitquad_benchmark.cpp`](../../../benchmarks/sibgrapi2026/bitquad_benchmark.cpp) measures the proposed method and the two paper baselines;
- [`run_experiment.py`](../../../benchmarks/sibgrapi2026/run_experiment.py) validates inputs and provenance and executes the fixed protocol;
- [`analyze_results.py`](../../../benchmarks/sibgrapi2026/analyze_results.py) validates the raw measurements and generates the camera-ready table;
- [`benchmarks/CMakeLists.txt`](../../../benchmarks/CMakeLists.txt) defines the exactness target and the optional target that compiles the unmodified reference [5,6] snapshot.

For max-trees and min-trees, the comparison baseline is the specialized component-tree bitquad algorithm of Silva et al. [6]. For the tree of shapes, it is the original public implementation of da Silva et al. [5].

## Dependencies

The build requires a C++20 compiler, CMake 3.20 or later, and Python 3.9 or later. The Python scripts use only the standard library.

The published implementations used as references [5] and [6] are compiled directly from the unmodified [MorphoTree repository](https://github.com/dennisjosesilva/morphotree) at commit:

```text
da32cf7666a774c25d11dc0200a63ebb3f1fe574
```

Prepare that source snapshot separately:

```bash
git clone https://github.com/dennisjosesilva/morphotree.git morphotree-reference5
git -C morphotree-reference5 checkout da32cf7666a774c25d11dc0200a63ebb3f1fe574
git -C morphotree-reference5 status --porcelain
```

The last command must produce no output. The benchmark build uses the original source files directly; it does not copy, patch, or replace either algorithm. In the MorphoTree filenames, `dt-max-tree-8c.dat` and `dt-min-tree-8c.dat` are decision tables used by reference [6].

## Exactness validation

The validator enumerates all 512 binary `3x3` images. For each image it builds an 8-connected max-tree, an 8-connected min-tree, and a 4/8 tree of shapes. At every live node, it reconstructs the node support, scans every framed `2x2` cell, classifies its binary pattern as `Q1`, `Q2`, `QD`, `Q3`, or `Q4`, and compares the result with the proposed finite-window computation.

Build and run the validator in checked mode:

```bash
cmake -S . -B build-sibgrapi2026-check \
  -DCMAKE_BUILD_TYPE=Release \
  -DMMCFILTERS_BUILD_PYTHON=OFF \
  -DMMCFILTERS_BUILD_BENCHMARKS=ON \
  -DMMCFILTERS_CONTRACT_MODE=CHECKED
cmake --build build-sibgrapi2026-check \
  --target mmcfilters_sibgrapi2026_bitquad_exactness
./build-sibgrapi2026-check/benchmarks/mmcfilters_sibgrapi2026_bitquad_exactness
```

A successful run reports `binary_patterns=512`, `hierarchy_instances=1536`, `validated_nodes=3887`, `mismatched_nodes=0`, and `maximum_absolute_family_error=0`. This is a correctness experiment, not a performance benchmark: its direct support reconstruction is deliberately exhaustive and suitable only for small validation inputs.

## Timing build

Configure an optimized unchecked build, supplying the unmodified MorphoTree snapshot:

```bash
cmake -S . -B build-sibgrapi2026 \
  -DCMAKE_BUILD_TYPE=Release \
  -DMMCFILTERS_BUILD_PYTHON=OFF \
  -DMMCFILTERS_BUILD_BENCHMARKS=ON \
  -DMMCFILTERS_CONTRACT_MODE=UNCHECKED \
  -DMMCFILTERS_SIBGRAPI2026_MORPHOTREE_SOURCE=/path/to/morphotree-reference5
cmake --build build-sibgrapi2026 \
  --target mmcfilters_sibgrapi2026_bitquad_benchmark
```

## Dataset and fixed protocol

The experiment uses the first 100 validation images, named `val_000.png` through `val_099.png`, from the [Occluded RoadText Challenge](https://rrc.cvc.uab.es/?ch=29) at the following supplied resolutions:

- `480p`: `853 x 480`;
- `720p`: `1280 x 720`;
- `1080p`: `1920 x 1080`.

The runner never resizes images. It checks every PNG dimension and records each input SHA-256 digest in `dataset-manifest.csv`.

For every image and hierarchy, both methods operate on the same tree topology. Tree construction, topology import, and decision-table loading are outside the timed region. Each method is executed once as an untimed warm-up and then exactly three times. Within each paired repetition, execution order alternates deterministically; the starting method is balanced across images and hierarchies. The table first averages the three repetitions for each image and then averages the 100 image means. No outlier removal, replacement run, or adaptive resampling is performed.

Run the complete campaign in a new output directory:

```bash
python3 benchmarks/sibgrapi2026/run_experiment.py \
  --runner build-sibgrapi2026/benchmarks/mmcfilters_sibgrapi2026_bitquad_benchmark \
  --morphotree-source /path/to/morphotree-reference5 \
  --dataset-480p /path/to/icdar_480p \
  --dataset-720p /path/to/icdar_720p \
  --dataset-1080p /path/to/icdar_1080p \
  --output-dir benchmarks/results/sibgrapi2026
```

The mmcfilters checkout must descend from `v4.3.0` and must be clean. The output directory must be absent or empty. These checks prevent untracked code changes and accidental mixing of campaigns. `--allow-dirty` exists only for explicitly marked development runs.

## Camera-ready benchmark results

The camera-ready table comes from [`reference-results/table.tex`](../../../benchmarks/sibgrapi2026/reference-results/table.tex). The values below are mean warm attribute-API times over 100 images, in milliseconds. Each image mean contains exactly three timed repetitions. Tree construction, topology import, and decision-table loading are excluded.

| Tree | Method | `853 x 480` | `1280 x 720` | `1920 x 1080` |
|---|---|---:|---:|---:|
| Max-tree | Proposed | 43.33 | 95.60 | 226.53 |
| Max-tree | Silva et al. [6] | **41.60** | **87.09** | **215.13** |
| Min-tree | Proposed | 43.24 | 94.82 | 223.69 |
| Min-tree | Silva et al. [6] | **41.02** | **85.74** | **211.36** |
| Tree of Shapes | Proposed | **49.62** | **103.62** | **266.53** |
| Tree of Shapes | da Silva et al. [5], original | 603.59 | 1338.33 | 3077.76 |

For component trees, the proposed generic method took `4.15%`--`10.59%` more time than the specialized [6] implementation and agreed exactly at every compared node. For the Tree of Shapes, it was `12.16`, `12.92`, and `11.55` times faster at `480p`, `720p`, and `1080p`, respectively. As explained above, [5] is a computational-time baseline only.

## Generated artifacts

- `experiment.json`: code versions, protocol, environment, and symbolic commands;
- `dataset-manifest.csv`: dimensions and SHA-256 digests without absolute dataset paths;
- `raw-480p.csv`, `raw-720p.csv`, and `raw-1080p.csv`: all three timed repetitions and audit fields;
- `run-480p.log`, `run-720p.log`, and `run-1080p.log`: progress logs;
- `summary.csv`: paper means and between-image variation;
- `speedups.csv`: baseline/proposed ratios and relative differences;
- `correctness.csv`: exact component-tree agreement and the recorded [5] `Q2/QD` differences;
- `table.tex`: a camera-ready LaTeX table generated only after the raw campaign passes all checks.

## Complexity

For a fixed `2x2` bitquad window, the proposed attribute phase has `O(P + N)` time and `O(N)` auxiliary storage after tree construction, where `P` is the number of pixels and `N` the number of tree nodes. On the tree of shapes, this bound assumes the constant-time lowest-common-ancestor support prepared by the method.
