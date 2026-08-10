# Validation-sensitive scientific benchmarks

This focused protocol is public for reproducibility, but it is a scientific
maintenance tool rather than part of the public API contract.

`mmcfilters_validation_sensitive_algorithms_benchmark` measures public scientific
operations whose implementations traverse an already established hierarchy. Its
purpose is to compare `CHECKED` and `UNCHECKED` builds without changing the
public API or giving either build a different mathematical problem.

## Established input state

Before any timed sample, the benchmark deterministically creates:

- one synthetic `uint8` image;
- its max-tree with the standard radius-1.5 grid adjacency;
- one dense `AREA` attribute buffer indexed by internal `NodeId`.

Tree construction and the shared `AREA` computation are setup costs. They are
excluded because the measured question is the cost of consuming a valid,
already established hierarchy, not the cost of producing it.

## Timed scenarios

Each sample invokes the public API for one complete operation:

- Viterbi attribute filtering;
- ultimate attribute opening and both result images;
- extinction-value computation followed by top-k reconstruction;
- MSER selection from an altitude window;
- stable-region selection from a topological-depth window;
- delta-augmented `AREA` materialisation;
- joint bitquad area, perimeter, and circularity materialisation;
- compact Higra hierarchy export.

Each scenario has one untimed warm-up and reports the median of the requested
number of repetitions. Object construction that belongs to the public operation
is included in the timed sample. Result hashing and result destruction after the
public call are excluded.

## Correctness condition

Every numeric result is hashed bit for bit with FNV-1a. Attribute results also
include their semantic layout and node-id space in the hash. A scenario fails if
repetitions inside one executable produce different hashes. A valid comparison
additionally requires every scenario hash from the `CHECKED` executable to equal
the corresponding hash from the `UNCHECKED` executable.

Build both configurations in release mode, then use the comparison runner:

```text
cmake -S . -B build-benchmark-checked -DCMAKE_BUILD_TYPE=Release \
  -DMMCFILTERS_BUILD_BENCHMARKS=ON -DMMCFILTERS_CONTRACT_MODE=CHECKED
cmake --build build-benchmark-checked --target mmcfilters_validation_sensitive_algorithms_benchmark

cmake -S . -B build-benchmark-unchecked -DCMAKE_BUILD_TYPE=Release \
  -DMMCFILTERS_BUILD_BENCHMARKS=ON -DMMCFILTERS_CONTRACT_MODE=UNCHECKED
cmake --build build-benchmark-unchecked --target mmcfilters_validation_sensitive_algorithms_benchmark

python3 benchmarks/compare_validation_modes.py \
  build-benchmark-checked/benchmarks/mmcfilters_validation_sensitive_algorithms_benchmark \
  build-benchmark-unchecked/benchmarks/mmcfilters_validation_sensitive_algorithms_benchmark \
  --rows 256 --cols 256 --repetitions 9 --process-runs 5
```

The runner alternates executable order, rejects any checksum mismatch, and
reports the median of the per-process medians. A negative relative percentage
means that `UNCHECKED` was faster; a positive value means it was slower. For
reported experiments, record
the compiler, flags, CPU, operating system, image dimensions, repetition count,
process-run count, timings, and hashes.
