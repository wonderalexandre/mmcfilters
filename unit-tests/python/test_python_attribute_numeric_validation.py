#!/usr/bin/env python3

import importlib.util
import os
import pathlib
import sys

import numpy as np


def load_native_module(build_dir: pathlib.Path):
    package_dir = build_dir / "python" / "mmcfilters"
    package_init = package_dir / "__init__.py"
    if not package_init.is_file():
        raise RuntimeError(f"package init not found: {package_init}")

    for name in list(sys.modules):
        if name == "mmcfilters" or name.startswith("mmcfilters."):
            sys.modules.pop(name, None)

    def handles_mmcfilters(finder):
        known_modules = {}
        known_modules.update(getattr(finder, "known_source_files", {}))
        known_modules.update(getattr(finder, "known_wheel_files", {}))
        return any(name == "mmcfilters" or name.startswith("mmcfilters.") for name in known_modules)

    sys.meta_path = [finder for finder in sys.meta_path if not handles_mmcfilters(finder)]

    spec = importlib.util.spec_from_file_location(
        "mmcfilters",
        package_init,
        submodule_search_locations=[str(package_dir)],
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules["mmcfilters"] = module
    spec.loader.exec_module(module)
    return module


def require(condition: bool, message: str):
    if not condition:
        raise RuntimeError(message)


def require_raises(fn, message: str):
    try:
        fn()
    except Exception:
        return
    raise RuntimeError(message)


def read_pgm(path: pathlib.Path) -> np.ndarray:
    def next_token(handle):
        token = bytearray()
        while True:
            ch = handle.read(1)
            if not ch:
                raise RuntimeError(f"unexpected end of file while reading {path}")
            if ch == b"#":
                handle.readline()
                continue
            if ch.isspace():
                continue
            token.extend(ch)
            break

        while True:
            ch = handle.read(1)
            if not ch or ch.isspace():
                break
            if ch == b"#":
                handle.readline()
                break
            token.extend(ch)
        return bytes(token)

    with path.open("rb") as handle:
        magic = next_token(handle)
        cols = int(next_token(handle))
        rows = int(next_token(handle))
        max_value = int(next_token(handle))

        if magic == b"P5":
            dtype = np.uint8 if max_value < 256 else ">u2"
            expected = rows * cols * np.dtype(dtype).itemsize
            payload = handle.read(expected)
            require(len(payload) == expected, f"{path.name}: incomplete PGM payload")
            image = np.frombuffer(payload, dtype=dtype).reshape(rows, cols)
        elif magic == b"P2":
            values = []
            try:
                while True:
                    values.append(int(next_token(handle)))
            except RuntimeError:
                pass
            require(len(values) == rows * cols, f"{path.name}: unexpected PGM sample count")
            image = np.asarray(values, dtype=np.uint16).reshape(rows, cols)
        else:
            raise RuntimeError(f"{path.name}: unsupported PGM magic {magic!r}")

    if max_value != 255:
        image = np.rint(image.astype(np.float64) * (255.0 / float(max_value)))
    return np.ascontiguousarray(image.astype(np.uint8))


def sampled(image: np.ndarray, limit: int = 40) -> np.ndarray:
    row_step = max(1, image.shape[0] // limit)
    col_step = max(1, image.shape[1] // limit)
    return np.ascontiguousarray(image[::row_step, ::col_step][:limit, :limit])


def synthetic_images():
    unit = np.array([[9]], dtype=np.uint8)

    constant = np.full((8, 8), 37, dtype=np.uint8)

    concentric = np.full((11, 11), 20, dtype=np.uint8)
    concentric[2:9, 2:9] = 100
    concentric[4:7, 4:7] = 180

    horizontal = np.zeros((9, 13), dtype=np.uint8)
    horizontal[4, 1:12] = 220
    horizontal[:, 6] = np.maximum(horizontal[:, 6], 140)

    diagonal = np.zeros((12, 12), dtype=np.uint8)
    np.fill_diagonal(diagonal, 240)
    diagonal[np.arange(11), np.arange(1, 12)] = 120

    ring = np.zeros((14, 14), dtype=np.uint8)
    ring[2:12, 2:12] = 180
    ring[5:9, 5:9] = 30

    checker = np.indices((10, 10)).sum(axis=0) % 2
    checker = (checker * 190 + 30).astype(np.uint8)

    yield "unit", unit
    yield "constant", constant
    yield "concentric", concentric
    yield "cross-line", horizontal
    yield "diagonal", diagonal
    yield "ring", ring
    yield "checker", checker


def real_images(data_dir: pathlib.Path):
    for filename in ("lena.pgm", "brain2.pgm", "wrist.pgm"):
        path = data_dir / filename
        require(path.is_file(), f"missing real test image: {path}")
        yield path.stem, sampled(read_pgm(path))


def finite_examples(values: np.ndarray, layout=None) -> str:
    bad = np.argwhere(~np.isfinite(values))
    if bad.size == 0:
        return ""

    inverse_layout = {}
    if layout is not None:
        inverse_layout = {int(column): str(name) for name, column in layout.items()}

    examples = []
    for index in bad[:8]:
        if values.ndim == 1:
            row = int(index[0])
            examples.append(f"row {row}={values[row]!r}")
        else:
            row = int(index[0])
            col = int(index[1])
            name = inverse_layout.get(col, f"column {col}")
            examples.append(f"{name}@row {row}={values[row, col]!r}")
    return ", ".join(examples)


def require_dtype_finite(values, label: str, layout=None, expected_dtype=np.float32):
    array = np.asarray(values)
    dtype = np.dtype(expected_dtype)
    require(array.dtype == dtype, f"{label}: expected {dtype}, got {array.dtype}")
    examples = finite_examples(array, layout)
    require(not examples, f"{label}: non-finite values found: {examples}")


def require_float32_finite(values, label: str, layout=None):
    require_dtype_finite(values, label, layout, np.float32)


def require_same_layout(lhs, rhs, label: str):
    require(dict(lhs) == dict(rhs), f"{label}: layouts differ")


def require_float32_matches_float64_cast(float32_values, float64_values, label: str):
    lhs = np.asarray(float32_values)
    rhs64 = np.asarray(float64_values)
    require(lhs.dtype == np.float32, f"{label}: expected float32 lhs, got {lhs.dtype}")
    require(rhs64.dtype == np.float64, f"{label}: expected float64 rhs, got {rhs64.dtype}")
    require(lhs.shape == rhs64.shape, f"{label}: shapes differ: {lhs.shape} != {rhs64.shape}")
    rhs = rhs64.astype(np.float32)
    if np.array_equal(lhs, rhs):
        return

    mismatch = np.argwhere(lhs != rhs)
    examples = []
    for index in mismatch[:8]:
        key = tuple(int(i) for i in index)
        examples.append(f"{key}: {lhs[key]!r} != {rhs[key]!r}")
    raise RuntimeError(f"{label}: float32 output differs from cast float64 output: {', '.join(examples)}")


def require_finite_matrix(values, label: str, layout=None, expected_dtype=np.float32):
    require(values.ndim == 2, f"{label}: expected a 2D attribute matrix")
    require(values.shape[0] > 0, f"{label}: expected at least one output row")
    if layout is not None:
        require(values.shape[1] == len(layout), f"{label}: layout size must match matrix width")
    require_dtype_finite(values, label, layout, expected_dtype)


def require_attribute_matrix(tree, layout, values, label: str, expected_dtype=np.float32):
    require_finite_matrix(values, label, layout, expected_dtype)

    for key in (
        "ECCENTRICITY",
        "BITQUADS_CIRCULARITY",
        "BITQUADS_PERIMETER_AVERAGE",
        "BITQUADS_LENGTH_AVERAGE",
        "BITQUADS_WIDTH_AVERAGE",
    ):
        require(key in layout, f"{label}: missing {key}")

    expected_rows = tree.numInternalNodeSlots
    require(values.shape[0] == expected_rows, f"{label}: expected {expected_rows} internal-node rows")


def validate_tree_dtype(mmcfilters, tree, label: str, dtype):
    dtype_label = np.dtype(dtype).name
    all_layout, all_values = mmcfilters.Attribute.computeAttributes(
        tree,
        [mmcfilters.Attribute.Group.ALL],
        dtype=dtype,
    )
    require_attribute_matrix(tree, all_layout, all_values, f"{label}: {dtype_label}: Attribute.Group.ALL", dtype)

    topology_layout, topology_values = mmcfilters.Attribute.computeTopologyAttributes(
        tree,
        [
            mmcfilters.Attribute.Group.MOMENTS,
            mmcfilters.Attribute.Group.BOUNDARY,
            mmcfilters.Attribute.Group.TREE_TOPOLOGY,
        ],
        dtype=dtype,
    )
    require_attribute_matrix(tree, topology_layout, topology_values, f"{label}: {dtype_label}: topology groups", dtype)

    for attr_name in (
        "ECCENTRICITY",
        "BITQUADS_CIRCULARITY",
        "BITQUADS_PERIMETER_AVERAGE",
        "BITQUADS_LENGTH_AVERAGE",
        "BITQUADS_WIDTH_AVERAGE",
    ):
        attr = getattr(mmcfilters.Attribute, attr_name)
        single = mmcfilters.Attribute.computeSingleAttribute(tree, attr, dtype=dtype)
        require(single.shape == (tree.numInternalNodeSlots,), f"{label}: {attr_name} single shape")
        require_dtype_finite(single, f"{label}: {dtype_label}: single {attr_name}", expected_dtype=dtype)

        topology_single = mmcfilters.Attribute.computeSingleTopologyAttribute(tree, attr, dtype=dtype)
        require(topology_single.shape == (tree.numInternalNodeSlots,), f"{label}: {attr_name} topology single shape")
        require_dtype_finite(topology_single, f"{label}: {dtype_label}: topology single {attr_name}", expected_dtype=dtype)

    for attr_name in ("ECCENTRICITY", "BITQUADS_WIDTH_AVERAGE", "LEVEL"):
        attr = getattr(mmcfilters.Attribute, attr_name)
        delta_layout, delta_values = mmcfilters.Attribute.computeSingleAttributeWithDelta(
            tree,
            attr,
            1,
            "last-padding",
            dtype=dtype,
        )
        require_finite_matrix(delta_values, f"{label}: {dtype_label}: delta {attr_name}", delta_layout, dtype)
        require(delta_values.shape[0] == tree.numInternalNodeSlots, f"{label}: delta {attr_name} row count")

    for attr_name in ("ECCENTRICITY", "BITQUADS_WIDTH_AVERAGE", "MAX_DIST"):
        attr = getattr(mmcfilters.Attribute, attr_name)
        mapped = mmcfilters.Attribute.computeAttributeMapping(tree, attr, dtype=dtype)
        require(mapped.shape == (tree.numRows, tree.numCols), f"{label}: mapped {attr_name} shape")
        require_dtype_finite(mapped, f"{label}: {dtype_label}: mapped {attr_name}", expected_dtype=dtype)

    projection_layout, projection_values = mmcfilters.Attribute.computeAttributes(
        tree,
        [mmcfilters.Attribute.AREA, mmcfilters.Attribute.MAX_DIST],
        dtype=dtype,
    )
    projected = tree.project_node_values_to_exported_higra(
        projection_values,
        [mmcfilters.Attribute.AREA, mmcfilters.Attribute.MAX_DIST],
    )
    require(projected.ndim == 2, f"{label}: {dtype_label}: exported projection must be 2D")
    require(projected.shape[1] == len(projection_layout), f"{label}: {dtype_label}: exported projection width")
    require_dtype_finite(projected, f"{label}: {dtype_label}: exported projection", expected_dtype=dtype)


def validate_tree(mmcfilters, tree, label: str):
    validate_tree_dtype(mmcfilters, tree, label, np.float32)
    validate_tree_dtype(mmcfilters, tree, label, np.float64)

    default_single = mmcfilters.Attribute.computeSingleAttribute(tree, mmcfilters.Attribute.ECCENTRICITY)
    require_float32_finite(default_single, f"{label}: default dtype")

    all32_layout, all32 = mmcfilters.Attribute.computeAttributes(
        tree,
        [mmcfilters.Attribute.Group.ALL],
        dtype=np.float32,
    )
    all64_layout, all64 = mmcfilters.Attribute.computeAttributes(
        tree,
        [mmcfilters.Attribute.Group.ALL],
        dtype=np.float64,
    )
    require_same_layout(all32_layout, all64_layout, f"{label}: Attribute.Group.ALL cast contract")
    require_float32_matches_float64_cast(all32, all64, f"{label}: Attribute.Group.ALL cast contract")

    topology32_layout, topology32 = mmcfilters.Attribute.computeTopologyAttributes(
        tree,
        [
            mmcfilters.Attribute.Group.MOMENTS,
            mmcfilters.Attribute.Group.BOUNDARY,
            mmcfilters.Attribute.Group.TREE_TOPOLOGY,
        ],
        dtype=np.float32,
    )
    topology64_layout, topology64 = mmcfilters.Attribute.computeTopologyAttributes(
        tree,
        [
            mmcfilters.Attribute.Group.MOMENTS,
            mmcfilters.Attribute.Group.BOUNDARY,
            mmcfilters.Attribute.Group.TREE_TOPOLOGY,
        ],
        dtype=np.float64,
    )
    require_same_layout(topology32_layout, topology64_layout, f"{label}: topology cast contract")
    require_float32_matches_float64_cast(topology32, topology64, f"{label}: topology cast contract")

    for attr_name in ("ECCENTRICITY", "BITQUADS_WIDTH_AVERAGE", "LEVEL", "MAX_DIST"):
        attr = getattr(mmcfilters.Attribute, attr_name)

        single32 = mmcfilters.Attribute.computeSingleAttribute(tree, attr, dtype=np.float32)
        single64 = mmcfilters.Attribute.computeSingleAttribute(tree, attr, dtype=np.float64)
        require_float32_matches_float64_cast(single32, single64, f"{label}: single {attr_name} cast contract")

        mapped32 = mmcfilters.Attribute.computeAttributeMapping(tree, attr, dtype=np.float32)
        mapped64 = mmcfilters.Attribute.computeAttributeMapping(tree, attr, dtype=np.float64)
        require_float32_matches_float64_cast(mapped32, mapped64, f"{label}: mapped {attr_name} cast contract")

    delta32_layout, delta32 = mmcfilters.Attribute.computeSingleAttributeWithDelta(
        tree,
        mmcfilters.Attribute.ECCENTRICITY,
        1,
        "last-padding",
        dtype=np.float32,
    )
    delta64_layout, delta64 = mmcfilters.Attribute.computeSingleAttributeWithDelta(
        tree,
        mmcfilters.Attribute.ECCENTRICITY,
        1,
        "last-padding",
        dtype=np.float64,
    )
    require_same_layout(delta32_layout, delta64_layout, f"{label}: delta cast contract")
    require_float32_matches_float64_cast(delta32, delta64, f"{label}: delta cast contract")

    projection32_layout, projection32_values = mmcfilters.Attribute.computeAttributes(
        tree,
        [mmcfilters.Attribute.AREA, mmcfilters.Attribute.MAX_DIST],
        dtype=np.float32,
    )
    projection64_layout, projection64_values = mmcfilters.Attribute.computeAttributes(
        tree,
        [mmcfilters.Attribute.AREA, mmcfilters.Attribute.MAX_DIST],
        dtype=np.float64,
    )
    require_same_layout(projection32_layout, projection64_layout, f"{label}: exported projection input layout")
    projected32 = tree.project_node_values_to_exported_higra(
        projection32_values,
        [mmcfilters.Attribute.AREA, mmcfilters.Attribute.MAX_DIST],
    )
    projected64 = tree.project_node_values_to_exported_higra(
        projection64_values,
        [mmcfilters.Attribute.AREA, mmcfilters.Attribute.MAX_DIST],
    )
    require_float32_matches_float64_cast(projected32, projected64, f"{label}: exported projection cast contract")


def validate_higra_roundtrip(mmcfilters, tree, label: str, kind):
    parent, altitude = tree.exportHigraHierarchy()
    imported = mmcfilters.MorphologicalTreeFactory.createFromHigraParent(
        parent,
        altitude,
        tree.numRows,
        tree.numCols,
        kind,
        1.5,
    )

    for dtype in (np.float32, np.float64):
        dtype_label = np.dtype(dtype).name
        layout, values = mmcfilters.Attribute.computeAttributes(
            imported,
            [mmcfilters.Attribute.Group.ALL],
            mmcfilters.NodeIdSpace.HIGRA,
            dtype=dtype,
        )
        require(values.ndim == 2, f"{label}: {dtype_label}: HIGRA ALL must be a 2D matrix")
        require(values.shape == (imported.numHigraNodes, len(layout)), f"{label}: {dtype_label}: HIGRA ALL shape")
        require_dtype_finite(values, f"{label}: {dtype_label}: HIGRA Attribute.Group.ALL", layout, dtype)


def validate_dtype_errors(mmcfilters):
    image = np.array([[0, 1], [2, 3]], dtype=np.uint8)
    tree = mmcfilters.MorphologicalTreeFactory.createMaxTree(image, 1.0)

    for bad_dtype in (np.float16, np.int32):
        require_raises(
            lambda bad_dtype=bad_dtype: mmcfilters.Attribute.computeSingleAttribute(
                tree, mmcfilters.Attribute.AREA, dtype=bad_dtype
            ),
            f"computeSingleAttribute must reject {bad_dtype}",
        )
        require_raises(
            lambda bad_dtype=bad_dtype: mmcfilters.Attribute.computeAttributes(
                tree, [mmcfilters.Attribute.AREA], dtype=bad_dtype
            ),
            f"computeAttributes must reject {bad_dtype}",
        )
        require_raises(
            lambda bad_dtype=bad_dtype: mmcfilters.Attribute.computeTopologyAttributes(
                tree, [mmcfilters.Attribute.AREA], dtype=bad_dtype
            ),
            f"computeTopologyAttributes must reject {bad_dtype}",
        )
        require_raises(
            lambda bad_dtype=bad_dtype: mmcfilters.Attribute.computeSingleAttributeWithDelta(
                tree, mmcfilters.Attribute.AREA, 1, dtype=bad_dtype
            ),
            f"computeSingleAttributeWithDelta must reject {bad_dtype}",
        )
        require_raises(
            lambda bad_dtype=bad_dtype: mmcfilters.Attribute.computeAttributeMapping(
                tree, mmcfilters.Attribute.AREA, dtype=bad_dtype
            ),
            f"computeAttributeMapping must reject {bad_dtype}",
        )


def validate_image(mmcfilters, name: str, image: np.ndarray, radius: float):
    require(image.dtype == np.uint8, f"{name}: image must be uint8")
    require(image.flags.c_contiguous, f"{name}: image must be contiguous")

    max_tree = mmcfilters.MorphologicalTreeFactory.createMaxTree(image, radius)
    validate_tree(mmcfilters, max_tree, f"{name}: max-tree r={radius:g}")
    validate_higra_roundtrip(
        mmcfilters,
        max_tree,
        f"{name}: max-tree r={radius:g}",
        mmcfilters.MorphologicalTreeKind.MAX_TREE,
    )

    min_tree = mmcfilters.MorphologicalTreeFactory.createMinTree(image, radius)
    validate_tree(mmcfilters, min_tree, f"{name}: min-tree r={radius:g}")
    validate_higra_roundtrip(
        mmcfilters,
        min_tree,
        f"{name}: min-tree r={radius:g}",
        mmcfilters.MorphologicalTreeKind.MIN_TREE,
    )


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: test_python_attribute_numeric_validation.py <build-dir>")

    build_dir = pathlib.Path(sys.argv[1]).resolve()
    data_dir = pathlib.Path(
        os.environ.get("MMCFILTERS_TEST_DATA_DIR", pathlib.Path(__file__).resolve().parents[2] / "dat")
    )

    mmcfilters = load_native_module(build_dir)
    validate_dtype_errors(mmcfilters)
    fixtures = list(synthetic_images()) + list(real_images(data_dir))
    require(len(fixtures) >= 10, "expected synthetic and real numeric validation fixtures")

    for name, image in fixtures:
        for radius in (1.0, 1.5):
            validate_image(mmcfilters, name, image, radius)

    print("python attribute numeric validation ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
