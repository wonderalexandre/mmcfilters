#!/usr/bin/env python3

import pathlib
import sys

import numpy as np

from test_python_attribute_numeric_validation import load_native_module


def main():
    mmcfilters = load_native_module(pathlib.Path(sys.argv[1]).resolve())
    attribute = mmcfilters.Attribute
    image = np.ones((5, 5), dtype=np.uint8)
    image[1, 1] = 0
    tree = mmcfilters.MorphologicalTreeFactory.create_max_tree(image, 1.0)
    expected_names = [
        "FILLED_AREA", "FILLED_CENTROID_ROW", "FILLED_CENTROID_COLUMN",
        "FILLED_LENGTH_MAJOR_AXIS", "FILLED_LENGTH_MINOR_AXIS",
        "FILLED_AXIS_ORIENTATION", "FILLED_ECCENTRICITY", "FILLED_INERTIA",
        "HOLE_AREA_FRACTION", "FILLED_CENTROID_DISPLACEMENT_NORMALIZED",
        "FILLED_COMPACTNESS", "FILLED_CIRCULARITY",
    ]
    for dtype in (np.float32, np.float64):
        names, values = attribute.compute_attributes(tree, [attribute.Group.FILLED_SHAPE], dtype=dtype)
        assert list(names) == expected_names
        assert values.dtype == dtype
        string_names, string_values = attribute.compute_attributes(tree, ["FILLED_SHAPE"], dtype=dtype)
        assert names == string_names
        np.testing.assert_array_equal(values, string_values)
        for name in expected_names:
            scalar = getattr(attribute, name)
            single = attribute.compute_single_attribute(tree, scalar, dtype=dtype)
            np.testing.assert_array_equal(single, values[:, names[name]])
            exported = tree.project_node_values_to_exported_higra(single, scalar)
            unit = np.zeros(25, dtype=dtype)
            if name in ("FILLED_AREA", "FILLED_ECCENTRICITY", "FILLED_CIRCULARITY"):
                unit[:] = 1
            elif name == "FILLED_CENTROID_ROW":
                unit[:] = np.arange(25) // 5
            elif name == "FILLED_CENTROID_COLUMN":
                unit[:] = np.arange(25) % 5
            np.testing.assert_array_equal(exported[:25], unit)
        np.testing.assert_array_equal(values[:, names["FILLED_AREA"]], 25)
        np.testing.assert_array_equal(values[:, names["FILLED_CENTROID_ROW"]], 2)
        np.testing.assert_array_equal(values[:, names["FILLED_CENTROID_COLUMN"]], 2)
        np.testing.assert_allclose(values[:, names["FILLED_COMPACTNESS"]], 1 / (8 * np.pi), rtol=1e-6)
        np.testing.assert_array_equal(values[:, names["FILLED_CIRCULARITY"]], 1)
        support_area = attribute.compute_single_attribute(tree, attribute.AREA, dtype=dtype)
        ring = np.flatnonzero(support_area == 24)[0]
        np.testing.assert_allclose(values[ring, names["HOLE_AREA_FRACTION"]], 1 / 25, rtol=1e-6)
        np.testing.assert_allclose(values[ring, names["FILLED_CENTROID_DISPLACEMENT_NORMALIZED"]],
                                   np.sqrt(2) / (24 * 5), rtol=1e-6)
        for group in (attribute.Group.SHAPE, attribute.Group.MOMENTS, attribute.Group.BOUNDARY):
            other_names, _ = attribute.compute_attributes(tree, [group], dtype=dtype)
            assert not set(expected_names).intersection(other_names)
    requirements = attribute.requirements(attribute.FILLED_AREA)
    assert requirements["grid_domain_2d"]
    assert requirements["canonical_4_or_8_adjacency"]
    assert requirements["altitude_for_directional_adjacency"]
    print("Filled-shape names, groups, values, and exported pixel rows passed.")


if __name__ == "__main__":
    main()
