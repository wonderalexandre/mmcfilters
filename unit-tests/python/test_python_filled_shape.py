#!/usr/bin/env python3

import pathlib
import sys

import numpy as np

from test_python_attribute_numeric_validation import load_native_module


def verify_signed_orientations(mmcfilters):
    attribute = mmcfilters.Attribute

    def fill(mask, foreground_four):
        occupied = np.pad(mask, 1)
        outside = np.zeros_like(occupied)
        outside[0, 0] = True
        queue = [(0, 0)]
        for y, x in queue:
            for dy in (-1, 0, 1):
                for dx in (-1, 0, 1):
                    if (dx == dy == 0) or (not foreground_four and dx != 0 and dy != 0):
                        continue
                    yy, xx = y + dy, x + dx
                    if 0 <= yy < occupied.shape[0] and 0 <= xx < occupied.shape[1]:
                        if not occupied[yy, xx] and not outside[yy, xx]:
                            outside[yy, xx] = True
                            queue.append((yy, xx))
        return ~outside[1:-1, 1:-1]

    def angle(mask):
        y, x = np.nonzero(mask)
        x = x.astype(float) - x.mean()
        y = y.astype(float) - y.mean()
        delta, cross = x @ x - y @ y, 2 * (x @ y)
        return np.degrees(0.5 * np.arctan2(cross, delta)) if delta != 0 or cross != 0 else 0.0

    diagonal = np.zeros((9, 9), dtype=np.uint8)
    for row in range(1, 8):
        diagonal[row, row - 1:row + 2] = 1
    rectangle = np.zeros((9, 13), dtype=np.uint8)
    rectangle[1:8, 1:12] = 1
    hole = rectangle.copy()
    hole[2:4, 2:5] = 0
    horizontal = np.zeros((5, 9), dtype=np.uint8)
    horizontal[2, 1:8] = 1
    square = np.pad(np.ones((3, 3), dtype=np.uint8), 1)
    sources = [diagonal, rectangle, hole, horizontal, square, np.ones((1, 1), dtype=np.uint8)]
    requested = ["AXIS_ORIENTATION_SIGNED", "AXIS_ORIENTATION",
                 "FILLED_AXIS_ORIENTATION_SIGNED", "FILLED_AXIS_ORIENTATION"]
    for source in sources:
        for rotation in range(4):
            image = np.ascontiguousarray(np.rot90(source, rotation))
            for foreground_four, radius in ((True, 1.0), (False, 1.5)):
                for maximum in (True, False):
                    factory = (mmcfilters.MorphologicalTreeFactory.create_max_tree if maximum
                               else mmcfilters.MorphologicalTreeFactory.create_min_tree)
                    tree = factory(image if maximum else 1 - image, radius)
                    for dtype in (np.float32, np.float64):
                        names, values = attribute.compute_attributes(tree, requested, dtype=dtype)
                        for signed, unsigned, group in ((requested[0], requested[1], attribute.Group.MOMENTS),
                                                        (requested[2], requested[3], attribute.Group.FILLED_SHAPE)):
                            scalar = attribute.compute_single_attribute(tree, getattr(attribute, signed), dtype=dtype)
                            np.testing.assert_array_equal(scalar, values[:, names[signed]])
                            np.testing.assert_array_equal(np.abs(scalar), values[:, names[unsigned]])
                            group_names, group_values = attribute.compute_attributes(tree, [group], dtype=dtype)
                            np.testing.assert_array_equal(scalar, group_values[:, group_names[signed]])
                            exported = tree.project_node_values_to_exported_higra(scalar, getattr(attribute, signed))
                            np.testing.assert_array_equal(exported[:image.size], 0)
                        for node in tree.alive_node_ids:
                            mask = np.zeros_like(image, dtype=bool)
                            mask.flat[np.fromiter(tree.node_support(node), dtype=int)] = True
                            for name, support in ((requested[0], mask), (requested[2], fill(mask, foreground_four))):
                                np.testing.assert_allclose(values[node, names[name]], angle(support), atol=1e-5, rtol=0)
    # The off-center hole tilts the support axis, but cannot tilt the filled rectangle.
    tree = mmcfilters.MorphologicalTreeFactory.create_max_tree(hole, 1.0)
    names, values = attribute.compute_attributes(tree, requested, dtype=np.float64)
    node = next(node for node in tree.alive_node_ids if tree.node_altitude(node) == 1)
    assert abs(values[node, names[requested[0]]]) > 1
    assert values[node, names[requested[2]]] == 0
    print("Signed support/filled orientations: raster oracles, rotations, holes, dtypes, groups and unit rows passed.")


def main():
    mmcfilters = load_native_module(pathlib.Path(sys.argv[1]).resolve())
    verify_signed_orientations(mmcfilters)
    attribute = mmcfilters.Attribute
    image = np.ones((5, 5), dtype=np.uint8)
    image[1, 1] = 0
    tree = mmcfilters.MorphologicalTreeFactory.create_max_tree(image, 1.0)
    expected_names = [
        "FILLED_AREA", "FILLED_CENTROID_ROW", "FILLED_CENTROID_COLUMN",
        "FILLED_LENGTH_MAJOR_AXIS", "FILLED_LENGTH_MINOR_AXIS",
        "FILLED_AXIS_ORIENTATION", "FILLED_ECCENTRICITY", "FILLED_INERTIA",
        "HOLE_AREA_FRACTION", "FILLED_CENTROID_DISPLACEMENT_NORMALIZED",
        "FILLED_COMPACTNESS", "FILLED_CIRCULARITY", "FILLED_AXIS_ORIENTATION_SIGNED",
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
