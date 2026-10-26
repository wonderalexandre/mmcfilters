#!/usr/bin/env python3

import pathlib
import sys

import numpy as np

from test_python_attribute_numeric_validation import load_native_module


def check_batch(attribute, tree, request, output_space, rows, dtype):
    layout, values = attribute.compute(tree, request, output_space, dtype=dtype)
    assert values.shape == (rows, len(layout))
    assert values.dtype == np.dtype(dtype)
    assert values.flags.c_contiguous
    assert list(layout.values()) == list(range(len(layout)))
    assert list(layout) == sorted(layout, key=lambda name: int(getattr(attribute, name)))
    for name, column in layout.items():
        scalar = attribute.compute(tree, name, output_space, dtype=dtype)
        assert scalar.shape == (rows,)
        assert scalar.dtype == np.dtype(dtype)
        assert scalar.flags.c_contiguous
        np.testing.assert_array_equal(values[:, column], scalar)
    return layout, values


def check_requests(mmcfilters, tree, output_space, rows):
    attribute = mmcfilters.Attribute
    for dtype in (np.float32, np.float64):
        for name in ("AREA", "BOUNDING_BOX_WIDTH", "VOLUME", "GRAY_LEVEL_HEIGHT"):
            named = attribute.compute(tree, name, output_space, dtype=dtype)
            typed = attribute.compute(tree, getattr(attribute, name), output_space, dtype=dtype)
            assert named.shape == (rows,)
            assert named.dtype == np.dtype(dtype)
            assert named.flags.c_contiguous
            np.testing.assert_array_equal(named, typed)

        # Sequences keep a column dimension. Repeated attributes and
        # overlapping groups use one column per scalar in canonical order.
        requests = [
            ["AREA"],
            ["BOUNDING_BOX_WIDTH", attribute.AREA, "AREA"],
            [attribute.Group.TREE_TOPOLOGY],
            ["GRAY_LEVEL"],
            ["AREA", attribute.Group.GRAY_LEVEL, "VOLUME"],
            ["SHAPE"],
            [],
        ]
        for request in requests:
            layout, values = check_batch(attribute, tree, request, output_space, rows, dtype)
            tuple_layout, tuple_values = attribute.compute(tree, tuple(request), output_space, dtype=dtype)
            assert list(layout.items()) == list(tuple_layout.items())
            np.testing.assert_array_equal(values, tuple_values)
        assert layout == {} and values.shape == (rows, 0)

        layout, _ = attribute.compute(tree, requests[1], output_space, dtype=dtype)
        assert list(layout) == ["AREA", "BOUNDING_BOX_WIDTH"]
        layout, _ = attribute.compute(tree, requests[4], output_space, dtype=dtype)
        assert list(layout) == ["AREA", "VOLUME", "RELATIVE_VOLUME", "GRAY_LEVEL_HEIGHT",
                                "MEAN_GRAY_LEVEL", "GRAY_LEVEL_VARIANCE"]

        for name in ("TREE_TOPOLOGY", "GRAY_LEVEL", "SHAPE"):
            layout, values = check_batch(attribute, tree, name, output_space, rows, dtype)
            for request in (getattr(attribute.Group, name), [name]):
                other_layout, other_values = attribute.compute(tree, request, output_space, dtype=dtype)
                assert list(layout.items()) == list(other_layout.items())
                np.testing.assert_array_equal(values, other_values)

    check_batch(attribute, tree, (name for name in ("AREA", "VOLUME")),
                output_space, rows, np.float32)
    assert attribute.compute(tree=tree, attributes="AREA").dtype == np.float32
    assert attribute.compute(tree, "AREA", dtype="float64").dtype == np.float64

    # Check topology and altitude results against the actual pixel support,
    # independently of the attribute pipeline, including after node removal.
    image = tree.reconstruct_from_node_altitudes()
    area = attribute.compute(tree, "AREA")
    volume = attribute.compute(tree, "VOLUME", dtype=np.float64)
    width = attribute.compute(tree, "BOUNDING_BOX_WIDTH")
    for node in tree.alive_node_ids:
        pixels = np.array(list(tree.node_support(node)), dtype=np.intp)
        assert area[node] == len(pixels)
        assert volume[node] == image.ravel()[pixels].sum()
        columns = pixels % image.shape[1]
        assert width[node] == columns.max() - columns.min() + 1


def main():
    if len(sys.argv) != 2:
        raise SystemExit("usage: test_python_attribute_compute.py <build-dir>")
    mmcfilters = load_native_module(pathlib.Path(sys.argv[1]).resolve())
    attribute = mmcfilters.Attribute
    factory = mmcfilters.MorphologicalTreeFactory
    image = np.array([[0, 0, 0, 0], [0, 3, 5, 0], [0, 3, 8, 0]], dtype=np.uint8)
    internal = mmcfilters.NodeIdSpace.MORPHOLOGICAL_TREE

    for old_method in ("compute_single_attribute", "compute_attributes",
                       "compute_single_topology_attribute", "compute_topology_attributes"):
        assert not hasattr(attribute, old_method), old_method

    for build in (factory.create_max_tree, factory.create_min_tree):
        tree = build(image, radius=1.5)
        check_requests(mmcfilters, tree, internal, tree.num_internal_node_slots)
        assert attribute.compute(tree, "AREA")[tree.root] == image.size

    # Projection to imported Higra ids preserves rows and scalar unit values.
    tree = factory.create_max_tree(image, radius=1.5)
    parent, altitudes = tree.export_higra_hierarchy()
    imported = factory.create_from_higra_parent(
        parent, np.asarray(altitudes, dtype=np.uint8), *image.shape,
        mmcfilters.MorphologicalTreeKind.MAX_TREE, 1.5,
    )
    check_requests(mmcfilters, imported, mmcfilters.NodeIdSpace.HIGRA, len(parent))
    np.testing.assert_array_equal(
        attribute.compute(imported, "AREA", output_space=mmcfilters.NodeIdSpace.HIGRA)[:image.size],
        np.ones(image.size, dtype=np.float32),
    )

    # Edits keep the dense slot domain, including dead slots. Returned arrays
    # own their values and remain usable after the source tree changes.
    saved = attribute.compute(tree, "AREA")
    saved_copy = saved.copy()
    slots = tree.num_internal_node_slots
    leaf = next(node for node in tree.leaves if node != tree.root)
    tree.merge_node_into_parent(leaf)
    assert tree.num_nodes < slots
    assert tree.num_internal_node_slots == slots
    check_requests(mmcfilters, tree, internal, slots)
    np.testing.assert_array_equal(saved, saved_copy)

    # Boundary attributes on trees of shapes need the valued-tree context
    # even when their registry requirements select the topology path.
    shapes = factory.create_tree_of_shapes(image)
    check_batch(attribute, shapes, attribute.Group.BOUNDARY, internal,
                shapes.num_internal_node_slots, np.float64)

    for request in (None, 1, object(), [object()]):
        with np.testing.assert_raises(TypeError):
            attribute.compute(tree, request)
    for request in ("area", "NOT_AN_ATTRIBUTE", ["AREA", "NOT_A_GROUP"]):
        with np.testing.assert_raises(ValueError):
            attribute.compute(tree, request)
    for request in ("AREA", "SHAPE", ["AREA", "VOLUME"]):
        with np.testing.assert_raises(ValueError):
            attribute.compute(tree, request, dtype=np.int32)
    with np.testing.assert_raises(TypeError):
        attribute.compute(None, "AREA")

    print("Unified attribute compute tests passed")


if __name__ == "__main__":
    main()
