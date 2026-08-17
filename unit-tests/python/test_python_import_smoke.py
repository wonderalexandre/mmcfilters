#!/usr/bin/env python3

import os
import sys

import numpy as np


def require(condition: bool, message: str):
    if not condition:
        raise RuntimeError(message)


def require_raises(fn, message: str):
    try:
        fn()
    except Exception:
        return
    raise RuntimeError(message)


def import_build_package():
    build_python = os.environ.get("MMCFILTERS_BUILD_PYTHON")
    if build_python:
        sys.path.insert(0, build_python)

    for name in list(sys.modules):
        if name == "mmcfilters" or name.startswith("mmcfilters."):
            sys.modules.pop(name, None)

    def handles_mmcfilters(finder):
        known_modules = {}
        known_modules.update(getattr(finder, "known_source_files", {}))
        known_modules.update(getattr(finder, "known_wheel_files", {}))
        return any(name == "mmcfilters" or name.startswith("mmcfilters.") for name in known_modules)

    sys.meta_path = [finder for finder in sys.meta_path if not handles_mmcfilters(finder)]

    import mmcfilters

    return mmcfilters


def main() -> int:
    mmcfilters = import_build_package()
    image = np.array([[1, 2], [3, 4]], dtype=np.uint8)
    tree = mmcfilters.MorphologicalTreeFactory.create_max_tree(image)

    require(tree.reconstruct_from_node_altitudes().shape == image.shape, "factory max-tree smoke")
    for public_name in (
        "NodePreservationMask",
        "NodePruningMask",
        "DirectReconstruction",
        "SubtractiveResidueModulation",
        "DirectAttributeFilter",
        "SubtractiveAttributeFilter",
        "SoftSubtractiveAttributeFilter",
        "IncompleteStabilityWindowPolicy",
        "compute_node_preservation_mask",
        "adjust_node_preservation_mask_by_altitude_stability",
        "adjust_node_preservation_mask_by_depth_stability",
        "to_node_pruning_mask",
        "to_node_preservation_mask",
    ):
        require(hasattr(mmcfilters, public_name), f"package must expose {public_name}")
    require(hasattr(tree, "reconstruct_from_node_contributions"), "valued tree must expose contribution reconstruction")
    require(hasattr(mmcfilters, "Attribute"), "package must expose Attribute")
    require(hasattr(mmcfilters, "NodeAttributeSamplingPolicy"), "package must expose NodeAttributeSamplingPolicy")
    require(hasattr(mmcfilters, "MissingNodeAttributeSamplePolicy"), "package must expose MissingNodeAttributeSamplePolicy")
    require(hasattr(mmcfilters.Attribute, "compute_sampled_node_attribute"), "Attribute must expose compute_sampled_node_attribute")
    require(
        not hasattr(mmcfilters, "AttributeOpeningPrimitivesFamily"),
        "package must not expose removed AttributeOpeningPrimitivesFamily",
    )
    require(hasattr(mmcfilters.Attribute, "CONTOUR_PIXELS"), "Attribute must expose CONTOUR_PIXELS")
    require(hasattr(mmcfilters.Attribute, "CONTOUR_SIDE_SOUTH"), "Attribute must expose CONTOUR_SIDE_SOUTH")
    require(hasattr(mmcfilters.Attribute.Group, "GRAY_LEVEL"), "Attribute.Group must expose GRAY_LEVEL")
    require(hasattr(mmcfilters.Attribute.Group, "SHAPE"), "Attribute.Group must expose SHAPE")
    require(hasattr(mmcfilters.Attribute.Group, "MOMENTS"), "Attribute.Group must expose MOMENTS")
    require(hasattr(mmcfilters.Attribute.Group, "BOUNDARY"), "Attribute.Group must expose BOUNDARY")
    for legacy_group in (
        "GEOMETRIC",
        "BOUNDING_BOX",
        "CENTRAL_MOMENTS",
        "HU_MOMENTS",
        "MOMENT_BASED",
        "TEXTURE",
        "BITQUADS",
        "CONTOUR",
    ):
        require(
            not hasattr(mmcfilters.Attribute.Group, legacy_group),
            f"Attribute.Group must not expose legacy {legacy_group}",
        )
    require(not hasattr(mmcfilters, "BitquadDeltas"), "package must not expose internal BitquadDeltas")

    area = mmcfilters.Attribute.compute_single_attribute(tree, mmcfilters.Attribute.AREA)
    require(area.shape[0] == tree.num_internal_node_slots, "single attribute smoke shape")

    topology_area = mmcfilters.Attribute.compute_single_topology_attribute(tree, mmcfilters.Attribute.AREA)
    require(topology_area.shape[0] == tree.num_internal_node_slots, "single topology attribute smoke shape")
    topology_names, topology_attrs = mmcfilters.Attribute.compute_topology_attributes(
        tree,
        [mmcfilters.Attribute.AREA, mmcfilters.Attribute.BOX_WIDTH],
    )
    require("AREA" in topology_names, "topology attribute facade must include AREA")
    require(topology_attrs.shape == (tree.num_internal_node_slots, len(topology_names)), "topology attributes smoke shape")
    require_raises(
        lambda: mmcfilters.Attribute.compute_single_topology_attribute(tree, mmcfilters.Attribute.MEAN_GRAY_LEVEL),
        "topology attribute facade must reject altitude-dependent attributes",
    )

    gray_names, gray_attrs = mmcfilters.Attribute.compute_attributes(tree, [mmcfilters.Attribute.Group.GRAY_LEVEL])
    require("VOLUME" in gray_names, "Attribute.Group.GRAY_LEVEL must include VOLUME")
    require("GRAY_LEVEL_VARIANCE" in gray_names, "Attribute.Group.GRAY_LEVEL must include GRAY_LEVEL_VARIANCE")
    require(gray_attrs.shape == (tree.num_internal_node_slots, len(gray_names)), "Attribute.Group.GRAY_LEVEL smoke shape")

    boundary_names, boundary_attrs = mmcfilters.Attribute.compute_topology_attributes(
        tree,
        [mmcfilters.Attribute.Group.BOUNDARY],
    )
    require("BITQUAD_AREA" in boundary_names, "Attribute.Group.BOUNDARY must include BITQUAD_AREA")
    require("CONTOUR_PIXELS" in boundary_names, "Attribute.Group.BOUNDARY must include CONTOUR_PIXELS")
    require("CONTOUR_SIDE_SOUTH" in boundary_names, "Attribute.Group.BOUNDARY must include CONTOUR_SIDE_SOUTH")
    require(
        boundary_attrs.shape == (tree.num_internal_node_slots, len(boundary_names)),
        "Attribute.Group.BOUNDARY topology smoke shape",
    )

    shape_names, shape_attrs = mmcfilters.Attribute.compute_attributes(
        tree,
        [mmcfilters.Attribute.Group.SHAPE],
    )
    require("AREA" in shape_names, "Attribute.Group.SHAPE must include AREA")
    require("MAX_DIST" in shape_names, "Attribute.Group.SHAPE must include MAX_DIST")
    require("CONTOUR_SIDE_SOUTH" in shape_names, "Attribute.Group.SHAPE must include CONTOUR_SIDE_SOUTH")
    require(
        shape_attrs.shape == (tree.num_internal_node_slots, len(shape_names)),
        "Attribute.Group.SHAPE smoke shape",
    )

    names, attrs = mmcfilters.Attribute.compute_attributes(tree, [mmcfilters.Attribute.ALL])
    require("AREA" in names, "Attribute.ALL must include AREA")
    require("CONTOUR_PIXELS" in names, "Attribute.ALL must include CONTOUR_PIXELS")
    require("CONTOUR_SIDE_SOUTH" in names, "Attribute.ALL must include CONTOUR_SIDE_SOUTH")
    require(attrs.shape == (tree.num_internal_node_slots, len(names)), "Attribute.ALL smoke shape")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
