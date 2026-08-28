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
    require(hasattr(mmcfilters.Attribute, "MAX_DIST_CENTER_ROW_EXACT"), "Attribute must expose MAX_DIST_CENTER_ROW_EXACT")
    require(hasattr(mmcfilters.Attribute, "MAX_DIST_CENTER_COLUMN"), "Attribute must expose MAX_DIST_CENTER_COLUMN")
    require(hasattr(mmcfilters.Attribute, "MAX_DIST_PLATEAU_AREA_EXACT"), "Attribute must expose MAX_DIST_PLATEAU_AREA_EXACT")
    require(
        hasattr(mmcfilters.Attribute, "MAX_DIST_PLATEAU_CENTROID_COLUMN"),
        "Attribute must expose MAX_DIST_PLATEAU_CENTROID_COLUMN",
    )
    require(hasattr(mmcfilters.Attribute, "DIST_SUM"), "Attribute must expose approximate DIST_SUM")
    require(hasattr(mmcfilters.Attribute, "DIST_SUM_EXACT"), "Attribute must expose exact DIST_SUM_EXACT")
    require(hasattr(mmcfilters.Attribute, "DIST_Q90"), "Attribute must expose approximate DIST_Q90")
    require(hasattr(mmcfilters.Attribute, "DIST_Q90_EXACT"), "Attribute must expose exact DIST_Q90_EXACT")
    require(
        hasattr(mmcfilters.Attribute, "DIST_WEIGHTED_ECCENTRICITY_EXACT"),
        "Attribute must expose the final exact distance-transform descriptor",
    )
    require(
        not hasattr(mmcfilters.Attribute, "MAX_DIST_SQUARED_APPROX"),
        "approximate MAX_DIST must not expose an APPROX-suffixed alias",
    )
    require(hasattr(mmcfilters.Attribute.Group, "GRAY_LEVEL"), "Attribute.Group must expose GRAY_LEVEL")
    require(hasattr(mmcfilters.Attribute.Group, "SHAPE"), "Attribute.Group must expose SHAPE")
    require(hasattr(mmcfilters.Attribute.Group, "MOMENTS"), "Attribute.Group must expose MOMENTS")
    require(hasattr(mmcfilters.Attribute.Group, "BOUNDARY"), "Attribute.Group must expose BOUNDARY")
    require(hasattr(mmcfilters.Attribute.Group, "DIST_TRANSF"), "Attribute.Group must expose DIST_TRANSF")
    require(hasattr(mmcfilters.Attribute.Group, "DIST_TRANSF_EXACT"), "Attribute.Group must expose DIST_TRANSF_EXACT")
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
    require("MAX_DIST_EXACT" in shape_names, "Attribute.Group.SHAPE must include MAX_DIST_EXACT")
    require(
        "MAX_DIST" in shape_names,
        "Attribute.Group.SHAPE must include MAX_DIST",
    )
    require("CONTOUR_SIDE_SOUTH" in shape_names, "Attribute.Group.SHAPE must include CONTOUR_SIDE_SOUTH")
    require(
        shape_attrs.shape == (tree.num_internal_node_slots, len(shape_names)),
        "Attribute.Group.SHAPE smoke shape",
    )

    distance_names, distance_attrs = mmcfilters.Attribute.compute_topology_attributes(
        tree,
        [mmcfilters.Attribute.Group.DIST_TRANSF],
    )
    ordered_distance_names = list(distance_names)
    require(len(distance_names) == 29, "Attribute.Group.DIST_TRANSF must expose all approximate distance-transform descriptors")
    require("MAX_DIST_EXACT" not in distance_names, "Attribute.Group.DIST_TRANSF must exclude MAX_DIST_EXACT")
    require(
        "MAX_DIST" in distance_names,
        "Attribute.Group.DIST_TRANSF must include MAX_DIST",
    )
    require("MAX_SQUARED_DIST_EXACT" not in distance_names, "Attribute.Group.DIST_TRANSF must exclude MAX_SQUARED_DIST_EXACT")
    require("MAX_SQUARED_DIST" in distance_names, "Attribute.Group.DIST_TRANSF must include MAX_SQUARED_DIST")
    require(
        all(not name.endswith("_EXACT") for name in ordered_distance_names),
        "DIST_TRANSF must contain the 29 unsuffixed approximate descriptors",
    )
    require(
        distance_attrs.shape == (tree.num_internal_node_slots, len(distance_names)),
        "Attribute.Group.DIST_TRANSF topology smoke shape",
    )
    string_distance_names, string_distance_attrs = mmcfilters.Attribute.compute_topology_attributes(
        tree,
        ["DIST_TRANSF"],
    )
    require(string_distance_names == distance_names, "concise string DIST_TRANSF request must preserve the group layout")
    require(
        np.array_equal(string_distance_attrs, distance_attrs, equal_nan=True),
        "concise string DIST_TRANSF request must preserve group values",
    )
    exact_distance_names, exact_distance_attrs = mmcfilters.Attribute.compute_topology_attributes(
        tree,
        [mmcfilters.Attribute.Group.DIST_TRANSF_EXACT],
    )
    require(len(exact_distance_names) == 29, "Attribute.Group.DIST_TRANSF_EXACT must expose all exact descriptors")
    require(all(name.endswith("_EXACT") for name in exact_distance_names), "DIST_TRANSF_EXACT descriptors must use _EXACT")
    require("MAX_DIST_EXACT" in exact_distance_names, "DIST_TRANSF_EXACT must include MAX_DIST_EXACT")
    require("MAX_DIST" not in exact_distance_names, "DIST_TRANSF_EXACT must exclude MAX_DIST")
    require(
        exact_distance_attrs.shape == (tree.num_internal_node_slots, len(exact_distance_names)),
        "Attribute.Group.DIST_TRANSF_EXACT topology smoke shape",
    )
    string_exact_names, string_exact_attrs = mmcfilters.Attribute.compute_topology_attributes(tree, ["DIST_TRANSF_EXACT"])
    require(string_exact_names == exact_distance_names, "concise DIST_TRANSF_EXACT request must preserve layout")
    require(np.array_equal(string_exact_attrs, exact_distance_attrs, equal_nan=True), "concise DIST_TRANSF_EXACT must preserve values")

    names, attrs = mmcfilters.Attribute.compute_attributes(tree, [mmcfilters.Attribute.ALL])
    require("AREA" in names, "Attribute.ALL must include AREA")
    require("DIST_SQUARED_MEAN_EXACT" in names, "Attribute.ALL must include DIST_SQUARED_MEAN_EXACT")
    require("DIST_SQUARED_MEAN" in names, "Attribute.ALL must include DIST_SQUARED_MEAN")
    require("MAX_DIST_CENTER_ROW_EXACT" in names, "Attribute.ALL must include MAX_DIST_CENTER_ROW_EXACT")
    require("MAX_DIST_CENTER_COLUMN" in names, "Attribute.ALL must include MAX_DIST_CENTER_COLUMN")
    require("MAX_DIST_PLATEAU_AREA_EXACT" in names, "Attribute.ALL must include MAX_DIST_PLATEAU_AREA_EXACT")
    require("DIST_ENTROPY" in names, "Attribute.ALL must include approximate DIST_ENTROPY")
    require("DIST_ENTROPY_EXACT" in names, "Attribute.ALL must include exact DIST_ENTROPY_EXACT")
    require(
        "DIST_WEIGHTED_ECCENTRICITY_EXACT" in names,
        "Attribute.ALL must include the final exact distance-transform descriptor",
    )
    require(
        "MAX_DIST_PLATEAU_CENTROID_COLUMN" in names,
        "Attribute.ALL must include MAX_DIST_PLATEAU_CENTROID_COLUMN",
    )
    require("CONTOUR_PIXELS" in names, "Attribute.ALL must include CONTOUR_PIXELS")
    require("CONTOUR_SIDE_SOUTH" in names, "Attribute.ALL must include CONTOUR_SIDE_SOUTH")
    require(attrs.shape == (tree.num_internal_node_slots, len(names)), "Attribute.ALL smoke shape")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
