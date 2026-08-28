#pragma once

namespace mmcfilters::attributes::computers {

/**
 * @brief Stable family id declared by each concrete attribute computer.
 *
 * @details
 * The scheduler uses this id to group scalar descriptors by the stateless
 * computer that materializes them. Public callers should still request
 * attributes, not families.
 */
enum class AttributeComputerFamily {
    Area,
    Volume,
    GrayLevelStats,
    MaxDist,
    MaxDistExact,
    BoundingBox,
    TreeTopology,
    CentralMoments,
    HuMoments,
    MomentDerived,
    Bitquad,
    ContourSide,
    Unsupported
};

} // namespace mmcfilters::attributes::computers
