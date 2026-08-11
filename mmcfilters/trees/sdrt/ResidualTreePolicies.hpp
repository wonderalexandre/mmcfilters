#pragma once

/**
 * @file ResidualTreePolicies.hpp
 * @brief Public construction policies and option objects for residual trees.
 */

#include "SdrtTiePolicy.hpp"

namespace mmcfilters::sdrt {

/** @brief Lowest-common-ancestor strategy used by saturated certification. */
enum class SaturatedMinMaxLcaPolicy {
    ParentClimb,   ///< Climb current parent paths.
    BlockedSnapshot, ///< Use periodically rebuilt Euler/RMQ snapshots.
    LinkCut        ///< Maintain the rooted forest with link-cut operations.
};

/** @brief Exact complement-connectivity fallback after an inconclusive certificate. */
enum class SaturatedMinMaxFallbackPolicy {
    SingleSourceDepthFirst, ///< Traverse the complement from one boundary pixel.
    BoundaryMultiSource     ///< Grow and merge all boundary components.
};

/** @brief Options that affect saturated residual-tree construction. */
struct SaturatedResidualTreeOptions {
    SdrtTiePolicy tiePolicy = SdrtTiePolicy::ContrastInvariantSpatial; ///< Equal-area event ordering.
    SaturatedMinMaxLcaPolicy lcaPolicy = SaturatedMinMaxLcaPolicy::ParentClimb; ///< Dynamic LCA strategy.
    SaturatedMinMaxFallbackPolicy fallbackPolicy = SaturatedMinMaxFallbackPolicy::BoundaryMultiSource; ///< Exact fallback strategy.
};

/** @brief Options that affect unrestricted residual-tree construction. */
struct UnrestrictedResidualTreeOptions {
    SdrtTiePolicy tiePolicy = SdrtTiePolicy::ContrastInvariantSpatial; ///< Equal-area event ordering.
};

} // namespace mmcfilters::sdrt
