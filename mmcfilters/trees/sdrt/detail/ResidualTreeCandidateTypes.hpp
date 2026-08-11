#pragma once

/**
 * @file ResidualTreeCandidateTypes.hpp
 * @brief Neutral candidate metadata shared by residual-tree construction components.
 */

#include "../../../utils/Common.hpp"

namespace mmcfilters::sdrt::detail {

/** @brief Derived properties of one current residual-tree candidate. */
struct ResidualTreeCandidateDescriptor {
    int area = 0;                              ///< Current flat-zone area.
    NodeId stableSpatialKey = InvalidNode;     ///< Smallest row-major pixel in the flat zone.
    bool containsExteriorSeed = false;         ///< Whether the flat zone contains the configured exterior seed.
};

} // namespace mmcfilters::sdrt::detail
