#pragma once

/**
 * @file ResidualTreeCandidateTypes.hpp
 * @brief Neutral candidate metadata shared by residual-tree construction components.
 */

#include "../../../utils/Common.hpp"

#include <cstddef>

namespace mmcfilters::sdrt::detail {

/** @brief Derived properties of one current residual-tree candidate. */
struct ResidualTreeCandidateDescriptor {
    std::size_t supportCardinality = 0;         ///< Cardinality of the current flat-zone support.
    PixelId spatialMinimum = InvalidPixel;      ///< Least support pixel in the configured spatial order.
    bool containsInfinityPixel = false;        ///< Whether the flat zone contains the declared infinity pixel.
};

} // namespace mmcfilters::sdrt::detail
