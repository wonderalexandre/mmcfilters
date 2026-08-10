#pragma once

/**
 * @file SdrtTiePolicy.hpp
 * @brief Deterministic ordering policies shared by SDRT backends.
 */

namespace mmcfilters::sdrt {

/**
 * @defgroup sdrt Self-dual residual trees
 * @brief Direct construction and reference implementations of the SDRT.
 *
 * An SDRT alternates regional maxima and minima in one hierarchy. Production
 * construction is provided by the saturated and unrestricted residual-tree
 * builders from a synchronized min-tree/max-tree pair. Independent
 * constructions used for differential validation remain confined to the test
 * suite.
 */

/**
 * @brief Deterministic tie policies shared by SDRT construction backends.
 *
 * Tie handling is part of the public SDRT semantics because it controls the
 * chronological residual-node numbering and may change the hierarchy when
 * interacting extrema have the same area.
 *
 * @ingroup sdrt
 */
enum class SdrtTiePolicy {
    /// Sort by area, then Max before Min, then the smallest row-major pixel.
    MaxBeforeMinThenSpatial,

    /// Sort by area and the smallest row-major pixel, independently of polarity.
    ContrastInvariantSpatial
};

} // namespace mmcfilters::sdrt
