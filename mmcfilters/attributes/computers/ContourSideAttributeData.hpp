#pragma once

namespace mmcfilters::attributes::computers {

/**
 * @brief Projected local contour-side counters.
 *
 * @details
 * `contourPixels` counts visible anchor pixels with at least one exposed side.
 * `exposedSides` is the 4-neighbour discrete perimeter contribution and equals
 * `north + west + east + south`.
 */
struct ContourSideCounts {
    /// Number of pixels with at least one exposed side.
    int contourPixels = 0;

    /// Total number of exposed 4-neighbour sides.
    int exposedSides = 0;

    /// Number of exposed north-facing sides.
    int north = 0;

    /// Number of exposed west-facing sides.
    int west = 0;

    /// Number of exposed east-facing sides.
    int east = 0;

    /// Number of exposed south-facing sides.
    int south = 0;
};

} // namespace mmcfilters::attributes::computers
