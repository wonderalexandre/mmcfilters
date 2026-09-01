#pragma once

#include "../utils/Common.hpp"

#include <cstdint>

namespace mmcfilters {

/** @brief Numerical backend selected for a public distance-transform product. */
enum class DistanceTransformAccuracy : std::uint8_t {
    /// Adaptive-A8 dynamic image foresting transform.
    Approximate,
    /// Exact squared Euclidean transform of the foreground A4 contour.
    Exact,
};

/** @brief One sparse squared-distance histogram level. */
struct DistanceTransformHistogramBin {
    /// Squared-distance level in non-negative pixel-squared units.
    std::int64_t squaredDistance = 0;
    /// Number of support samples at this level.
    std::uint64_t count = 0;
};

/** @brief One squared-distance threshold and its surviving support area. */
struct DistanceTransformSurvivalBin {
    /// Squared-distance threshold in non-negative pixel-squared units.
    std::int64_t squaredDistance = 0;
    /// Number of support samples whose squared distance is at least the threshold.
    std::uint64_t countAtOrAbove = 0;
};

} // namespace mmcfilters
