#pragma once

#include "../../../../trees/MorphologicalTree.hpp"
#include "../../../../utils/Common.hpp"

#include <algorithm>
#include <concepts>
#include <cstdint>

namespace mmcfilters::attributes::computers::detail::distance_transform {

/**
 * @brief Integer sample domain shared by exact and approximate providers.
 */
using SquaredDistance = std::int64_t;

/**
 * @brief Maximum squared distance and its deterministic row-major location.
 */
struct DistanceFieldExtremum {
    SquaredDistance squaredDistance = 0;
    PixelId pixel = InvalidPixel;
};

/**
 * @brief Maximum squared distance and geometry of its complete level-set plateau.
 */
struct DistanceFieldMaximumPlateau {
    SquaredDistance squaredDistance = 0;
    PixelId pixel = InvalidPixel;
    std::uint64_t count = 0;
    std::uint64_t rowSum = 0;
    std::uint64_t columnSum = 0;

    [[nodiscard]] long double centroidRow() const noexcept { return count == 0 ? 0.0L : static_cast<long double>(rowSum) / static_cast<long double>(count); }

    [[nodiscard]] long double centroidColumn() const noexcept {
        return count == 0 ? 0.0L : static_cast<long double>(columnSum) / static_cast<long double>(count);
    }
};

/**
 * @brief Tests whether a sample supersedes the current maximum and tie-break.
 */
[[nodiscard]] inline constexpr bool prefersDistanceFieldSample(SquaredDistance squaredDistance, PixelId pixel, const DistanceFieldExtremum& current) noexcept {
    return current.pixel == InvalidPixel || squaredDistance > current.squaredDistance || (squaredDistance == current.squaredDistance && pixel < current.pixel);
}

/**
 * @brief Inserts a sample into a maximum with smallest-row-major tie-breaking.
 */
inline constexpr void updateDistanceFieldExtremum(DistanceFieldExtremum& extremum, PixelId pixel, SquaredDistance squaredDistance) noexcept {
    if (prefersDistanceFieldSample(squaredDistance, pixel, extremum)) {
        extremum = DistanceFieldExtremum{squaredDistance, pixel};
    }
}

/**
 * @brief Inserts one sample into the maximum-distance plateau summary.
 */
inline constexpr void updateDistanceFieldMaximumPlateau(DistanceFieldMaximumPlateau& plateau, PixelId pixel, SquaredDistance squaredDistance,
                                                        int numColumns) noexcept {
    const std::uint64_t row = static_cast<std::uint64_t>(pixel / numColumns);
    const std::uint64_t column = static_cast<std::uint64_t>(pixel % numColumns);
    if (plateau.pixel == InvalidPixel || squaredDistance > plateau.squaredDistance) {
        plateau = DistanceFieldMaximumPlateau{squaredDistance, pixel, 1, row, column};
    } else if (squaredDistance == plateau.squaredDistance) {
        plateau.pixel = std::min(plateau.pixel, pixel);
        ++plateau.count;
        plateau.rowSum += row;
        plateau.columnSum += column;
    }
}

/**
 * @brief Merges a disjoint per-root plateau into a node-level maximum plateau.
 */
inline constexpr void mergeDistanceFieldMaximumPlateau(DistanceFieldMaximumPlateau& plateau, const DistanceFieldMaximumPlateau& other) noexcept {
    if (other.pixel == InvalidPixel || other.count == 0) {
        return;
    }
    if (plateau.pixel == InvalidPixel || other.squaredDistance > plateau.squaredDistance) {
        plateau = other;
    } else if (other.squaredDistance == plateau.squaredDistance) {
        plateau.pixel = std::min(plateau.pixel, other.pixel);
        plateau.count += other.count;
        plateau.rowSum += other.rowSum;
        plateau.columnSum += other.columnSum;
    }
}

/**
 * @brief Numerical contract advertised by a node distance-field provider.
 */
enum class DistanceFieldAccuracy : std::uint8_t {
    Exact,
    Approximate,
};

/**
 * @brief Protocol for a synchronous reduction over node distance samples.
 *
 * Reducers own their accumulated state and must not retain borrowed provider
 * storage beyond a callback.
 */
template <class Reducer>
concept NodeDistanceTransformReducer = requires(Reducer& reducer, NodeId node, PixelId pixel, SquaredDistance squaredDistance) {
    { reducer.beginNode(node) } -> std::same_as<void>;
    { reducer.consumeSample(pixel, squaredDistance) } -> std::same_as<void>;
    { reducer.endNode(node) } -> std::same_as<void>;
};

/**
 * @brief Tests whether a provider can feed one reducer from a tree.
 *
 * A provider owns scheduling and numerical field computation. The reducer owns
 * only the attribute-specific projection, which is the separation needed to
 * compare exact and approximate backends without duplicating `MAX_DIST_EXACT` logic.
 */
template <class Provider, class Reducer>
concept NodeDistanceFieldProviderFor = NodeDistanceTransformReducer<Reducer> && requires(const MorphologicalTree& tree, Reducer& reducer) {
    { Provider::accuracy } -> std::convertible_to<DistanceFieldAccuracy>;
    { Provider::reduce(tree, reducer) } -> std::same_as<void>;
};

} // namespace mmcfilters::attributes::computers::detail::distance_transform
