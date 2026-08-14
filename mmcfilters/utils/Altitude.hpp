#pragma once

/**
 * @file Altitude.hpp
 * @brief Shared altitude type contracts and buffer aliases.
 */

#include <concepts>
#include <cstdint>
#include <span>
#include <type_traits>
#include <vector>

namespace mmcfilters {

namespace detail {

/// True when `T` is an integral altitude type with safe 64-bit differences.
template <class T>
inline constexpr bool SupportedIntegralAltitude = std::is_integral_v<std::remove_cv_t<T>> && sizeof(std::remove_cv_t<T>) < sizeof(std::int64_t);

/// True when `T` is a floating-point altitude type.
template <class T> inline constexpr bool SupportedFloatingAltitude = std::is_floating_point_v<std::remove_cv_t<T>>;

} // namespace detail

/**
 * @brief Value types accepted by generic altitude helpers.
 *
 * The concept is intentionally conservative: floating-point values are accepted,
 * and integral values must fit safely inside the signed 64-bit difference type.
 * `bool` is rejected because it has no meaningful altitude interval semantics.
 *
 * Public valued-tree and attribute APIs use this concept to keep altitude
 * arithmetic independent of the original image pixel type while preserving the
 * monotone ordering required by max-trees and min-trees.
 */
template <class T>
concept AltitudeValue =
    std::totally_ordered<T> && !std::is_same_v<std::remove_cv_t<T>, bool> && (detail::SupportedFloatingAltitude<T> || detail::SupportedIntegralAltitude<T>);

/**
 * @brief Read-only contiguous view over node altitudes.
 *
 * The span is indexed by dense internal `NodeId`. Callers must provide at least
 * one value for every internal node slot in the associated topology.
 */
template <AltitudeValue T> using NodeAltitudeSpan = std::span<const T>;

namespace detail {

/**
 * @brief Selects a safe difference type for an altitude type.
 *
 * @tparam T Altitude type whose safe arithmetic difference type is selected.
 */
template <AltitudeValue T> struct AltitudeDifferenceSelector {
    /** @brief Defines the `type` alias used by the component. */
    using type = std::conditional_t<std::is_integral_v<T>, std::int64_t, T>;
};

} // namespace detail

/**
 * @brief Arithmetic result type for altitude differences.
 *
 * Integral altitudes use a signed 64-bit difference type to avoid narrowing
 * when comparing typed trees whose altitude domain is wider than `uint8_t`.
 * Floating-point altitudes keep their own precision so callers can express
 * real-valued differences such as `0.05f`.
 */
template <AltitudeValue T> using AltitudeDifference = typename detail::AltitudeDifferenceSelector<T>::type;

/**
 * @brief Owning dense altitude buffer indexed by internal node id.
 *
 * `ValuedMorphologicalTree<T>` owns one of these buffers. Read-only kernels
 * may instead accept a `NodeAltitudeSpan<T>` or `ValuedMorphologicalTreeView<T>` when the
 * topology and altitude storage have separate owners.
 */
template <AltitudeValue T> using NodeAltitudeBuffer = std::vector<T>;

} // namespace mmcfilters
