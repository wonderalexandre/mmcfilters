#pragma once

#include <algorithm>
#include <cmath>
#include <concepts>
#include <limits>

namespace mmcfilters::attributes::numeric {

/**
 * @brief Numerically safe division used by attribute kernels.
 *
 * Attribute descriptors often become undefined on degenerate supports. The
 * library policy is to materialize the documented finite fallback instead of
 * leaking NaN or infinity into public buffers.
 *
 * @param numerator Number represented by `numerator`.
 * @param denominator Division denominator.
 * @param fallback Value returned when division is not numerically valid.
 * @return Numerically safe division used by attribute kernels.
 */
template <std::floating_point Real> [[nodiscard]] Real safeDivide(Real numerator, Real denominator, Real fallback = Real{0}) noexcept {
    return std::abs(denominator) > std::numeric_limits<Real>::epsilon() ? numerator / denominator : fallback;
}

/**
 * @brief Clamps non negative.
 *
 * @param value Value used by the operation.
 * @return `value` clamped to the non-negative range.
 */
template <std::floating_point Real> [[nodiscard]] Real clampNonNegative(Real value) noexcept { return value > Real{0} ? value : Real{0}; }

/**
 * @brief Computes sqrt.
 *
 * @param value Value used by the operation.
 * @return Computed sqrt.
 */
template <std::floating_point Real> [[nodiscard]] Real safeSqrt(Real value) noexcept { return std::sqrt(clampNonNegative(value)); }

/**
 * @brief Clamps upper.
 *
 * @param value Value used by the operation.
 * @param upper Upper admissible priority.
 * @return Minimum of `value` and `upper`.
 */
template <std::floating_point Real> [[nodiscard]] Real clampUpper(Real value, Real upper) noexcept { return std::min(value, upper); }

} // namespace mmcfilters::attributes::numeric
