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
 */
template <std::floating_point Real>
[[nodiscard]] Real safeDivide(Real numerator, Real denominator, Real fallback = Real{0}) noexcept
{
    return std::abs(denominator) > std::numeric_limits<Real>::epsilon()
        ? numerator / denominator
        : fallback;
}

template <std::floating_point Real>
[[nodiscard]] Real clampNonNegative(Real value) noexcept
{
    return value > Real{0} ? value : Real{0};
}

template <std::floating_point Real>
[[nodiscard]] Real safeSqrt(Real value) noexcept
{
    return std::sqrt(clampNonNegative(value));
}

template <std::floating_point Real>
[[nodiscard]] Real clampUpper(Real value, Real upper) noexcept
{
    return std::min(value, upper);
}

} // namespace mmcfilters::attributes::numeric
