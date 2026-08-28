#pragma once

#include "NodeDistanceFieldProvider.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numbers>
#include <stdexcept>

namespace mmcfilters::attributes::computers::detail::distance_transform {

/**
 * @brief Additive spatial moments weighted by Euclidean contour distance.
 *
 * Coordinates use zero-based image rows and columns. The weight of a sample
 * is `sqrt(squaredDistance)`. Unweighted first moments are retained solely to
 * define the centroid of all-zero fields without emitting NaNs.
 */
class DistanceWeightedSpatialMoments {
  public:
    void clear() noexcept { *this = {}; }

    template <bool ValidateInternalOperations = true> void add(PixelId pixel, SquaredDistance squaredDistance, int numColumns) {
        update<ValidateInternalOperations>(pixel, squaredDistance, numColumns, 1.0L);
    }

    template <bool ValidateInternalOperations = true> void remove(PixelId pixel, SquaredDistance squaredDistance, int numColumns) {
        if constexpr (ValidateInternalOperations) {
            if (count_ == 0) {
                throw std::logic_error("Distance-weighted spatial removal requires an existing sample.");
            }
        }
        update<ValidateInternalOperations>(pixel, squaredDistance, numColumns, -1.0L);
        if (count_ == 0) {
            clear();
        }
    }

    template <bool ValidateInternalOperations = true> void merge(const DistanceWeightedSpatialMoments& other) {
        if constexpr (ValidateInternalOperations) {
            if (other.count_ > std::numeric_limits<std::uint64_t>::max() - count_) {
                throw std::overflow_error("Merged distance-weighted sample count exceeds the 64-bit domain.");
            }
        }
        count_ += other.count_;
        rowSum_ += other.rowSum_;
        columnSum_ += other.columnSum_;
        weightSum_ += other.weightSum_;
        weightedRowSum_ += other.weightedRowSum_;
        weightedColumnSum_ += other.weightedColumnSum_;
        weightedRowSquaredSum_ += other.weightedRowSquaredSum_;
        weightedColumnSquaredSum_ += other.weightedColumnSquaredSum_;
        weightedRowColumnSum_ += other.weightedRowColumnSum_;
    }

    [[nodiscard]] std::uint64_t count() const noexcept { return count_; }
    [[nodiscard]] long double weightSum() const noexcept { return weightSum_; }

    [[nodiscard]] long double centroidRow() const noexcept {
        if (weightSum_ > 0.0L) {
            return weightedRowSum_ / weightSum_;
        }
        return count_ == 0 ? 0.0L : rowSum_ / static_cast<long double>(count_);
    }

    [[nodiscard]] long double centroidColumn() const noexcept {
        if (weightSum_ > 0.0L) {
            return weightedColumnSum_ / weightSum_;
        }
        return count_ == 0 ? 0.0L : columnSum_ / static_cast<long double>(count_);
    }

    /**
     * @brief Returns the unnormalized weighted column central moment.
     */
    [[nodiscard]] long double centralMoment20() const noexcept {
        return weightSum_ <= 0.0L
                   ? 0.0L
                   : std::max(0.0L, weightedColumnSquaredSum_ - weightedColumnSum_ * weightedColumnSum_ / weightSum_);
    }

    /**
     * @brief Returns the unnormalized weighted row central moment.
     */
    [[nodiscard]] long double centralMoment02() const noexcept {
        return weightSum_ <= 0.0L ? 0.0L : std::max(0.0L, weightedRowSquaredSum_ - weightedRowSum_ * weightedRowSum_ / weightSum_);
    }

    /**
     * @brief Returns the unnormalized weighted row-column central moment.
     */
    [[nodiscard]] long double centralMoment11() const noexcept {
        return weightSum_ <= 0.0L ? 0.0L : weightedRowColumnSum_ - weightedRowSum_ * weightedColumnSum_ / weightSum_;
    }

    /**
     * @brief Returns principal-axis orientation in degrees.
     */
    [[nodiscard]] long double axisOrientationDegrees() const noexcept {
        const long double mu20 = centralMoment20();
        const long double mu02 = centralMoment02();
        const long double mu11 = centralMoment11();
        const long double scale = std::max({1.0L, std::abs(mu20), std::abs(mu02)});
        if (std::hypot(mu20 - mu02, 2.0L * mu11) <= 64.0L * std::numeric_limits<long double>::epsilon() * scale) {
            return 0.0L;
        }
        const long double radians = 0.5L * std::atan2(2.0L * mu11, mu20 - mu02);
        const long double degrees = radians * 180.0L / std::numbers::pi_v<long double>;
        return std::fmod(std::abs(degrees), 360.0L);
    }

    /**
     * @brief Returns the major/minor second-moment eigenvalue ratio.
     */
    [[nodiscard]] long double eccentricity() const noexcept {
        constexpr long double MaximumFiniteEccentricity = 1.0e6L;
        const long double mu20 = centralMoment20();
        const long double mu02 = centralMoment02();
        const long double mu11 = centralMoment11();
        const long double trace = mu20 + mu02;
        const long double discriminant = std::sqrt(std::max(0.0L, (mu20 - mu02) * (mu20 - mu02) + 4.0L * mu11 * mu11));
        const long double major = 0.5L * (trace + discriminant);
        const long double minor = std::max(0.0L, 0.5L * (trace - discriminant));
        if (major <= std::numeric_limits<long double>::epsilon()) {
            return 1.0L;
        }
        if (minor <= std::numeric_limits<long double>::epsilon()) {
            return MaximumFiniteEccentricity;
        }
        return std::min(MaximumFiniteEccentricity, major / minor);
    }

  private:
    template <bool ValidateInternalOperations>
    void update(PixelId pixel, SquaredDistance squaredDistance, int numColumns, long double sign) {
        if constexpr (ValidateInternalOperations) {
            if (pixel < 0 || squaredDistance < 0 || numColumns <= 0) {
                throw std::invalid_argument("Distance-weighted spatial moments require a valid sample and 2D domain.");
            }
            if (sign > 0.0L && count_ == std::numeric_limits<std::uint64_t>::max()) {
                throw std::overflow_error("Distance-weighted sample count exceeds the 64-bit domain.");
            }
        }
        const long double row = static_cast<long double>(pixel / numColumns);
        const long double column = static_cast<long double>(pixel % numColumns);
        const long double weight = std::sqrt(static_cast<long double>(squaredDistance));
        if (sign > 0.0L) {
            ++count_;
        } else {
            --count_;
        }
        rowSum_ += sign * row;
        columnSum_ += sign * column;
        weightSum_ += sign * weight;
        weightedRowSum_ += sign * weight * row;
        weightedColumnSum_ += sign * weight * column;
        weightedRowSquaredSum_ += sign * weight * row * row;
        weightedColumnSquaredSum_ += sign * weight * column * column;
        weightedRowColumnSum_ += sign * weight * row * column;
    }

    std::uint64_t count_ = 0;
    long double rowSum_ = 0.0L;
    long double columnSum_ = 0.0L;
    long double weightSum_ = 0.0L;
    long double weightedRowSum_ = 0.0L;
    long double weightedColumnSum_ = 0.0L;
    long double weightedRowSquaredSum_ = 0.0L;
    long double weightedColumnSquaredSum_ = 0.0L;
    long double weightedRowColumnSum_ = 0.0L;
};

} // namespace mmcfilters::attributes::computers::detail::distance_transform
