#pragma once

#include "NodeDistanceFieldProvider.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace mmcfilters::attributes::computers::detail::distance_transform {

/**
 * @brief Additive raw moments of one finite squared-distance field.
 *
 * The state deliberately supports insertion, removal, and component merging.
 * Exact reducers use only insertion; the adaptive DIFT observer additionally
 * removes invalidated labels and merges pairwise-incomparable components when
 * their propagation edge becomes active.
 */
class DistanceFieldMoments {
  public:
    /**
     * @brief Clears all samples.
     */
    void clear() noexcept {
        count_ = 0;
        sum_ = 0.0L;
        sumOfSquares_ = 0.0L;
        distanceSum_ = 0.0L;
    }

    /**
     * @brief Inserts one finite non-negative squared-distance sample.
     */
    template <bool ValidateInternalOperations = true> void add(SquaredDistance cost) {
        if constexpr (ValidateInternalOperations) {
            if (count_ == std::numeric_limits<std::uint64_t>::max()) {
                throw std::overflow_error("Distance-field sample count exceeds the 64-bit domain.");
            }
        }
        const long double value = static_cast<long double>(cost);
        ++count_;
        sum_ += value;
        sumOfSquares_ += value * value;
        distanceSum_ += std::sqrt(value);
    }

    /**
     * @brief Removes one previously inserted finite squared-distance sample.
     */
    template <bool ValidateInternalOperations = true> void remove(SquaredDistance cost) {
        if constexpr (ValidateInternalOperations) {
            if (count_ == 0) {
                throw std::logic_error("Distance-field moment removal requires a finite sample.");
            }
        }
        const long double value = static_cast<long double>(cost);
        --count_;
        sum_ -= value;
        sumOfSquares_ -= value * value;
        distanceSum_ -= std::sqrt(value);
        if (count_ == 0) {
            sum_ = 0.0L;
            sumOfSquares_ = 0.0L;
            distanceSum_ = 0.0L;
        }
    }

    /**
     * @brief Merges a disjoint component summary.
     */
    template <bool ValidateInternalOperations = true> void merge(const DistanceFieldMoments& other) {
        if constexpr (ValidateInternalOperations) {
            if (other.count_ > std::numeric_limits<std::uint64_t>::max() - count_) {
                throw std::overflow_error("Merged distance-field sample count exceeds the 64-bit domain.");
            }
        }
        count_ += other.count_;
        sum_ += other.sum_;
        sumOfSquares_ += other.sumOfSquares_;
        distanceSum_ += other.distanceSum_;
    }

    [[nodiscard]] std::uint64_t count() const noexcept { return count_; }
    [[nodiscard]] long double sum() const noexcept { return sum_; }
    [[nodiscard]] long double mean() const noexcept { return count_ == 0 ? 0.0L : sum_ / static_cast<long double>(count_); }
    [[nodiscard]] long double rms() const noexcept { return std::sqrt(std::max(0.0L, mean())); }
    [[nodiscard]] long double distanceSum() const noexcept { return distanceSum_; }
    [[nodiscard]] long double distanceMean() const noexcept {
        return count_ == 0 ? 0.0L : distanceSum_ / static_cast<long double>(count_);
    }

    /**
     * @brief Returns the population variance of the Euclidean distances.
     */
    [[nodiscard]] long double distancePopulationVariance() const noexcept {
        if (count_ == 0) {
            return 0.0L;
        }
        const long double average = distanceMean();
        const long double variance = mean() - average * average;
        return std::max(0.0L, variance);
    }

    /**
     * @brief Returns the population variance of the squared costs.
     */
    [[nodiscard]] long double populationVariance() const noexcept {
        if (count_ == 0) {
            return 0.0L;
        }
        const long double average = mean();
        const long double variance = sumOfSquares_ / static_cast<long double>(count_) - average * average;
        return std::max(0.0L, variance);
    }

  private:
    std::uint64_t count_ = 0;
    long double sum_ = 0.0L;
    long double sumOfSquares_ = 0.0L;
    long double distanceSum_ = 0.0L;
};

} // namespace mmcfilters::attributes::computers::detail::distance_transform
