#pragma once

#include "NodeDistanceFieldProvider.hpp"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>

namespace mmcfilters::attributes::computers::detail::distance_transform {

/**
 * @brief One sorted squared-distance level and its sample multiplicity.
 */
struct DistanceHistogramBin {
    SquaredDistance squaredDistance = 0;
    std::uint64_t count = 0;
};

/**
 * @brief Sparse deterministic histogram of one finite squared-distance field.
 *
 * The histogram supports the same insertion, removal, and disjoint-component
 * merge protocol as the dynamic approximate DIFT observer. Ordered keys give
 * deterministic modes, quantiles, and public profile materialization.
 */
class DistanceFieldHistogram {
  public:
    using Storage = std::map<SquaredDistance, std::uint64_t>;

    void clear() noexcept {
        bins_.clear();
        count_ = 0;
    }

    template <bool ValidateInternalOperations = true> void add(SquaredDistance squaredDistance) {
        if constexpr (ValidateInternalOperations) {
            if (squaredDistance < 0) {
                throw std::invalid_argument("Distance-field histogram requires a non-negative squared distance.");
            }
            if (count_ == std::numeric_limits<std::uint64_t>::max()) {
                throw std::overflow_error("Distance-field histogram sample count exceeds the 64-bit domain.");
            }
        }
        std::uint64_t& binCount = bins_[squaredDistance];
        if constexpr (ValidateInternalOperations) {
            if (binCount == std::numeric_limits<std::uint64_t>::max()) {
                throw std::overflow_error("Distance-field histogram bin count exceeds the 64-bit domain.");
            }
        }
        ++binCount;
        ++count_;
    }

    template <bool ValidateInternalOperations = true> void remove(SquaredDistance squaredDistance) {
        const auto it = bins_.find(squaredDistance);
        if constexpr (ValidateInternalOperations) {
            if (count_ == 0 || it == bins_.end() || it->second == 0) {
                throw std::logic_error("Distance-field histogram removal requires an existing sample.");
            }
        }
        --it->second;
        --count_;
        if (it->second == 0) {
            bins_.erase(it);
        }
    }

    template <bool ValidateInternalOperations = true> void merge(DistanceFieldHistogram& other) {
        if constexpr (ValidateInternalOperations) {
            if (other.count_ > std::numeric_limits<std::uint64_t>::max() - count_) {
                throw std::overflow_error("Merged distance-field histogram sample count exceeds the 64-bit domain.");
            }
        }
        if (bins_.size() < other.bins_.size()) {
            bins_.swap(other.bins_);
            std::swap(count_, other.count_);
        }
        for (const auto& [squaredDistance, otherCount] : other.bins_) {
            std::uint64_t& binCount = bins_[squaredDistance];
            if constexpr (ValidateInternalOperations) {
                if (otherCount > std::numeric_limits<std::uint64_t>::max() - binCount) {
                    throw std::overflow_error("Merged distance-field histogram bin count exceeds the 64-bit domain.");
                }
            }
            binCount += otherCount;
        }
        count_ += other.count_;
        other.clear();
    }

    [[nodiscard]] std::uint64_t count() const noexcept { return count_; }
    [[nodiscard]] std::size_t levelCount() const noexcept { return bins_.size(); }
    [[nodiscard]] const Storage& bins() const noexcept { return bins_; }

    [[nodiscard]] std::uint64_t positiveArea() const noexcept {
        const auto zero = bins_.find(SquaredDistance{0});
        return zero == bins_.end() ? count_ : count_ - zero->second;
    }

    /**
     * @brief Returns the lower quantile as a Euclidean distance in pixels.
     */
    [[nodiscard]] long double quantile(long double probability) const {
        if (bins_.empty() || count_ == 0) {
            return 0.0L;
        }
        if (!(probability >= 0.0L && probability <= 1.0L)) {
            throw std::invalid_argument("Distance-field quantile probability must belong to [0, 1].");
        }
        const long double scaledRank = probability * static_cast<long double>(count_);
        const std::uint64_t target = probability <= 0.0L
                                         ? std::uint64_t{1}
                                         : static_cast<std::uint64_t>(std::ceil(std::max(1.0L, scaledRank)));
        std::uint64_t cumulative = 0;
        for (const auto& [squaredDistance, binCount] : bins_) {
            cumulative += binCount;
            if (cumulative >= target) {
                return std::sqrt(static_cast<long double>(squaredDistance));
            }
        }
        return std::sqrt(static_cast<long double>(bins_.rbegin()->first));
    }

    /**
     * @brief Returns the smallest Euclidean distance among equally frequent modes.
     */
    [[nodiscard]] long double mode() const noexcept {
        SquaredDistance bestDistance = 0;
        std::uint64_t bestCount = 0;
        for (const auto& [squaredDistance, binCount] : bins_) {
            if (binCount > bestCount) {
                bestDistance = squaredDistance;
                bestCount = binCount;
            }
        }
        return std::sqrt(static_cast<long double>(bestDistance));
    }

    /**
     * @brief Returns Shannon entropy of the normalized histogram in bits.
     */
    [[nodiscard]] long double entropyBits() const noexcept {
        if (count_ == 0) {
            return 0.0L;
        }
        const long double denominator = static_cast<long double>(count_);
        long double entropy = 0.0L;
        for (const auto& [squaredDistance, binCount] : bins_) {
            static_cast<void>(squaredDistance);
            const long double probability = static_cast<long double>(binCount) / denominator;
            entropy -= probability * std::log2(probability);
        }
        return entropy;
    }

  private:
    Storage bins_;
    std::uint64_t count_ = 0;
};

} // namespace mmcfilters::attributes::computers::detail::distance_transform
