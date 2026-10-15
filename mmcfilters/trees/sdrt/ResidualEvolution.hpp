#pragma once

/**
 * @file ResidualEvolution.hpp
 * @brief Normative values used by self-dual residual-tree evolution.
 */

#include "../../utils/Altitude.hpp"
#include "../../utils/Common.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <set>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters::sdrt {

/** @brief Polarity of a regional extremum in synchronized residual evolution. */
enum class Polarity { Maximum, Minimum };

/**
 * @brief Total order over the pixels of one construction domain.
 *
 * An explicit order is represented by a permutation of all pixel identifiers.
 * The default row-major specialization avoids storing that permutation and uses
 * the dense pixel identifier itself as its rank.
 */
class SpatialOrder {
  public:
    /**
     * @brief Creates an explicit total spatial order.
     * @param pixelsInOrder Pixel identifiers from first to last.
     * @throws std::invalid_argument if the sequence is not a permutation of `[0,n)`.
     */
    explicit SpatialOrder(std::vector<PixelId> pixelsInOrder) : rankByPixel_(pixelsInOrder.size(), invalidRank()) {
        for (std::size_t rank = 0; rank < pixelsInOrder.size(); ++rank) {
            const PixelId pixel = pixelsInOrder[rank];
            if (pixel < 0 || static_cast<std::size_t>(pixel) >= pixelsInOrder.size()) {
                throw std::invalid_argument("A spatial order must be a permutation of the pixel domain.");
            }
            auto& storedRank = rankByPixel_[static_cast<std::size_t>(pixel)];
            if (storedRank != invalidRank()) {
                throw std::invalid_argument("A spatial order cannot contain a pixel more than once.");
            }
            storedRank = rank;
        }
    }

    /**
     * @brief Returns whether the order uses implicit row-major ranks.
     * @return `true` for the implicit row-major order.
     */
    [[nodiscard]] bool isRowMajor() const noexcept { return rankByPixel_.empty(); }

    /**
     * @brief Validates that this order covers a construction domain.
     * @param numPixels Number of pixels in the domain.
     */
    void validateDomainSize(std::size_t numPixels) const {
        if (!isRowMajor() && rankByPixel_.size() != numPixels) {
            throw std::invalid_argument("The explicit spatial order size differs from the image domain.");
        }
    }

    /**
     * @brief Returns whether `lhs` precedes `rhs`.
     * @param lhs First pixel.
     * @param rhs Second pixel.
     * @return `true` when `lhs` has a lower spatial rank than `rhs`.
     */
    [[nodiscard]] bool precedes(PixelId lhs, PixelId rhs) const { return rank(lhs) < rank(rhs); }

    /**
     * @brief Returns the least pixel of a non-empty support.
     * @param support Pixel support whose spatial minimum is requested.
     * @return Pixel of `support` with the lowest configured spatial rank.
     */
    [[nodiscard]] PixelId spatialMinimum(std::span<const PixelId> support) const {
        if (support.empty()) {
            throw std::invalid_argument("A spatial minimum requires a non-empty support.");
        }
        PixelId minimum = support.front();
        for (PixelId pixel : support.subspan(1)) {
            if (precedes(pixel, minimum)) {
                minimum = pixel;
            }
        }
        return minimum;
    }

  protected:
    /** @brief Creates the implicit row-major order. */
    SpatialOrder() = default;

  private:
    /** @return Sentinel used for a rank that has not been assigned. */
    [[nodiscard]] static constexpr std::size_t invalidRank() noexcept { return std::numeric_limits<std::size_t>::max(); }

    /**
     * @brief Returns the configured rank of a pixel.
     * @param pixel Pixel whose rank is requested.
     * @return Zero-based spatial rank.
     */
    [[nodiscard]] std::size_t rank(PixelId pixel) const {
        if (pixel < 0) {
            throw std::out_of_range("A spatial-order query received a negative pixel identifier.");
        }
        if (isRowMajor()) {
            return static_cast<std::size_t>(pixel);
        }
        if (static_cast<std::size_t>(pixel) >= rankByPixel_.size()) {
            throw std::out_of_range("A spatial-order query lies outside the configured domain.");
        }
        return rankByPixel_[static_cast<std::size_t>(pixel)];
    }

    std::vector<std::size_t> rankByPixel_; ///< Explicit rank indexed by pixel; empty means row-major.
};

/** @brief Dense row-major spatial order used by default. */
class RowMajorSpatialOrder final : public SpatialOrder {
  public:
    RowMajorSpatialOrder() = default;
};

/** @brief Contrast-invariant key of one current residual extremum. */
struct SelfDualResidualKey {
    std::size_t supportCardinality = 0;    ///< Cardinality of the current extremum support.
    PixelId spatialMinimum = InvalidPixel; ///< Least support pixel in the configured spatial order.

    /** @brief Compares both coordinates of two residual keys. */
    friend bool operator==(const SelfDualResidualKey&, const SelfDualResidualKey&) = default;
};

/** @brief Strict order induced by `(supportCardinality, spatialMinimum)`. */
class SelfDualResidualOrder {
  public:
    /**
     * @brief Creates the canonical residual-key order.
     * @param spatialOrder Total order used by the second key coordinate.
     */
    explicit SelfDualResidualOrder(SpatialOrder spatialOrder = RowMajorSpatialOrder{}) : spatialOrder_(std::move(spatialOrder)) {}

    /**
     * @brief Compares two extremum keys without consulting polarity or altitude.
     * @param lhs Left extremum key.
     * @param rhs Right extremum key.
     * @return `true` when `lhs` precedes `rhs`.
     */
    [[nodiscard]] bool compareExtremumKeys(const SelfDualResidualKey& lhs, const SelfDualResidualKey& rhs) const {
        if (lhs.supportCardinality != rhs.supportCardinality) {
            return lhs.supportCardinality < rhs.supportCardinality;
        }
        if (lhs.spatialMinimum == rhs.spatialMinimum) {
            return false;
        }
        return spatialOrder_.precedes(lhs.spatialMinimum, rhs.spatialMinimum);
    }

    /**
     * @brief Applies the canonical strict comparison.
     * @param lhs Left extremum key.
     * @param rhs Right extremum key.
     * @return `true` when `lhs` precedes `rhs`.
     */
    [[nodiscard]] bool operator()(const SelfDualResidualKey& lhs, const SelfDualResidualKey& rhs) const { return compareExtremumKeys(lhs, rhs); }

    /**
     * @brief Returns the total spatial order used by this comparison.
     * @return Configured total pixel order.
     */
    [[nodiscard]] const SpatialOrder& spatialOrder() const noexcept { return spatialOrder_; }

  private:
    SpatialOrder spatialOrder_; ///< Total pixel order defining the second key coordinate.
};

/** @brief Selects the next residual extremum from the canonical self-dual order. */
class SelfDualResidualSchedule {
  public:
    /**
     * @brief Creates the unique canonical residual schedule.
     * @param spatialOrder Total order used by equal-cardinality extrema.
     */
    explicit SelfDualResidualSchedule(SpatialOrder spatialOrder = RowMajorSpatialOrder{}) : residualOrder_(std::move(spatialOrder)) {}

    /**
     * @brief Returns the index of the first key in canonical residual order.
     * @param residualKeys Non-empty extremum-key sequence.
     * @return Index of the least extremum key.
     */
    [[nodiscard]] std::size_t selectSmallestExtremumIndex(std::span<const SelfDualResidualKey> residualKeys) const {
        if (residualKeys.empty()) {
            throw std::invalid_argument("Residual-extremum selection requires a non-empty schedule.");
        }
        std::set<std::pair<std::size_t, PixelId>> seenKeys;
        std::size_t selected = 0;
        for (std::size_t index = 0; index < residualKeys.size(); ++index) {
            if (!seenKeys.emplace(residualKeys[index].supportCardinality, residualKeys[index].spatialMinimum).second) {
                throw std::invalid_argument("A self-dual residual schedule cannot contain duplicate extremum keys.");
            }
            if (index != 0 && residualOrder_.compareExtremumKeys(residualKeys[index], residualKeys[selected])) {
                selected = index;
            }
        }
        return selected;
    }

    /**
     * @brief Returns the canonical key order used by this schedule.
     * @return Residual-key order owned by the schedule.
     */
    [[nodiscard]] const SelfDualResidualOrder& residualOrder() const noexcept { return residualOrder_; }

  private:
    SelfDualResidualOrder residualOrder_; ///< Canonical contrast-invariant extremum order.
};

/** @brief Options that affect extremum residual-tree construction. */
struct ExtremumResidualTreeOptions {
    SpatialOrder spatialOrder = RowMajorSpatialOrder{}; ///< Total order defining each support's spatial minimum.
};

} // namespace mmcfilters::sdrt
