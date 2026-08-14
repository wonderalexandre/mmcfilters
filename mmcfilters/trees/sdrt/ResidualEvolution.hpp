#pragma once

/**
 * @file ResidualEvolution.hpp
 * @brief Normative values used by self-dual residual-tree evolution.
 */

#include "../../utils/Altitude.hpp"
#include "../../utils/Common.hpp"

#include <cstddef>
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
    explicit SpatialOrder(std::vector<PixelId> pixelsInOrder)
        : rankByPixel_(pixelsInOrder.size(), invalidRank()) {
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
    [[nodiscard]] bool precedes(PixelId lhs, PixelId rhs) const {
        return rank(lhs) < rank(rhs);
    }

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
    [[nodiscard]] static constexpr std::size_t invalidRank() noexcept {
        return std::numeric_limits<std::size_t>::max();
    }

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

/** @brief Contrast-invariant key of one current residual candidate. */
struct SelfDualResidualKey {
    std::size_t supportCardinality = 0; ///< Cardinality of the current candidate support.
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
    explicit SelfDualResidualOrder(SpatialOrder spatialOrder = RowMajorSpatialOrder{})
        : spatialOrder_(std::move(spatialOrder)) {}

    /**
     * @brief Compares two candidate keys without consulting polarity or altitude.
     * @param lhs Left candidate key.
     * @param rhs Right candidate key.
     * @return `true` when `lhs` precedes `rhs`.
     */
    [[nodiscard]] bool compareResidualCandidates(const SelfDualResidualKey& lhs, const SelfDualResidualKey& rhs) const {
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
     * @param lhs Left candidate key.
     * @param rhs Right candidate key.
     * @return `true` when `lhs` precedes `rhs`.
     */
    [[nodiscard]] bool operator()(const SelfDualResidualKey& lhs, const SelfDualResidualKey& rhs) const {
        return compareResidualCandidates(lhs, rhs);
    }

    /**
     * @brief Returns the total spatial order used by this comparison.
     * @return Configured total pixel order.
     */
    [[nodiscard]] const SpatialOrder& spatialOrder() const noexcept { return spatialOrder_; }

  private:
    SpatialOrder spatialOrder_; ///< Total pixel order defining the second key coordinate.
};

/** @brief Selects the next residual candidate from the canonical self-dual order. */
class SelfDualResidualSchedule {
  public:
    /**
     * @brief Creates the unique canonical residual schedule.
     * @param spatialOrder Total order used by equal-cardinality candidates.
     */
    explicit SelfDualResidualSchedule(SpatialOrder spatialOrder = RowMajorSpatialOrder{})
        : residualOrder_(std::move(spatialOrder)) {}

    /**
     * @brief Returns the index of the first key in canonical residual order.
     * @param residualKeys Non-empty candidate-key sequence.
     * @return Index of the least candidate key.
     */
    [[nodiscard]] std::size_t selectResidualCandidate(std::span<const SelfDualResidualKey> residualKeys) const {
        if (residualKeys.empty()) {
            throw std::invalid_argument("Residual-candidate selection requires a non-empty schedule.");
        }
        std::set<std::pair<std::size_t, PixelId>> seenKeys;
        std::size_t selected = 0;
        for (std::size_t index = 0; index < residualKeys.size(); ++index) {
            if (!seenKeys.emplace(residualKeys[index].supportCardinality, residualKeys[index].spatialMinimum).second) {
                throw std::invalid_argument("A self-dual residual schedule cannot contain duplicate candidate keys.");
            }
            if (index != 0 && residualOrder_.compareResidualCandidates(residualKeys[index], residualKeys[selected])) {
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
    SelfDualResidualOrder residualOrder_; ///< Canonical contrast-invariant candidate order.
};

/** @brief Immutable view of an eligible residual candidate before leveling. */
template <AltitudeValue T> struct ResidualCandidate {
    std::span<const PixelId> support; ///< Current candidate support.
    Polarity polarity = Polarity::Maximum; ///< Regional-extremum polarity.
    T candidateAltitude{}; ///< Constant altitude on the candidate before leveling.
    T firstMergingLevel{}; ///< First adjacent level reached by elementary leveling.
    SelfDualResidualKey residualKey; ///< Canonical support-only scheduling key.
};

/** @brief Immutable record captured immediately before one elementary leveling. */
template <AltitudeValue T> struct ResidualEvent {
    std::size_t eventIndex = 0; ///< Zero-based chronological event index.
    std::span<const PixelId> support; ///< Candidate support in the pre-leveling state.
    Polarity polarity = Polarity::Maximum; ///< Polarity in the pre-leveling state.
    T nodeAltitude{}; ///< Candidate altitude before leveling.
    T firstMergingLevel{}; ///< Level to which the support is moved.
    AltitudeDifference<T> signedResidualValue{}; ///< `nodeAltitude - firstMergingLevel`.
};

/**
 * @brief Records the immutable event associated with one prepared candidate.
 * @param eventIndex Zero-based chronological index.
 * @param candidate Eligible candidate observed before mutation.
 * @return Complete residual event captured from the pre-leveling state.
 */
template <AltitudeValue T>
[[nodiscard]] ResidualEvent<T> recordResidualEvent(std::size_t eventIndex, const ResidualCandidate<T>& candidate) {
    if (candidate.support.empty()) {
        throw std::invalid_argument("A residual event requires a non-empty candidate support.");
    }
    if (candidate.residualKey.supportCardinality != candidate.support.size()) {
        throw std::invalid_argument("A residual candidate key has an inconsistent support cardinality.");
    }
    const auto signedResidualValue = static_cast<AltitudeDifference<T>>(candidate.candidateAltitude) -
                                     static_cast<AltitudeDifference<T>>(candidate.firstMergingLevel);
    if ((candidate.polarity == Polarity::Maximum && signedResidualValue <= AltitudeDifference<T>{}) ||
        (candidate.polarity == Polarity::Minimum && signedResidualValue >= AltitudeDifference<T>{})) {
        throw std::invalid_argument("A residual candidate polarity is inconsistent with its signed residual value.");
    }
    return ResidualEvent<T>{eventIndex, candidate.support, candidate.polarity, candidate.candidateAltitude,
                            candidate.firstMergingLevel, signedResidualValue};
}

} // namespace mmcfilters::sdrt
