#pragma once

#include "../distance_transform/DistanceFieldHistogram.hpp"
#include "../distance_transform/DistanceFieldMoments.hpp"
#include "../distance_transform/DistanceWeightedSpatialMoments.hpp"
#include "../../../../utils/Common.hpp"

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace mmcfilters::attributes::computers::detail::distance_transform_approx {

/**
 * @brief Selective union-find summaries for approximate distance fields.
 *
 * The three boolean policies remove unused state and update work at compile
 * time. At least one statistic must be enabled. This lets profile-only,
 * spatial-only, and mixed attribute requests share the same dynamic DIFT
 * traversal without paying for unrelated reducers.
 *
 * @tparam TrackMoments Retain scalar distance moments.
 * @tparam TrackHistogram Retain the sparse squared-distance histogram.
 * @tparam TrackSpatialMoments Retain distance-weighted spatial moments.
 * @tparam ValidateInternalOperations Check trusted union-find operations.
 */
template <bool TrackMoments, bool TrackHistogram, bool TrackSpatialMoments, bool ValidateInternalOperations = true>
class BasicDynamicDistanceFieldStatistics2D {
    static_assert(TrackMoments || TrackHistogram || TrackSpatialMoments,
                  "A dynamic distance-field observer must retain at least one statistic.");

    using MomentsStorage =
        std::conditional_t<TrackMoments, std::vector<distance_transform::DistanceFieldMoments>, std::monostate>;
    using HistogramStorage =
        std::conditional_t<TrackHistogram, std::vector<distance_transform::DistanceFieldHistogram>, std::monostate>;
    using SpatialStorage =
        std::conditional_t<TrackSpatialMoments, std::vector<distance_transform::DistanceWeightedSpatialMoments>, std::monostate>;

    [[nodiscard]] static std::size_t allocationSize(int numPixels) {
        if constexpr (ValidateInternalOperations) {
            if (numPixels <= 0) {
                throw std::invalid_argument("Dynamic distance-field statistics require a non-empty pixel domain.");
            }
        }
        return static_cast<std::size_t>(numPixels);
    }

  public:
    /// Compile-time statement of the retained scalar-moment state.
    inline static constexpr bool tracksMoments = TrackMoments;
    /// Compile-time statement of the retained histogram state.
    inline static constexpr bool tracksHistogram = TrackHistogram;
    /// Compile-time statement of the retained spatial-moment state.
    inline static constexpr bool tracksSpatialMoments = TrackSpatialMoments;

    /**
     * @brief Allocates the selected summaries over a dense pixel domain.
     * @param numPixels Number of addressable pixels.
     * @param numColumns Image width, required only for spatial moments.
     */
    BasicDynamicDistanceFieldStatistics2D(int numPixels, int numColumns)
        : numColumns_(numColumns), parent_(allocationSize(numPixels), InvalidPixel),
          supportSizes_(allocationSize(numPixels), 0) {
        if constexpr (TrackMoments) {
            moments_.resize(static_cast<std::size_t>(numPixels));
        }
        if constexpr (TrackHistogram) {
            histograms_.resize(static_cast<std::size_t>(numPixels));
        }
        if constexpr (TrackSpatialMoments) {
            spatialMoments_.resize(static_cast<std::size_t>(numPixels));
        }
        if constexpr (ValidateInternalOperations) {
            if constexpr (TrackSpatialMoments) {
                if (numColumns <= 0) {
                    throw std::invalid_argument("Spatial distance-field statistics require a positive column count.");
                }
            }
        }
    }

    /**
     * @brief Inserts one pixel as a singleton dynamic component.
     */
    void insertPixel(PixelId pixel) {
        const std::size_t index = checkedIndex(pixel);
        if constexpr (ValidateInternalOperations) {
            if (parent_[index] != InvalidPixel) {
                throw std::logic_error("Dynamic distance-field statistics received a duplicate support pixel.");
            }
        }
        parent_[index] = pixel;
        supportSizes_[index] = 1;
        if constexpr (TrackMoments) {
            moments_[index].clear();
        }
        if constexpr (TrackHistogram) {
            histograms_[index].clear();
        }
        if constexpr (TrackSpatialMoments) {
            spatialMoments_[index].clear();
        }
    }

    /**
     * @brief Merges two active pixel components by support size.
     */
    void mergePixels(PixelId first, PixelId second) {
        PixelId firstRoot = findMutable(first);
        PixelId secondRoot = findMutable(second);
        if (firstRoot == secondRoot) {
            return;
        }
        std::size_t firstIndex = static_cast<std::size_t>(firstRoot);
        std::size_t secondIndex = static_cast<std::size_t>(secondRoot);
        if (supportSizes_[firstIndex] < supportSizes_[secondIndex]) {
            std::swap(firstRoot, secondRoot);
            std::swap(firstIndex, secondIndex);
        }
        parent_[secondIndex] = firstRoot;
        supportSizes_[firstIndex] += supportSizes_[secondIndex];
        supportSizes_[secondIndex] = 0;
        if constexpr (TrackMoments) {
            moments_[firstIndex].template merge<ValidateInternalOperations>(moments_[secondIndex]);
            moments_[secondIndex].clear();
        }
        if constexpr (TrackHistogram) {
            histograms_[firstIndex].template merge<ValidateInternalOperations>(histograms_[secondIndex]);
        }
        if constexpr (TrackSpatialMoments) {
            spatialMoments_[firstIndex].template merge<ValidateInternalOperations>(spatialMoments_[secondIndex]);
            spatialMoments_[secondIndex].clear();
        }
    }

    /**
     * @brief Removes one finite DIFT label from its active component.
     */
    void removeFiniteCost(PixelId pixel, std::int64_t cost) {
        const std::size_t root = static_cast<std::size_t>(findMutable(pixel));
        if constexpr (TrackMoments) {
            moments_[root].template remove<ValidateInternalOperations>(cost);
        }
        if constexpr (TrackHistogram) {
            histograms_[root].template remove<ValidateInternalOperations>(cost);
        }
        if constexpr (TrackSpatialMoments) {
            spatialMoments_[root].template remove<ValidateInternalOperations>(pixel, cost, numColumns_);
        }
    }

    /**
     * @brief Adds one finite DIFT label to its active component.
     */
    void addFiniteCost(PixelId pixel, std::int64_t cost) {
        const std::size_t root = static_cast<std::size_t>(findMutable(pixel));
        if constexpr (TrackMoments) {
            moments_[root].template add<ValidateInternalOperations>(cost);
        }
        if constexpr (TrackHistogram) {
            histograms_[root].template add<ValidateInternalOperations>(cost);
        }
        if constexpr (TrackSpatialMoments) {
            spatialMoments_[root].template add<ValidateInternalOperations>(pixel, cost, numColumns_);
        }
    }

    /**
     * @brief Returns scalar moments for the component containing a pixel.
     */
    [[nodiscard]] const distance_transform::DistanceFieldMoments& momentsFor(PixelId pixel) const
        requires TrackMoments
    {
        return moments_[static_cast<std::size_t>(find(pixel))];
    }

    /**
     * @brief Returns the histogram for the component containing a pixel.
     */
    [[nodiscard]] const distance_transform::DistanceFieldHistogram& histogramFor(PixelId pixel) const
        requires TrackHistogram
    {
        return histograms_[static_cast<std::size_t>(find(pixel))];
    }

    /**
     * @brief Returns spatial moments for the component containing a pixel.
     */
    [[nodiscard]] const distance_transform::DistanceWeightedSpatialMoments& spatialMomentsFor(PixelId pixel) const
        requires TrackSpatialMoments
    {
        return spatialMoments_[static_cast<std::size_t>(find(pixel))];
    }

    /**
     * @brief Returns the support cardinality of a pixel component.
     */
    [[nodiscard]] std::uint64_t supportSizeFor(PixelId pixel) const {
        return supportSizes_[static_cast<std::size_t>(find(pixel))];
    }

  private:
    [[nodiscard]] std::size_t checkedIndex(PixelId pixel) const {
        if constexpr (ValidateInternalOperations) {
            if (pixel < 0 || static_cast<std::size_t>(pixel) >= parent_.size()) {
                throw std::out_of_range("Dynamic distance-field statistics received an invalid pixel id.");
            }
        }
        return static_cast<std::size_t>(pixel);
    }

    [[nodiscard]] PixelId find(PixelId pixel) const {
        PixelId root = pixel;
        std::size_t index = checkedIndex(root);
        if constexpr (ValidateInternalOperations) {
            if (parent_[index] == InvalidPixel) {
                throw std::logic_error("Dynamic distance-field statistics require an inserted support pixel.");
            }
        }
        while (parent_[index] != root) {
            root = parent_[index];
            index = static_cast<std::size_t>(root);
        }
        return root;
    }

    [[nodiscard]] PixelId findMutable(PixelId pixel) {
        const PixelId root = find(pixel);
        PixelId current = pixel;
        while (current != root) {
            const std::size_t index = static_cast<std::size_t>(current);
            const PixelId next = parent_[index];
            parent_[index] = root;
            current = next;
        }
        return root;
    }

    int numColumns_ = 0;
    std::vector<PixelId> parent_;
    std::vector<std::uint64_t> supportSizes_;
    MomentsStorage moments_;
    HistogramStorage histograms_;
    SpatialStorage spatialMoments_;
};

/**
 * @brief Checked complete observer retained for diagnostics and compatibility.
 */
using DynamicDistanceFieldStatistics2D = BasicDynamicDistanceFieldStatistics2D<true, true, true, true>;
/**
 * @brief Production complete observer retained for compatibility.
 */
using UncheckedDynamicDistanceFieldStatistics2D = BasicDynamicDistanceFieldStatistics2D<true, true, true, false>;

} // namespace mmcfilters::attributes::computers::detail::distance_transform_approx
