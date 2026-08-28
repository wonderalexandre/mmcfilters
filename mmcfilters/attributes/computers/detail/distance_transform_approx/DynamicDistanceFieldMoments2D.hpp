#pragma once

#include "../distance_transform/DistanceFieldMoments.hpp"
#include "../../../../utils/Common.hpp"

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcfilters::attributes::computers::detail::distance_transform_approx {

/**
 * @brief Zero-cost observer used when only the paper's Bedt maximum is requested.
 */
class NoopDistanceFieldObserver2D {
  public:
    explicit NoopDistanceFieldObserver2D(int, int = 0) noexcept {}
    void insertPixel(PixelId) noexcept {}
    void mergePixels(PixelId, PixelId) noexcept {}
    void removeFiniteCost(PixelId, std::int64_t) noexcept {}
    void addFiniteCost(PixelId, std::int64_t) noexcept {}
};

/**
 * @brief Union-find summaries of finite DIFT labels in active A8 components.
 *
 * Nodes processed in one topological level are pairwise incomparable. Their
 * support pixels coexist in the global DIFT workspace, while LCA-gated edge
 * activation keeps their union-find components separate until the parent is
 * processed. Cost replacement is additive and therefore remains O(1), apart
 * from inverse-Ackermann union-find operations.
 */
template <bool ValidateInternalOperations = true> class BasicDynamicDistanceFieldMoments2D {
  public:
    explicit BasicDynamicDistanceFieldMoments2D(int numPixels, int = 0)
        : parent_(static_cast<std::size_t>(numPixels), InvalidPixel), supportSizes_(static_cast<std::size_t>(numPixels), 0),
          moments_(static_cast<std::size_t>(numPixels)) {}

    void insertPixel(PixelId pixel) {
        const std::size_t index = checkedIndex(pixel);
        if constexpr (ValidateInternalOperations) {
            if (parent_[index] != InvalidPixel) {
                throw std::logic_error("Dynamic distance-field observer received a duplicate support pixel.");
            }
        }
        parent_[index] = pixel;
        supportSizes_[index] = 1;
        moments_[index].clear();
    }

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
        moments_[firstIndex].template merge<ValidateInternalOperations>(moments_[secondIndex]);
        moments_[secondIndex].clear();
    }

    void removeFiniteCost(PixelId pixel, std::int64_t cost) {
        moments_[static_cast<std::size_t>(findMutable(pixel))].template remove<ValidateInternalOperations>(cost);
    }

    void addFiniteCost(PixelId pixel, std::int64_t cost) {
        moments_[static_cast<std::size_t>(findMutable(pixel))].template add<ValidateInternalOperations>(cost);
    }

    [[nodiscard]] const distance_transform::DistanceFieldMoments& momentsFor(PixelId pixel) const {
        return moments_[static_cast<std::size_t>(find(pixel))];
    }

    [[nodiscard]] std::uint64_t supportSizeFor(PixelId pixel) const {
        return supportSizes_[static_cast<std::size_t>(find(pixel))];
    }

  private:
    [[nodiscard]] std::size_t checkedIndex(PixelId pixel) const {
        if constexpr (ValidateInternalOperations) {
            if (pixel < 0 || static_cast<std::size_t>(pixel) >= parent_.size()) {
                throw std::out_of_range("Dynamic distance-field observer received an invalid pixel id.");
            }
        }
        return static_cast<std::size_t>(pixel);
    }

    [[nodiscard]] PixelId find(PixelId pixel) const {
        PixelId root = pixel;
        std::size_t index = checkedIndex(root);
        if constexpr (ValidateInternalOperations) {
            if (parent_[index] == InvalidPixel) {
                throw std::logic_error("Dynamic distance-field observer requires an inserted support pixel.");
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

    std::vector<PixelId> parent_;
    std::vector<std::uint64_t> supportSizes_;
    std::vector<distance_transform::DistanceFieldMoments> moments_;
};

using DynamicDistanceFieldMoments2D = BasicDynamicDistanceFieldMoments2D<true>;
using UncheckedDynamicDistanceFieldMoments2D = BasicDynamicDistanceFieldMoments2D<false>;

} // namespace mmcfilters::attributes::computers::detail::distance_transform_approx
