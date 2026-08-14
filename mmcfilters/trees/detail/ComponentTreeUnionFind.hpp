#pragma once

#include "../../utils/Common.hpp"
#include "../../utils/Image.hpp"
#include "../../utils/RegularGridAdjacency2D.hpp"
#include "../../utils/CommittedGridAccess.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <vector>

namespace mmcfilters::detail {

/**
 * @brief Union-find construction used by the concrete max-tree/min-tree producer.
 *
 * This is an implementation detail rather than a polymorphic morphological-tree
 * builder. The construction owns no hierarchy state and borrows its adjacency
 * for the duration of one build.
 */
class ComponentTreeUnionFind {
    /** @brief Adjacency. */
    const RegularGridAdjacency2D* adjacency_;
    /** @brief Indicates whether the union-find builds a max-tree. */
    bool isMaxTree_;

  public:
    /**
     * @brief Constructs `ComponentTreeUnionFind` from the supplied inputs.
     *
     * @param adjacency Adjacency relation.
     * @param isMaxTree Flag controlling is max tree.
     */
    explicit ComponentTreeUnionFind(const RegularGridAdjacency2D* adjacency, bool isMaxTree) noexcept : adjacency_(adjacency), isMaxTree_(isMaxTree) {}

    /**
     * @brief Sorts image pixels in the order required by component-tree construction.
     *
     * @param image Image.
     * @return Values produced by the operation.
     */
    template <typename PixelType> [[nodiscard]] std::vector<PixelId> sort(const ImagePtr<PixelType>& image) const {
        const int numPixels = image->getSize();
        std::vector<PixelId> orderedPixels(static_cast<std::size_t>(numPixels));
        const PixelType* values = image->rawData();

        if constexpr (!std::is_same_v<PixelType, std::uint8_t>) {
            if (PRINT_LOG) {
                std::cout << "Sorting image with comparison sort, size: " << numPixels << std::endl;
            }
            std::iota(orderedPixels.begin(), orderedPixels.end(), 0);
            if (isMaxTree_) {
                std::stable_sort(orderedPixels.begin(), orderedPixels.end(), [&](PixelId lhs, PixelId rhs) { return values[lhs] < values[rhs]; });
            } else {
                std::stable_sort(orderedPixels.begin(), orderedPixels.end(), [&](PixelId lhs, PixelId rhs) { return values[lhs] > values[rhs]; });
            }
        } else {
            if (PRINT_LOG) {
                std::cout << "Sorting uint8 image with counting sort, size: " << numPixels << std::endl;
            }

            int maxValue = static_cast<int>(values[0]);
            for (int index = 1; index < numPixels; ++index) {
                if (maxValue < values[index]) {
                    maxValue = values[index];
                }
            }

            std::vector<std::uint32_t> counter(static_cast<std::size_t>(maxValue) + 1, 0);
            if (isMaxTree_) {
                for (int index = 0; index < numPixels; ++index) {
                    ++counter[values[index]];
                }
                for (int level = 1; level <= maxValue; ++level) {
                    counter[level] += counter[level - 1];
                }
                for (int index = numPixels - 1; index >= 0; --index) {
                    orderedPixels[--counter[values[index]]] = index;
                }
            } else {
                for (int index = 0; index < numPixels; ++index) {
                    ++counter[maxValue - values[index]];
                }
                for (int level = 1; level <= maxValue; ++level) {
                    counter[level] += counter[level - 1];
                }
                for (int index = numPixels - 1; index >= 0; --index) {
                    orderedPixels[--counter[maxValue - values[index]]] = index;
                }
            }
        }
        return orderedPixels;
    }

    template <typename PixelType>
    [[nodiscard]] std::tuple<std::vector<PixelId>, std::vector<PixelId>, int>
    /**
     * @brief Builds a component tree from the input image.
     *
     * @param image Image.
     * @return Component tree built from the input image.
     */
    build(const ImagePtr<PixelType>& image) const {
        std::vector<PixelId> orderedPixels = sort(image);
        const PixelType* values = image->rawData();

        const int numPixels = image->getSize();
        std::vector<PixelId> unionParent(static_cast<std::size_t>(numPixels), InvalidPixel);
        std::vector<PixelId> pixelParent(static_cast<std::size_t>(numPixels), InvalidPixel);
        auto findRoot = [&](PixelId pixel) {
            while (unionParent[pixel] != pixel) {
                unionParent[pixel] = unionParent[unionParent[pixel]];
                pixel = unionParent[pixel];
            }
            return pixel;
        };

        for (int index = numPixels - 1; index >= 0; --index) {
            const PixelId pixel = orderedPixels[static_cast<std::size_t>(index)];
            pixelParent[pixel] = pixel;
            unionParent[pixel] = pixel;
            for (PixelId neighbor : CommittedGridAccess::neighbors(*adjacency_, pixel)) {
                if (unionParent[neighbor] == InvalidPixel) {
                    continue;
                }
                const PixelId neighborRoot = findRoot(neighbor);
                if (pixel != neighborRoot) {
                    pixelParent[neighborRoot] = pixel;
                    unionParent[neighborRoot] = pixel;
                }
            }
        }

        int numNodes = 0;
        for (int index = 0; index < numPixels; ++index) {
            const PixelId pixel = orderedPixels[static_cast<std::size_t>(index)];
            const PixelId parent = pixelParent[pixel];
            if (values[pixelParent[parent]] == values[parent]) {
                pixelParent[pixel] = pixelParent[parent];
            }
            if (pixelParent[pixel] == pixel || values[pixelParent[pixel]] != values[pixel]) {
                ++numNodes;
            }
        }

        return {std::move(pixelParent), std::move(orderedPixels), numNodes};
    }
};

} // namespace mmcfilters::detail
