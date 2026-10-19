#pragma once

#include "MorphologicalTreeSemantics.hpp"
#include "NativeHierarchyOwned.hpp"
#include "../utils/Altitude.hpp"
#include "../utils/Common.hpp"
#include "../utils/Image.hpp"
#include "../utils/RegularGridAdjacency2D.hpp"
#include "../utils/CommittedGridAccess.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace mmcfilters {

/**
 * @brief Builds max-trees and min-trees with union-find.
 *
 * The builder owns its adjacency. `build()` constructs the pixel-parent
 * representation, computes the intrinsic attributes, and transfers the
 * native hierarchy buffers and adjacency to the result, consuming the builder.
 *
 * Construction and canonicalization follow [BGLWBB07], Section 3.1, Figure 5,
 * also presented in [CG14], Section III-A, Figure 3. This is the union-find
 * variant without union by rank. The min-tree uses the dual gray-level order.
 * Dense node IDs and ownership transfer are native representation choices.
 *
 * @par Scientific references
 * - [BGLWBB07] Christophe Berger, Thierry Géraud, Roland Levillain,
 *   Nicolas Widynski, Anthony Baillard, and Emmanuel Bertin,
 *   "Effective Component Tree Computation with Application to
 *   Pattern Recognition in Astronomical Imaging," 2007 IEEE International
 *   Conference on Image Processing (ICIP), vol. 4, pp. IV-41-IV-44, 2007.
 *   [DOI: 10.1109/ICIP.2007.4379949](https://doi.org/10.1109/ICIP.2007.4379949).
 * - [CG14] Edwin Carlinet and Thierry Géraud, "A Comparative Review of
 *   Component Tree Computation Algorithms," IEEE Transactions on Image
 *   Processing, vol. 23, no. 9, pp. 3885-3895, 2014.
 *   [DOI: 10.1109/TIP.2014.2336551](https://doi.org/10.1109/TIP.2014.2336551).
 */
class ComponentTreeBuilder {
    /** @brief Owned adjacency transferred to a completed native hierarchy. */
    RegularGridAdjacency2D adjacency_;
    /** @brief Indicates whether the union-find builds a max-tree. */
    bool isMaxTree_;

    /**
     * @brief Sorts image pixels in the order required by component-tree construction.
     *
     * orderedPixels is increasing for a max-tree and decreasing for a min-tree;
     * buildTree() traverses it backward during union-find. This is the S
     * orientation of [CG14], Figure 3, opposite to R in [BGLWBB07], Figure 5.
     * The local sorting policy uses counting sort for uint8 and stable comparison
     * sort for other types. Equal values retain row-major order in both cases.
     *
     * @param image Image.
     * @return Pixels ordered by gray level with stable ties.
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

    /**
     * @brief Builds and canonizes the pixel-parent representation with union-find.
     *
     * Implements COMPUTE-TREE and CANONIZE-TREE from [BGLWBB07], Section 3.1,
     * Figure 5, using the traversal orientation of [CG14], Section III-A1
     * and III-A3, Figure 3. unionParent is the auxiliary zpar forest;
     * pixelParent stores the hierarchy. Root lookup uses iterative path halving
     * on unionParent, while canonicalization updates pixelParent in forward
     * order. There is no union by rank or level-compression optimization.
     *
     * @param image Image.
     * @return Pixel parents, ordered pixels, and the number of canonical nodes.
     */
    template <typename PixelType>
    [[nodiscard]] std::tuple<std::vector<PixelId>, std::vector<PixelId>, int> buildTree(const ImagePtr<PixelType>& image) const {
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
            for (PixelId neighbor : detail::CommittedGridAccess::neighbors(adjacency_, pixel)) {
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

    /**
     * @brief Materializes the three native arrays and transfers ownership with the adjacency.
     *
     * Stable sorting visits equal-valued pixels in row-major order. The first
     * pixel in each direct proper part is canonical and supplies the node
     * altitude, preserving equal floating-point values with different representations.
     *
     * Canonical pixels satisfy the representative criterion in [CG14],
     * Section II-B: the root, or a pixel whose value differs from its parent.
     * Converting these representatives to dense node IDs, mapping each source
     * pixel to its node, and transferring adjacency ownership are local steps.
     *
     * @param image Source image.
     * @param pixelParent Canonized pixel parents.
     * @param orderedPixels Pixels ordered with each canonical parent before its descendants.
     * @param numNodes Number of canonical nodes.
     * @return Owning hierarchy independent of the builder and source image.
     */
    template <AltitudeValue T>
    [[nodiscard]] NativeHierarchyOwned<T> computeIntrinsicAttributes(const ImagePtr<T>& image, const std::vector<PixelId>& pixelParent,
                                                                   const std::vector<PixelId>& orderedPixels, int numNodes) {
        const MorphologicalTreeKind kind = isMaxTree_ ? MorphologicalTreeKind::MaxTree : MorphologicalTreeKind::MinTree;
        const GridDomain2D gridDomain{image->getNumRows(), image->getNumColumns()};
        const auto* imageValues = image->rawData();
        std::vector<NodeId> nodeParent(static_cast<std::size_t>(numNodes));
        std::vector<T> nodeValuation(static_cast<std::size_t>(numNodes));
        std::vector<NodeId> smallestNode(static_cast<std::size_t>(image->getSize()));

        NodeId idNode = 0;
        for (const PixelId pixel : orderedPixels) {
            const PixelId parent = pixelParent[pixel];
            if (pixel == parent || imageValues[pixel] != imageValues[parent]) {
                assert(pixel != parent || idNode == 0);
                smallestNode[pixel] = idNode;
                nodeParent[idNode] = pixel == parent ? idNode : smallestNode[parent];
                nodeValuation[idNode] = imageValues[pixel];
                ++idNode;
            } else {
                smallestNode[pixel] = smallestNode[parent];
            }
        }

        assert(idNode == numNodes);
        return makeNativeHierarchyOwned<T>(std::move(nodeParent), std::move(smallestNode), std::move(nodeValuation), 0, gridDomain,
                                           makeMorphologicalTreeSemantics(kind, SharedAdjacencyContext{std::move(adjacency_)}));
    }

  public:
    /**
     * @brief Configures max-tree or min-tree construction with an owned adjacency.
     * @param adjacency Adjacency relation taken by value.
     * @param isMaxTree Whether to build a max-tree rather than a min-tree.
     */
    explicit ComponentTreeBuilder(RegularGridAdjacency2D adjacency, bool isMaxTree) noexcept
        : adjacency_(std::move(adjacency)), isMaxTree_(isMaxTree) {}

    /**
     * @brief Builds and returns the native hierarchy, consuming the builder.
     *
     * The result owns the adjacency after construction; named builders must be
     * passed through `std::move` to make that transfer explicit.
     *
     * @param image Non-null, non-empty image with finite values and the same domain as the connected adjacency.
     * @return Owning hierarchy independent of the builder and source image.
     */
    template <AltitudeValue T> [[nodiscard]] NativeHierarchyOwned<T> build(const ImagePtr<T>& image) && {
        assert(adjacency_.getNumRows() == image->getNumRows() && adjacency_.getNumColumns() == image->getNumColumns());
        auto [pixelParent, orderedPixels, numNodes] = buildTree(image);
        return computeIntrinsicAttributes(image, pixelParent, orderedPixels, numNodes);
    }
};

} // namespace mmcfilters
