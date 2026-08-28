#pragma once

#include "component_tree_dift/EdtDIFT.hpp"
#include "mmcfilters/contours/ContoursComputedIncrementally.hpp"
#include "mmcfilters/trees/MorphologicalTree.hpp"
#include "mmcfilters/utils/Altitude.hpp"
#include "mmcfilters/utils/Common.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace mmcfilters::unit_tests::maxdist_oracle {

/** @brief Exact test-only EDT samples for one projected node support. */
struct ExactNodeDistanceTransform {
    std::vector<PixelId> supportPixels;
    std::vector<PixelId> contourPixels;
    std::vector<std::int64_t> squaredDistances;
};

/**
 * @brief Reconstructs one complete node EDT by independent brute force.
 *
 * Support pixels and squared distances use matching indices. This richer
 * oracle is the reference for future reducers that need more than the maximum.
 */
inline ExactNodeDistanceTransform exactNodeSquaredDistanceTransform(const MorphologicalTree& tree, NodeId node) {
    if (!tree.hasGridDomain2D()) {
        throw std::invalid_argument("Exact node EDT oracle requires a regular 2D domain.");
    }
    if (!tree.isAlive(node)) {
        throw std::out_of_range("Exact node EDT oracle requires a live node.");
    }

    const int rows = tree.numRows();
    const int columns = tree.numColumns();
    const int numPixels = rows * columns;
    std::vector<std::uint8_t> supportMask(static_cast<std::size_t>(numPixels), 0);
    ExactNodeDistanceTransform result;

    for (PixelId pixel : tree.nodeSupport(node)) {
        result.supportPixels.push_back(pixel);
        supportMask[static_cast<std::size_t>(pixel)] = 1;
    }

    for (PixelId pixel : result.supportPixels) {
        const int row = pixel / columns;
        const int column = pixel % columns;
        const bool touchesDomainBoundary = row == 0 || row == rows - 1 || column == 0 || column == columns - 1;
        const bool touchesComplement =
            (!touchesDomainBoundary &&
             (supportMask[static_cast<std::size_t>(pixel - columns)] == 0 || supportMask[static_cast<std::size_t>(pixel + columns)] == 0 ||
              supportMask[static_cast<std::size_t>(pixel - 1)] == 0 || supportMask[static_cast<std::size_t>(pixel + 1)] == 0));
        if (touchesDomainBoundary || touchesComplement) {
            result.contourPixels.push_back(pixel);
        }
    }

    if (result.supportPixels.empty() || result.contourPixels.empty()) {
        throw std::logic_error("Exact node EDT oracle requires non-empty support and foreground contour.");
    }

    result.squaredDistances.reserve(result.supportPixels.size());
    for (PixelId pixel : result.supportPixels) {
        const std::int64_t row = pixel / columns;
        const std::int64_t column = pixel % columns;
        std::int64_t minSquaredDistance = std::numeric_limits<std::int64_t>::max();
        for (PixelId boundaryPixel : result.contourPixels) {
            const std::int64_t boundaryRow = boundaryPixel / columns;
            const std::int64_t boundaryColumn = boundaryPixel % columns;
            const std::int64_t deltaRow = row - boundaryRow;
            const std::int64_t deltaColumn = column - boundaryColumn;
            minSquaredDistance = std::min(minSquaredDistance, deltaRow * deltaRow + deltaColumn * deltaColumn);
        }
        result.squaredDistances.push_back(minSquaredDistance);
    }

    return result;
}

/**
 * @brief Independent brute-force squared Euclidean MAX_DIST_EXACT oracle.
 *
 * This implementation intentionally reconstructs support masks and foreground
 * A4 contours without using the production distance-transform backend or the
 * incremental contour implementation.
 */
inline std::vector<std::int64_t> exactMaxSquaredDistance(const MorphologicalTree& tree) {
    if (!tree.hasGridDomain2D()) {
        throw std::invalid_argument("MAX_DIST_EXACT oracle requires a regular 2D domain.");
    }

    std::vector<std::int64_t> result(static_cast<std::size_t>(tree.numInternalNodeSlots()), 0);

    for (NodeId node : tree.aliveNodeIds()) {
        const ExactNodeDistanceTransform nodeTransform = exactNodeSquaredDistanceTransform(tree, node);
        result[static_cast<std::size_t>(node)] = *std::max_element(nodeTransform.squaredDistances.begin(), nodeTransform.squaredDistances.end());
    }

    return result;
}

/**
 * @brief Test-only optimized component-tree DIFT oracle.
 *
 * This preserves the former altitude-level sweep solely for regression and
 * cross-validation of max-trees/min-trees. Production code never dispatches to
 * this implementation.
 */
template <AltitudeValue T> inline std::vector<int> componentTreeDiftMaxSquaredDistance(const MorphologicalTree& tree, std::span<const T> altitude) {
    if (tree.nodeAltitudeOrder() == NodeAltitudeOrder::Unconstrained) {
        throw std::invalid_argument("Component-tree DIFT oracle requires monotone node altitudes.");
    }
    if (altitude.size() != static_cast<std::size_t>(tree.numInternalNodeSlots())) {
        throw std::invalid_argument("Component-tree DIFT oracle received an altitude buffer with invalid size.");
    }

    using EdtDIFT = component_tree_dift::EdtDIFT;
    const int numNodes = tree.numInternalNodeSlots();
    const int numPixels = tree.numRows() * tree.numColumns();
    EdtDIFT transform(tree.numRows(), tree.numColumns());
    const auto contourDeltas = ContoursComputedIncrementally::extractContourDeltas(tree);
    std::vector<std::vector<PixelId>> contours(static_cast<std::size_t>(numNodes));
    std::vector<std::uint8_t> removalMark(static_cast<std::size_t>(numPixels), 0);
    std::vector<std::uint8_t> additionMark(static_cast<std::size_t>(numPixels), 0);
    std::vector<int> result(static_cast<std::size_t>(numNodes), 0);

    std::vector<NodeId> nodes;
    nodes.reserve(static_cast<std::size_t>(tree.numNodes()));
    for (NodeId node : tree.postOrder()) {
        nodes.push_back(node);
    }
    std::stable_sort(nodes.begin(), nodes.end(), [&](NodeId lhs, NodeId rhs) {
        const T lhsAltitude = altitude[static_cast<std::size_t>(lhs)];
        const T rhsAltitude = altitude[static_cast<std::size_t>(rhs)];
        return tree.nodeAltitudeOrder() == NodeAltitudeOrder::Increasing ? lhsAltitude > rhsAltitude : lhsAltitude < rhsAltitude;
    });

    std::vector<PixelId> removedPixels;
    std::size_t groupBegin = 0;
    while (groupBegin < nodes.size()) {
        std::size_t groupEnd = groupBegin + 1;
        while (groupEnd < nodes.size() && altitude[static_cast<std::size_t>(nodes[groupBegin])] == altitude[static_cast<std::size_t>(nodes[groupEnd])]) {
            ++groupEnd;
        }

        for (std::size_t index = groupBegin; index < groupEnd; ++index) {
            const NodeId node = nodes[index];
            std::vector<PixelId>& nodeContour = contours[static_cast<std::size_t>(node)];
            nodeContour.clear();

            const auto removals = contourDeltas.removals(node);
            removedPixels.assign(removals.begin(), removals.end());
            for (PixelId pixel : removals) {
                removalMark[static_cast<std::size_t>(pixel)] = 1;
            }

            for (NodeId child : tree.children(node)) {
                std::vector<PixelId>& childContour = contours[static_cast<std::size_t>(child)];
                for (PixelId pixel : childContour) {
                    if (removalMark[static_cast<std::size_t>(pixel)] == 0) {
                        nodeContour.push_back(pixel);
                    }
                }
                std::vector<PixelId>().swap(childContour);
            }
            for (PixelId pixel : removals) {
                removalMark[static_cast<std::size_t>(pixel)] = 0;
            }
            if (!removedPixels.empty()) {
                transform.treeRemoval(removedPixels);
            }

            const auto additions = contourDeltas.additions(node);
            for (PixelId pixel : additions) {
                additionMark[static_cast<std::size_t>(pixel)] = 1;
            }
            for (PixelId pixel : tree.properPart(node)) {
                transform.addPixelToBinaryImage(pixel);
                if (additionMark[static_cast<std::size_t>(pixel)] != 0) {
                    nodeContour.push_back(pixel);
                    transform.seed(pixel);
                } else {
                    transform.open(pixel);
                    transform.insertNeighborsPQueue(pixel);
                }
            }
            for (PixelId pixel : additions) {
                additionMark[static_cast<std::size_t>(pixel)] = 0;
            }
        }

        transform.run();
        for (std::size_t index = groupBegin; index < groupEnd; ++index) {
            const NodeId node = nodes[index];
            result[static_cast<std::size_t>(node)] = transform.maxBedt(contours[static_cast<std::size_t>(node)]);
        }
        groupBegin = groupEnd;
    }

    return result;
}

} // namespace mmcfilters::unit_tests::maxdist_oracle
