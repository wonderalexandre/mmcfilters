#pragma once

#include "AttributeComputer.hpp"
#include "AttributeComputedIncrementally.hpp"
#include "../trees/MorphologicalTree.hpp"

namespace mmcfilters {

namespace detail {
inline NodeId boundingBoxSlotOf(const MorphologicalTree&, NodeId nodeId) noexcept {
    return nodeId;
}
}

/**
 * @brief Computes descriptors derived from the axis-aligned bounding box of
 * the node support.
 *
 * @details
 * For each live node, the computer determines the smallest axis-aligned box
 * that contains the complete subtree support. The support is defined as the
 * union of the node's proper parts and the supports of all descendants. From
 * the box extrema, the computer derives:
 * - width and height;
 * - minimum and maximum row/column coordinates;
 * - diagonal length;
 * - aspect ratio (`RATIO_WH`);
 * - rectangularity, defined as `AREA / (width * height)`.
 *
 * The extrema are accumulated in a post-order traversal by initialising each
 * node from its own proper parts and then merging child boxes into the parent.
 * Only `RECTANGULARITY` depends on the external `AREA` attribute; all other
 * outputs come directly from the tracked extrema.
 *
 * @note The computed box is purely image-domain based. It assumes that the
 * tree exposes a valid original image domain through `rows` and `cols`.
 */
class BoundingBoxComputer : public AttributeComputer {
public:
    /**
     * @brief Returns the full family of bounding-box descriptors produced by
     * this computer.
     */
    std::vector<Attribute> attributes() const override {
        return {BOX_WIDTH, BOX_HEIGHT, RECTANGULARITY, RATIO_WH,BOX_COL_MIN, BOX_COL_MAX, BOX_ROW_MIN, BOX_ROW_MAX,DIAGONAL_LENGTH};
    }

    /**
     * @brief Declares the dependencies required by the derived descriptors.
     */
    std::vector<AttributeOrGroup> requiredAttributes() const override {
        return {AREA};
    }

    /**
     * @brief Computes the requested bounding-box descriptors.
     */
    void compute(MorphologicalTree& tree, const AltitudeBuffer*, std::span<float> buffer, const AttributeNames& attrNames, std::span<const Attribute> requestedAttributes, std::span<const DependencySource> dependencySources) const override {
        if (PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing BOUNDING_BOX group" << std::endl;

        auto indexOfWidth  = [&](NodeId idx) { return attrNames.linearIndex(idx, BOX_WIDTH); };
        auto indexOfHeight = [&](NodeId idx) { return attrNames.linearIndex(idx, BOX_HEIGHT); };
        auto indexOfRectangularity = [&](NodeId idx) { return attrNames.linearIndex(idx, RECTANGULARITY); };
        auto indexOfRatioWH = [&](NodeId idx) { return attrNames.linearIndex(idx, RATIO_WH); };
        auto indexOfColMin = [&](NodeId idx) { return attrNames.linearIndex(idx, BOX_COL_MIN); };
        auto indexOfColMax = [&](NodeId idx) { return attrNames.linearIndex(idx, BOX_COL_MAX); };
        auto indexOfRowMin = [&](NodeId idx) { return attrNames.linearIndex(idx, BOX_ROW_MIN); };
        auto indexOfRowMax = [&](NodeId idx) { return attrNames.linearIndex(idx, BOX_ROW_MAX); };
        auto indexOfDiagonalLength = [&](NodeId idx) { return attrNames.linearIndex(idx, DIAGONAL_LENGTH); };

        bool computeWidth  = std::find(requestedAttributes.begin(), requestedAttributes.end(), BOX_WIDTH)  != requestedAttributes.end();
        bool computeHeight = std::find(requestedAttributes.begin(), requestedAttributes.end(), BOX_HEIGHT) != requestedAttributes.end();
        bool computeRectangularity = std::find(requestedAttributes.begin(), requestedAttributes.end(), RECTANGULARITY) != requestedAttributes.end();
        bool computeRatioWH = std::find(requestedAttributes.begin(), requestedAttributes.end(), RATIO_WH) != requestedAttributes.end();
        bool computeColMin = std::find(requestedAttributes.begin(), requestedAttributes.end(), BOX_COL_MIN) != requestedAttributes.end();
        bool computeColMax = std::find(requestedAttributes.begin(), requestedAttributes.end(), BOX_COL_MAX) != requestedAttributes.end();
        bool computeRowMin = std::find(requestedAttributes.begin(), requestedAttributes.end(), BOX_ROW_MIN) != requestedAttributes.end();
        bool computeRowMax = std::find(requestedAttributes.begin(), requestedAttributes.end(), BOX_ROW_MAX) != requestedAttributes.end();
        bool computeDiagonalLength = std::find(requestedAttributes.begin(), requestedAttributes.end(), DIAGONAL_LENGTH) != requestedAttributes.end();
        const auto& dependencyArea = dependencySources[0];
        auto indexOfArea = [&](NodeId idx) { return dependencyArea.attrNames->linearIndex(idx, AREA); };

        int n = tree.getNumInternalNodeSlots();
        int numCols = tree.getNumColsOfImage();
        int numRows = tree.getNumRowsOfImage();

        std::vector<int> xmin(n, numCols);
        std::vector<int> xmax(n, 0);
        std::vector<int> ymin(n, numRows);
        std::vector<int> ymax(n, 0);

        AttributeComputedIncrementally::traversePostOrder(
            tree,
            tree.getRoot(),
            [&](NodeId nodeId) {
                const NodeId idx = detail::boundingBoxSlotOf(tree, nodeId);
                xmin[idx] = numCols;
                xmax[idx] = 0;
                ymin[idx] = numRows;
                ymax[idx] = 0;

                for (int p : tree.getProperParts(nodeId)) {
                    auto [y, x] = ImageUtils::to2D(p, numCols);
                    xmin[idx] = std::min(xmin[idx], x);
                    xmax[idx] = std::max(xmax[idx], x);
                    ymin[idx] = std::min(ymin[idx], y);
                    ymax[idx] = std::max(ymax[idx], y);
                }
            },
            [&](NodeId parentNodeId, NodeId childNodeId) {
                const NodeId pid = detail::boundingBoxSlotOf(tree, parentNodeId);
                const NodeId cid = detail::boundingBoxSlotOf(tree, childNodeId);
                xmin[pid] = std::min(xmin[pid], xmin[cid]);
                xmax[pid] = std::max(xmax[pid], xmax[cid]);
                ymin[pid] = std::min(ymin[pid], ymin[cid]);
                ymax[pid] = std::max(ymax[pid], ymax[cid]);
            },
            [&](NodeId nodeId) {
                const NodeId idx = detail::boundingBoxSlotOf(tree, nodeId);
                if(computeWidth)
                    buffer[indexOfWidth(idx)]  = xmax[idx] - xmin[idx] + 1;
                if(computeHeight)
                    buffer[indexOfHeight(idx)] = ymax[idx] - ymin[idx] + 1;

                if(computeRectangularity) {
                    float area = dependencyArea.buffer[indexOfArea(idx)];
                    float width = xmax[idx] - xmin[idx] + 1;
                    float height = ymax[idx] - ymin[idx] + 1;
                    float denom = width * height;
                    buffer[indexOfRectangularity(idx)] = (denom > 0.0f) ? (area / denom) : 0.0f;
                }
                if(computeRatioWH) {
                    float width  = xmax[idx] - xmin[idx] + 1;
                    float height = ymax[idx] - ymin[idx] + 1;
                    if (width > 0 && height > 0) {
                        buffer[indexOfRatioWH(idx)] = std::max(width, height) / std::min(width, height);
                    } else {
                        buffer[indexOfRatioWH(idx)] = 0.0f;
                    }
                }
                if(computeColMin)
                    buffer[indexOfColMin(idx)]  = xmin[idx];
                if(computeColMax)
                    buffer[indexOfColMax(idx)]  = xmax[idx];
                if(computeRowMin)
                    buffer[indexOfRowMin(idx)]  = ymin[idx];
                if(computeRowMax)
                    buffer[indexOfRowMax(idx)]  = ymax[idx];
                if(computeDiagonalLength) {
                    float width  = xmax[idx] - xmin[idx] + 1;
                    float height = ymax[idx] - ymin[idx] + 1;
                    buffer[indexOfDiagonalLength(idx)] = std::sqrt(width*width + height*height);
                }
            }
        );
    }
};

} // namespace mmcfilters
