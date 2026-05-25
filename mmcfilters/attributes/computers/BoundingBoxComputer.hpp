#pragma once

#include "../AttributeComputer.hpp"
#include "../../trees/detail/TreeTraversalDetail.hpp"
#include "../../trees/MorphologicalTree.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <span>
#include <vector>

namespace mmcfilters::attributes::computers {

namespace detail {
inline NodeId boundingBoxSlotOf(const MorphologicalTree&, NodeId nodeId) noexcept {
    return nodeId;
}

inline constexpr std::array<Attribute, 9> BOUNDING_BOX_ATTRIBUTES{
    BOX_WIDTH,
    BOX_HEIGHT,
    RECTANGULARITY,
    RATIO_WH,
    BOX_COL_MIN,
    BOX_COL_MAX,
    BOX_ROW_MIN,
    BOX_ROW_MAX,
    DIAGONAL_LENGTH};

struct BoundingBoxRequest {
    bool width = false;
    bool height = false;
    bool rectangularity = false;
    bool ratioWH = false;
    bool colMin = false;
    bool colMax = false;
    bool rowMin = false;
    bool rowMax = false;
    bool diagonalLength = false;

    [[nodiscard]] bool any() const noexcept {
        return width || height || rectangularity || ratioWH || colMin ||
               colMax || rowMin || rowMax || diagonalLength;
    }

    [[nodiscard]] bool needsAreaDependency() const noexcept {
        return rectangularity;
    }

    [[nodiscard]] static BoundingBoxRequest from(std::span<const Attribute> requestedAttributes) {
        return {
            .width = contains(requestedAttributes, BOX_WIDTH),
            .height = contains(requestedAttributes, BOX_HEIGHT),
            .rectangularity = contains(requestedAttributes, RECTANGULARITY),
            .ratioWH = contains(requestedAttributes, RATIO_WH),
            .colMin = contains(requestedAttributes, BOX_COL_MIN),
            .colMax = contains(requestedAttributes, BOX_COL_MAX),
            .rowMin = contains(requestedAttributes, BOX_ROW_MIN),
            .rowMax = contains(requestedAttributes, BOX_ROW_MAX),
            .diagonalLength = contains(requestedAttributes, DIAGONAL_LENGTH)};
    }

private:
    [[nodiscard]] static bool contains(std::span<const Attribute> requestedAttributes, Attribute attribute) {
        return std::find(requestedAttributes.begin(), requestedAttributes.end(), attribute) != requestedAttributes.end();
    }
};
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
    using AttributeComputer::compute;
    using AttributeComputer::computeUnitAttributes;

    /**
     * @brief Returns the full family of bounding-box descriptors produced by
     * this computer.
     */
    [[nodiscard]] std::vector<Attribute> attributes() const override {
        return {detail::BOUNDING_BOX_ATTRIBUTES.begin(), detail::BOUNDING_BOX_ATTRIBUTES.end()};
    }

    /**
     * @brief Computes the requested bounding-box descriptors.
     *
     * The output buffer is indexed by dense internal node id. `RECTANGULARITY`
     * requires dependency source `0` containing `AREA`; all other descriptors
     * are computed directly from row-major proper-part coordinates.
     */
    void compute(const MorphologicalTree& tree, AttributeAltitudeView, std::span<float> buffer, const AttributeNames& attrNames, std::span<const Attribute> requestedAttributes, std::span<const DependencySource> dependencySources) const override {
        requireAttributeBufferShape(tree, buffer, attrNames);

        auto indexOfWidth  = [&](NodeId idx) { return attrNames.linearIndex(idx, BOX_WIDTH); };
        auto indexOfHeight = [&](NodeId idx) { return attrNames.linearIndex(idx, BOX_HEIGHT); };
        auto indexOfRectangularity = [&](NodeId idx) { return attrNames.linearIndex(idx, RECTANGULARITY); };
        auto indexOfRatioWH = [&](NodeId idx) { return attrNames.linearIndex(idx, RATIO_WH); };
        auto indexOfColMin = [&](NodeId idx) { return attrNames.linearIndex(idx, BOX_COL_MIN); };
        auto indexOfColMax = [&](NodeId idx) { return attrNames.linearIndex(idx, BOX_COL_MAX); };
        auto indexOfRowMin = [&](NodeId idx) { return attrNames.linearIndex(idx, BOX_ROW_MIN); };
        auto indexOfRowMax = [&](NodeId idx) { return attrNames.linearIndex(idx, BOX_ROW_MAX); };
        auto indexOfDiagonalLength = [&](NodeId idx) { return attrNames.linearIndex(idx, DIAGONAL_LENGTH); };

        const detail::BoundingBoxRequest request = detail::BoundingBoxRequest::from(requestedAttributes);
        if (!request.any()) {
            return;
        }

        const DependencySource* dependencyArea = request.needsAreaDependency()
            ? &requireDependencySource(dependencySources, 0, AREA)
            : nullptr;
        auto indexOfArea = [&](NodeId idx) { return dependencyArea->attrNames->linearIndex(idx, AREA); };

        int n = tree.getNumInternalNodeSlots();
        int numCols = tree.getNumColsOfImage();
        int numRows = tree.getNumRowsOfImage();

        std::vector<int> xmin(n, numCols);
        std::vector<int> xmax(n, 0);
        std::vector<int> ymin(n, numRows);
        std::vector<int> ymax(n, 0);

        ::mmcfilters::detail::traversePostOrder(
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
                if(request.width)
                    buffer[indexOfWidth(idx)]  = xmax[idx] - xmin[idx] + 1;
                if(request.height)
                    buffer[indexOfHeight(idx)] = ymax[idx] - ymin[idx] + 1;

                if(request.rectangularity) {
                    float area = dependencyArea->buffer[indexOfArea(idx)];
                    float width = xmax[idx] - xmin[idx] + 1;
                    float height = ymax[idx] - ymin[idx] + 1;
                    float denom = width * height;
                    buffer[indexOfRectangularity(idx)] = (denom > 0.0f) ? (area / denom) : 0.0f;
                }
                if(request.ratioWH) {
                    float width  = xmax[idx] - xmin[idx] + 1;
                    float height = ymax[idx] - ymin[idx] + 1;
                    if (width > 0 && height > 0) {
                        buffer[indexOfRatioWH(idx)] = std::max(width, height) / std::min(width, height);
                    } else {
                        buffer[indexOfRatioWH(idx)] = 0.0f;
                    }
                }
                if(request.colMin)
                    buffer[indexOfColMin(idx)]  = xmin[idx];
                if(request.colMax)
                    buffer[indexOfColMax(idx)]  = xmax[idx];
                if(request.rowMin)
                    buffer[indexOfRowMin(idx)]  = ymin[idx];
                if(request.rowMax)
                    buffer[indexOfRowMax(idx)]  = ymax[idx];
                if(request.diagonalLength) {
                    float width  = xmax[idx] - xmin[idx] + 1;
                    float height = ymax[idx] - ymin[idx] + 1;
                    buffer[indexOfDiagonalLength(idx)] = std::sqrt(width*width + height*height);
                }
            }
        );
    }

    /**
     * @brief Materializes bounding-box descriptors for one-pixel unit supports.
     *
     * A unit support has width/height `1`, rectangularity `1`, ratio `1`, and
     * min/max coordinates equal to the proper part coordinate.
     */
    void computeUnitAttributes(
        const MorphologicalTree& tree,
        AttributeAltitudeView,
        std::span<const NodeId> unitProperParts,
        std::span<float> buffer,
        const AttributeNames& attrNames,
        std::span<const Attribute> requestedAttributes) const override
    {
        requireUnitAttributeBufferShape(tree, unitProperParts, buffer, attrNames);

        const detail::BoundingBoxRequest request = detail::BoundingBoxRequest::from(requestedAttributes);
        if (!request.any()) {
            return;
        }

        const int numCols = tree.getNumColsOfImage();
        for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(unitProperParts.size()); ++leafIndex) {
            const NodeId properPart = unitProperParts[static_cast<size_t>(leafIndex)];
            const auto [row, col] = ImageUtils::to2D(properPart, numCols);
            if (request.width) {
                buffer[attrNames.linearIndex(leafIndex, BOX_WIDTH)] = 1.0f;
            }
            if (request.height) {
                buffer[attrNames.linearIndex(leafIndex, BOX_HEIGHT)] = 1.0f;
            }
            if (request.rectangularity) {
                buffer[attrNames.linearIndex(leafIndex, RECTANGULARITY)] = 1.0f;
            }
            if (request.ratioWH) {
                buffer[attrNames.linearIndex(leafIndex, RATIO_WH)] = 1.0f;
            }
            if (request.colMin) {
                buffer[attrNames.linearIndex(leafIndex, BOX_COL_MIN)] = static_cast<float>(col);
            }
            if (request.colMax) {
                buffer[attrNames.linearIndex(leafIndex, BOX_COL_MAX)] = static_cast<float>(col);
            }
            if (request.rowMin) {
                buffer[attrNames.linearIndex(leafIndex, BOX_ROW_MIN)] = static_cast<float>(row);
            }
            if (request.rowMax) {
                buffer[attrNames.linearIndex(leafIndex, BOX_ROW_MAX)] = static_cast<float>(row);
            }
            if (request.diagonalLength) {
                buffer[attrNames.linearIndex(leafIndex, DIAGONAL_LENGTH)] = std::sqrt(2.0f);
            }
        }
    }
};

} // namespace mmcfilters::attributes::computers
