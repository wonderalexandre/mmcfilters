#pragma once

#include "AttributeComputerDomain.hpp"
#include "AttributeComputerFamily.hpp"
#include "../detail/AttributeKernelSupport.hpp"
#include "../../trees/detail/TreeTraversalDetail.hpp"
#include "../../trees/MorphologicalTree.hpp"
#include "../../utils/Image.hpp"

#include <algorithm>
#include <array>
#include <concepts>
#include <cmath>
#include <span>
#include <string_view>
#include <vector>

namespace mmcfilters::attributes::computers {

namespace detail {
/**
 * @brief Returns the bounding-box buffer slot for a tree node.
 *
 * @param nodeId Identifier of the node used by the operation.
 * @return Dense buffer slot for the node bounding box.
 */
inline NodeId boundingBoxSlotOf(const MorphologicalTree&, NodeId nodeId) noexcept { return nodeId; }

/** @brief Describes the bounding-box attributes requested from the attribute pipeline. */
struct BoundingBoxRequest {
    /** @brief Indicates whether the width attribute was requested. */
    bool width = false;
    /** @brief Indicates whether the height attribute was requested. */
    bool height = false;
    /** @brief Indicates whether the rectangularity attribute was requested. */
    bool rectangularity = false;
    /** @brief Indicates whether the width-to-height ratio was requested. */
    bool ratioWH = false;
    /** @brief Indicates whether the minimum-column coordinate was requested. */
    bool colMin = false;
    /** @brief Indicates whether the maximum-column coordinate was requested. */
    bool colMax = false;
    /** @brief Indicates whether the minimum-row coordinate was requested. */
    bool rowMin = false;
    /** @brief Indicates whether the maximum-row coordinate was requested. */
    bool rowMax = false;
    /** @brief Indicates whether the bounding-box diagonal length was requested. */
    bool diagonalLength = false;

    /**
     * @brief Tests whether any requested feature is enabled.
     *
     * @return True when any requested feature is enabled; otherwise false.
     */
    [[nodiscard]] bool any() const noexcept { return width || height || rectangularity || ratioWH || colMin || colMax || rowMin || rowMax || diagonalLength; }

    /**
     * @brief Tests whether area dependency holds.
     *
     * @return True when area dependency; otherwise false.
     */
    [[nodiscard]] bool needsAreaDependency() const noexcept { return rectangularity; }

    /**
     * @brief Builds a request descriptor from the requested attributes.
     *
     * @param requestedAttributes Requested attribute subset.
     * @return Resulting request descriptor from the requested attributes.
     */
    [[nodiscard]] static BoundingBoxRequest from(std::span<const Attribute> requestedAttributes) {
        return {.width = contains(requestedAttributes, BOX_WIDTH),
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
    /**
     * @brief Tests whether contains holds.
     *
     * @param requestedAttributes Requested attribute subset.
     * @param attribute Attribute requested by the operation.
     * @return True when contains; otherwise false.
     */
    [[nodiscard]] static bool contains(std::span<const Attribute> requestedAttributes, Attribute attribute) {
        return std::find(requestedAttributes.begin(), requestedAttributes.end(), attribute) != requestedAttributes.end();
    }
};
} // namespace detail

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
class BoundingBoxComputer {
  public:
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "bounding-box";

    /// Stable family id used by the scheduler.
    static constexpr AttributeComputerFamily family = AttributeComputerFamily::BoundingBox;

    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Topology;

    /**
     * @brief Canonical list of bounding-box descriptors produced by this computer.
     */
    inline static constexpr std::array<Attribute, 9> producedAttributes{BOX_WIDTH,   BOX_HEIGHT,  DIAGONAL_LENGTH, RECTANGULARITY, RATIO_WH,
                                                                        BOX_COL_MIN, BOX_COL_MAX, BOX_ROW_MIN,     BOX_ROW_MAX};

    /**
     * @brief Computes the requested bounding-box descriptors.
     *
     * @details
     * The output buffer is indexed by dense internal node id and interpreted by
     * `context.attrNames`. `context.requestedAttributes` selects which columns
     * are written. `RECTANGULARITY` requires an `AREA` dependency available
     * through `context.dependencies`; all other descriptors are computed
     * directly from row-major proper-part coordinates.
     *
     * @param context Operation context or diagnostic label.
     */
    template <std::floating_point Real> static void compute(const AttributeComputeContext<Real>& context) {
        const MorphologicalTree& tree = context.tree;
        std::span<Real> buffer = context.buffer;
        const AttributeNames& attrNames = context.attrNames;
        std::span<const Attribute> requestedAttributes = context.requestedAttributes;

        requireAttributeBufferShape(tree, buffer, attrNames);

        auto indexOfWidth = [&](NodeId idx) { return attrNames.linearIndex(idx, BOX_WIDTH); };
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

        const DependencySourceT<Real>* dependencyArea = request.needsAreaDependency() ? &context.dependencies.require(AREA) : nullptr;
        auto indexOfArea = [&](NodeId idx) { return dependencyArea->attrNames->linearIndex(idx, AREA); };

        int n = tree.getNumInternalNodeSlots();
        int numCols = tree.getNumColsOfGridDomain2D();
        int numRows = tree.getNumRowsOfGridDomain2D();

        std::vector<int> xmin(n, numCols);
        std::vector<int> xmax(n, 0);
        std::vector<int> ymin(n, numRows);
        std::vector<int> ymax(n, 0);

        ::mmcfilters::detail::traversePostOrder(
            tree, tree.getRoot(),
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
                if (request.width)
                    buffer[indexOfWidth(idx)] = static_cast<Real>(xmax[idx] - xmin[idx] + 1);
                if (request.height)
                    buffer[indexOfHeight(idx)] = static_cast<Real>(ymax[idx] - ymin[idx] + 1);

                if (request.rectangularity) {
                    Real area = dependencyArea->buffer[indexOfArea(idx)];
                    Real width = static_cast<Real>(xmax[idx] - xmin[idx] + 1);
                    Real height = static_cast<Real>(ymax[idx] - ymin[idx] + 1);
                    Real denom = width * height;
                    buffer[indexOfRectangularity(idx)] = ::mmcfilters::attributes::numeric::safeDivide(area, denom);
                }
                if (request.ratioWH) {
                    Real width = static_cast<Real>(xmax[idx] - xmin[idx] + 1);
                    Real height = static_cast<Real>(ymax[idx] - ymin[idx] + 1);
                    buffer[indexOfRatioWH(idx)] = ::mmcfilters::attributes::numeric::safeDivide(std::max(width, height), std::min(width, height));
                }
                if (request.colMin)
                    buffer[indexOfColMin(idx)] = static_cast<Real>(xmin[idx]);
                if (request.colMax)
                    buffer[indexOfColMax(idx)] = static_cast<Real>(xmax[idx]);
                if (request.rowMin)
                    buffer[indexOfRowMin(idx)] = static_cast<Real>(ymin[idx]);
                if (request.rowMax)
                    buffer[indexOfRowMax(idx)] = static_cast<Real>(ymax[idx]);
                if (request.diagonalLength) {
                    Real width = static_cast<Real>(xmax[idx] - xmin[idx] + 1);
                    Real height = static_cast<Real>(ymax[idx] - ymin[idx] + 1);
                    buffer[indexOfDiagonalLength(idx)] = ::mmcfilters::attributes::numeric::safeSqrt(width * width + height * height);
                }
            });
    }

    /**
     * @brief Materializes bounding-box descriptors for one-pixel unit supports.
     *
     * A unit support has width/height `1`, rectangularity `1`, ratio `1`, and
     * min/max coordinates equal to the proper part coordinate.
     *
     * @param context Operation context or diagnostic label.
     */
    template <std::floating_point Real> static void computeUnitRows(const UnitAttributeComputeContext<Real>& context) {
        requireUnitAttributeBufferShape(context.tree, context.unitProperParts, context.buffer, context.attrNames);

        const detail::BoundingBoxRequest request = detail::BoundingBoxRequest::from(context.requestedAttributes);
        if (!request.any()) {
            return;
        }

        const int numCols = context.tree.getNumColsOfGridDomain2D();
        for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(context.unitProperParts.size()); ++leafIndex) {
            const NodeId properPart = context.unitProperParts[static_cast<size_t>(leafIndex)];
            const auto [row, col] = ImageUtils::to2D(properPart, numCols);
            if (request.width) {
                context.buffer[context.attrNames.linearIndex(leafIndex, BOX_WIDTH)] = Real{1};
            }
            if (request.height) {
                context.buffer[context.attrNames.linearIndex(leafIndex, BOX_HEIGHT)] = Real{1};
            }
            if (request.rectangularity) {
                context.buffer[context.attrNames.linearIndex(leafIndex, RECTANGULARITY)] = Real{1};
            }
            if (request.ratioWH) {
                context.buffer[context.attrNames.linearIndex(leafIndex, RATIO_WH)] = Real{1};
            }
            if (request.colMin) {
                context.buffer[context.attrNames.linearIndex(leafIndex, BOX_COL_MIN)] = static_cast<Real>(col);
            }
            if (request.colMax) {
                context.buffer[context.attrNames.linearIndex(leafIndex, BOX_COL_MAX)] = static_cast<Real>(col);
            }
            if (request.rowMin) {
                context.buffer[context.attrNames.linearIndex(leafIndex, BOX_ROW_MIN)] = static_cast<Real>(row);
            }
            if (request.rowMax) {
                context.buffer[context.attrNames.linearIndex(leafIndex, BOX_ROW_MAX)] = static_cast<Real>(row);
            }
            if (request.diagonalLength) {
                context.buffer[context.attrNames.linearIndex(leafIndex, DIAGONAL_LENGTH)] = std::sqrt(Real{2});
            }
        }
    }
};

} // namespace mmcfilters::attributes::computers
