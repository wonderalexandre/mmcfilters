#pragma once

#include "AttributeComputerDomain.hpp"
#include "AttributeComputerFamily.hpp"
#include "../detail/AttributeKernelSupport.hpp"
#include "../../trees/detail/TreeTraversalDetail.hpp"
#include "../../trees/detail/CommittedTreeAccess.hpp"
#include "../../trees/MorphologicalTree.hpp"
#include "../../utils/Contract.hpp"
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

namespace kernel {

/**
 * @brief Computes requested bounding-box descriptors over an established tree.
 * @param context Established tree, output layout, and output buffer.
 * @param request Bounding-box columns to materialize.
 * @param areaDependency Optional established area dependency.
 */
template <std::floating_point Real>
inline void computeBoundingBox(const AttributeComputeContext<Real>& context, const BoundingBoxRequest& request,
                               const DependencySourceT<Real>* areaDependency) {
    if (!request.any()) {
        return;
    }

    const int stride = context.attrNames.NUM_ATTRIBUTES;
    const auto offsetOf = [&](Attribute attribute) { return context.attrNames.indexMap.find(attribute)->second; };
    const int widthOffset = request.width ? offsetOf(BOX_WIDTH) : 0;
    const int heightOffset = request.height ? offsetOf(BOX_HEIGHT) : 0;
    const int rectangularityOffset = request.rectangularity ? offsetOf(RECTANGULARITY) : 0;
    const int ratioOffset = request.ratioWH ? offsetOf(RATIO_WH) : 0;
    const int colMinOffset = request.colMin ? offsetOf(BOX_COL_MIN) : 0;
    const int colMaxOffset = request.colMax ? offsetOf(BOX_COL_MAX) : 0;
    const int rowMinOffset = request.rowMin ? offsetOf(BOX_ROW_MIN) : 0;
    const int rowMaxOffset = request.rowMax ? offsetOf(BOX_ROW_MAX) : 0;
    const int diagonalOffset = request.diagonalLength ? offsetOf(DIAGONAL_LENGTH) : 0;
    auto outputIndex = [&](NodeId node, int offset) { return static_cast<std::size_t>(node * stride + offset); };

    const int areaStride = areaDependency != nullptr ? areaDependency->attrNames->NUM_ATTRIBUTES : 0;
    const int areaOffset = areaDependency != nullptr ? areaDependency->attrNames->indexMap.find(AREA)->second : 0;
    auto areaIndex = [&](NodeId node) { return static_cast<std::size_t>(node * areaStride + areaOffset); };

    const int numNodes = context.tree.getNumInternalNodeSlots();
    const GridDomain2D& domain = ::mmcfilters::detail::CommittedTreeAccess::gridDomain2D(context.tree);
    std::vector<int> columnMin(static_cast<std::size_t>(numNodes), domain.cols);
    std::vector<int> columnMax(static_cast<std::size_t>(numNodes), 0);
    std::vector<int> rowMin(static_cast<std::size_t>(numNodes), domain.rows);
    std::vector<int> rowMax(static_cast<std::size_t>(numNodes), 0);

    ::mmcfilters::detail::kernel::traversePostOrder(
        context.tree, context.tree.getRoot(),
        [&](NodeId node) {
            for (int properPart : ::mmcfilters::detail::CommittedTreeAccess::properParts(context.tree, node)) {
                const int row = properPart / domain.cols;
                const int column = properPart % domain.cols;
                columnMin[static_cast<std::size_t>(node)] = std::min(columnMin[static_cast<std::size_t>(node)], column);
                columnMax[static_cast<std::size_t>(node)] = std::max(columnMax[static_cast<std::size_t>(node)], column);
                rowMin[static_cast<std::size_t>(node)] = std::min(rowMin[static_cast<std::size_t>(node)], row);
                rowMax[static_cast<std::size_t>(node)] = std::max(rowMax[static_cast<std::size_t>(node)], row);
            }
        },
        [&](NodeId parent, NodeId child) {
            columnMin[static_cast<std::size_t>(parent)] =
                std::min(columnMin[static_cast<std::size_t>(parent)], columnMin[static_cast<std::size_t>(child)]);
            columnMax[static_cast<std::size_t>(parent)] =
                std::max(columnMax[static_cast<std::size_t>(parent)], columnMax[static_cast<std::size_t>(child)]);
            rowMin[static_cast<std::size_t>(parent)] = std::min(rowMin[static_cast<std::size_t>(parent)], rowMin[static_cast<std::size_t>(child)]);
            rowMax[static_cast<std::size_t>(parent)] = std::max(rowMax[static_cast<std::size_t>(parent)], rowMax[static_cast<std::size_t>(child)]);
        },
        [&](NodeId node) {
            const Real width = static_cast<Real>(columnMax[static_cast<std::size_t>(node)] - columnMin[static_cast<std::size_t>(node)] + 1);
            const Real height = static_cast<Real>(rowMax[static_cast<std::size_t>(node)] - rowMin[static_cast<std::size_t>(node)] + 1);
            if (request.width)
                context.buffer[outputIndex(node, widthOffset)] = width;
            if (request.height)
                context.buffer[outputIndex(node, heightOffset)] = height;
            if (request.rectangularity)
                context.buffer[outputIndex(node, rectangularityOffset)] =
                    ::mmcfilters::attributes::numeric::safeDivide(areaDependency->buffer[areaIndex(node)], width * height);
            if (request.ratioWH)
                context.buffer[outputIndex(node, ratioOffset)] =
                    ::mmcfilters::attributes::numeric::safeDivide(std::max(width, height), std::min(width, height));
            if (request.colMin)
                context.buffer[outputIndex(node, colMinOffset)] = static_cast<Real>(columnMin[static_cast<std::size_t>(node)]);
            if (request.colMax)
                context.buffer[outputIndex(node, colMaxOffset)] = static_cast<Real>(columnMax[static_cast<std::size_t>(node)]);
            if (request.rowMin)
                context.buffer[outputIndex(node, rowMinOffset)] = static_cast<Real>(rowMin[static_cast<std::size_t>(node)]);
            if (request.rowMax)
                context.buffer[outputIndex(node, rowMaxOffset)] = static_cast<Real>(rowMax[static_cast<std::size_t>(node)]);
            if (request.diagonalLength)
                context.buffer[outputIndex(node, diagonalOffset)] = ::mmcfilters::attributes::numeric::safeSqrt(width * width + height * height);
        });
}

} // namespace kernel

template <std::floating_point Real>
inline void validateBoundingBoxContext(const AttributeComputeContext<Real>& context) {
    requireAttributeBufferShape(context.tree, context.buffer, context.attrNames);
    requireRequestedAttributeColumns(context);
    static_cast<void>(context.tree.requireGridDomain2D("BoundingBoxComputer"));
}
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
        const detail::BoundingBoxRequest request = detail::BoundingBoxRequest::from(context.requestedAttributes);
        MMCFILTERS_CONTRACT_CHECKED_ONLY(detail::validateBoundingBoxContext(context));

        const DependencySourceT<Real>* areaDependency = nullptr;
        if (request.needsAreaDependency()) {
            if constexpr (contract::validationsEnabled) {
                areaDependency = &context.dependencies.require(AREA);
            } else {
                areaDependency = ::mmcfilters::findDependencySource(context.dependencySources, AREA);
            }
        }
        detail::kernel::computeBoundingBox(context, request, areaDependency);
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
