#pragma once

#include "../detail/AttributeKernelSupport.hpp"
#include "../../trees/detail/TreeTraversalDetail.hpp"
#include "../../trees/TreeAltitudeAlgorithms.hpp"
#include "../../utils/Altitude.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <span>
#include <stdexcept>
#include <vector>

namespace mmcfilters::attributes::computers {

namespace detail {
inline NodeId grayStatsSlotOf(const MorphologicalTree&, NodeId nodeId) noexcept {
    return nodeId;
}

inline constexpr std::array<Attribute, 4> GRAY_LEVEL_STATS_ATTRIBUTES{
    LEVEL,
    MEAN_LEVEL,
    VARIANCE_LEVEL,
    GRAY_HEIGHT};

struct GrayLevelStatsRequest {
    bool level = false;
    bool meanLevel = false;
    bool varianceLevel = false;
    bool grayHeight = false;

    [[nodiscard]] bool any() const noexcept {
        return level || meanLevel || varianceLevel || grayHeight;
    }

    [[nodiscard]] bool needsAggregateDependencies() const noexcept {
        return meanLevel || varianceLevel;
    }

    [[nodiscard]] static GrayLevelStatsRequest from(std::span<const Attribute> requestedAttributes) {
        return {
            .level = containsGrayStatsAttribute(requestedAttributes, LEVEL),
            .meanLevel = containsGrayStatsAttribute(requestedAttributes, MEAN_LEVEL),
            .varianceLevel = containsGrayStatsAttribute(requestedAttributes, VARIANCE_LEVEL),
            .grayHeight = containsGrayStatsAttribute(requestedAttributes, GRAY_HEIGHT)};
    }

private:
    [[nodiscard]] static bool containsGrayStatsAttribute(
        std::span<const Attribute> requestedAttributes,
        Attribute attribute)
    {
        return std::find(requestedAttributes.begin(), requestedAttributes.end(), attribute) != requestedAttributes.end();
    }
};

/**
 * @brief Typed kernel for grey-level statistics.
 *
 * `LEVEL` and `GRAY_HEIGHT` use only the altitude span. `MEAN_LEVEL` and
 * `VARIANCE_LEVEL` require dependency source `0` containing `VOLUME` and
 * source `1` containing `AREA`.
 *
 * @param tree Tree whose dense internal node ids index `buffer`.
 * @param altitude Dense altitude span indexed by internal node id.
 * @param buffer Flat internal-node output buffer.
 * @param attrNames Layout containing requested grey-level attributes.
 * @param requestedAttributes Requested subset of grey-level descriptors.
 * @param dependencySources Required dependency buffers for aggregate statistics.
 */
template<std::floating_point Real, AltitudeValue T>
void computeGrayLevelStatsAttributeKernel(
    const MorphologicalTree& tree,
    std::span<const T> altitude,
    std::span<Real> buffer,
    const AttributeNames& attrNames,
    std::span<const Attribute> requestedAttributes,
    std::span<const DependencySourceT<Real>> dependencySources)
{
    TreeAltitudeAlgorithms::validateAltitudeBufferShape(tree, altitude);

    const GrayLevelStatsRequest request = GrayLevelStatsRequest::from(requestedAttributes);
    if (!request.any()) {
        return;
    }

    const bool needsAggregateDependencies = request.needsAggregateDependencies();

    auto indexOfMean = [&](NodeId idx) { return attrNames.linearIndex(idx, MEAN_LEVEL); };
    auto indexOfLevel = [&](NodeId idx) { return attrNames.linearIndex(idx, LEVEL); };
    auto indexOfVariance = [&](NodeId idx) { return attrNames.linearIndex(idx, VARIANCE_LEVEL); };
    auto indexOfGrayHeight = [&](NodeId idx) { return attrNames.linearIndex(idx, GRAY_HEIGHT); };

    const DependencyResolver<Real> dependencies{dependencySources};
    const DependencySourceT<Real>* dependencyVol = needsAggregateDependencies
        ? &dependencies.require(VOLUME)
        : nullptr;
    const DependencySourceT<Real>* dependencyArea = needsAggregateDependencies
        ? &dependencies.require(AREA)
        : nullptr;
    auto indexOfVol = [&](NodeId idx) { return dependencyVol->attrNames->linearIndex(idx, VOLUME); };
    auto indexOfArea = [&](NodeId idx) { return dependencyArea->attrNames->linearIndex(idx, AREA); };

    std::vector<double> sumGrayLevelSquare;
    if (request.varianceLevel) {
        sumGrayLevelSquare.assign(tree.getNumInternalNodeSlots(), 0.0);
    }

    ::mmcfilters::detail::traversePostOrder(
        tree,
        tree.getRoot(),
        [&](NodeId nodeId) {
            const NodeId node = detail::grayStatsSlotOf(tree, nodeId);
            const T nodeAltitude = TreeAltitudeAlgorithms::getAltitude(altitude, nodeId);
            const Real nodeAltitudeAsReal = static_cast<Real>(nodeAltitude);
            if (request.varianceLevel) {
                const double nodeAltitudeAsDouble = static_cast<double>(nodeAltitude);
                sumGrayLevelSquare[node] =
                    static_cast<double>(tree.getNumProperParts(nodeId)) *
                    nodeAltitudeAsDouble *
                    nodeAltitudeAsDouble;
            }
            if (request.level)
                buffer[indexOfLevel(node)] = nodeAltitudeAsReal;
            if (request.grayHeight)
                buffer[indexOfGrayHeight(node)] = nodeAltitudeAsReal;
        },
        [&](NodeId parentNodeId, NodeId childNodeId) {
            const NodeId parent = detail::grayStatsSlotOf(tree, parentNodeId);
            const NodeId child = detail::grayStatsSlotOf(tree, childNodeId);
            if (request.varianceLevel)
                sumGrayLevelSquare[parent] += sumGrayLevelSquare[child];
            if (request.grayHeight) {
                Real childValue = buffer[indexOfGrayHeight(child)];
                Real& parentValue = buffer[indexOfGrayHeight(parent)];
                if (tree.getTreeType() == MorphologicalTreeKind::MAX_TREE)
                    parentValue = std::max(parentValue, childValue);
                else
                    parentValue = std::min(parentValue, childValue);
            }
        },
        [&](NodeId nodeId) {
            const NodeId node = detail::grayStatsSlotOf(tree, nodeId);
            Real area = needsAggregateDependencies ? dependencyArea->buffer[indexOfArea(node)] : Real{0};
            if (request.meanLevel)
                buffer[indexOfMean(node)] =
                    ::mmcfilters::attributes::numeric::safeDivide(dependencyVol->buffer[indexOfVol(node)], area);
            if (request.varianceLevel) {
                Real meanGrayLevel =
                    ::mmcfilters::attributes::numeric::safeDivide(dependencyVol->buffer[indexOfVol(node)], area);
                double meanGrayLevelSquare =
                    ::mmcfilters::attributes::numeric::safeDivide(sumGrayLevelSquare[node], static_cast<double>(area));
                Real var = static_cast<Real>(meanGrayLevelSquare - (static_cast<double>(meanGrayLevel) * static_cast<double>(meanGrayLevel)));
                buffer[indexOfVariance(node)] = ::mmcfilters::attributes::numeric::clampNonNegative(var);
            }
        });

    if (request.grayHeight) {
        for (NodeId nodeId : tree.getPostOrderNodes()) {
            const NodeId node = detail::grayStatsSlotOf(tree, nodeId);
            if (tree.isLeaf(nodeId)) {
                buffer[indexOfGrayHeight(node)] = 0.0f;
            } else {
                buffer[indexOfGrayHeight(node)] =
                    std::abs(static_cast<Real>(TreeAltitudeAlgorithms::getAltitude(altitude, nodeId)) -
                             buffer[indexOfGrayHeight(node)]);
            }
        }
    }
}

}

/**
 * @brief Computes grey-level statistics derived from subtree aggregation.
 *
 * @details
 * This computer provides four complementary descriptors:
 * - `LEVEL`: the altitude of the node itself;
 * - `MEAN_LEVEL`: the average grey level over the full subtree support;
 * - `VARIANCE_LEVEL`: the grey-level variance over the full subtree support;
 * - `GRAY_HEIGHT`: the grey-level span between the node altitude and the
 *   most extreme descendant altitude, measured according to the tree type.
 *
 * `MEAN_LEVEL` is obtained from the ratio `VOLUME / AREA`, while
 * `VARIANCE_LEVEL` additionally requires the accumulated sum of squared grey
 * levels over the subtree support. `GRAY_HEIGHT` is computed in two phases:
 * the traversal first propagates the most extreme descendant altitude upward
 * and a second pass then converts that extremal value into a span measured
 * from the current node altitude.
 *
 * The semantics of `GRAY_HEIGHT` depend on the hierarchy polarity:
 * - in a max-tree, the relevant descendant is the one with maximum altitude;
 * - in a min-tree, the relevant descendant is the one with minimum altitude.
 */
class GrayLevelStatsComputer {
public:

    /**
     * @brief Returns the grey-level descriptors naturally produced together.
     */
    [[nodiscard]] std::vector<Attribute> attributes() const {
        return {detail::GRAY_LEVEL_STATS_ATTRIBUTES.begin(), detail::GRAY_LEVEL_STATS_ATTRIBUTES.end()};
    }

    /**
     * @brief Computes the requested grey-level statistics.
     *
     * @details
     * Requires a typed altitude span in dense internal-node order. `LEVEL` and
     * `GRAY_HEIGHT` use altitude directly. `MEAN_LEVEL` and `VARIANCE_LEVEL`
     * additionally require `VOLUME` and `AREA` dependencies available through
     * `context.dependencies`.
     */
    template <std::floating_point Real, AltitudeValue T>
    static void compute(const AltitudeAttributeComputeContext<Real, T>& context)
    {
        requireAttributeBufferShape(context.tree, context.buffer, context.attrNames);
        detail::computeGrayLevelStatsAttributeKernel(
            context.tree,
            context.altitude,
            context.buffer,
            context.attrNames,
            context.requestedAttributes,
            context.dependencySources);
    }

    /**
     * @brief Materializes grey-level statistics for one-pixel unit supports.
     *
     * Unit `LEVEL` and `MEAN_LEVEL` equal the owner-node altitude. Unit
     * `VARIANCE_LEVEL` and `GRAY_HEIGHT` are zero.
     */
    template <std::floating_point Real, AltitudeValue T>
    static void computeUnitRows(const AltitudeUnitAttributeComputeContext<Real, T>& context)
    {
        requireUnitAttributeBufferShape(context.tree, context.unitProperParts, context.buffer, context.attrNames);
        TreeAltitudeAlgorithms::validateAltitudeBufferShape(context.tree, context.altitude);

        const detail::GrayLevelStatsRequest request = detail::GrayLevelStatsRequest::from(context.requestedAttributes);
        if (!request.any()) {
            return;
        }

        for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(context.unitProperParts.size()); ++leafIndex) {
            const NodeId properPart = context.unitProperParts[static_cast<size_t>(leafIndex)];
            const Real altitudeValue =
                (request.level || request.meanLevel)
                    ? static_cast<Real>(unitAltitude(context.tree, context.altitude, properPart))
                    : Real{0};

            if (request.level) {
                context.buffer[context.attrNames.linearIndex(leafIndex, LEVEL)] = altitudeValue;
            }
            if (request.meanLevel) {
                context.buffer[context.attrNames.linearIndex(leafIndex, MEAN_LEVEL)] = altitudeValue;
            }
            if (request.varianceLevel) {
                context.buffer[context.attrNames.linearIndex(leafIndex, VARIANCE_LEVEL)] = Real{0};
            }
            if (request.grayHeight) {
                context.buffer[context.attrNames.linearIndex(leafIndex, GRAY_HEIGHT)] = Real{0};
            }
        }
    }

};

} // namespace mmcfilters::attributes::computers
