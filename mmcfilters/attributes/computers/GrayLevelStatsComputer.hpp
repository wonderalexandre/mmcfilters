#pragma once

#include "AttributeComputerDomain.hpp"
#include "AttributeComputerFamily.hpp"
#include "../detail/AttributeKernelSupport.hpp"
#include "../../trees/detail/TreeTraversalDetail.hpp"
#include "../../trees/TreeAltitudeAlgorithms.hpp"
#include "../../utils/Altitude.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <span>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace mmcfilters::attributes::computers {

namespace detail {
/**
 * @brief Returns the gray-statistics buffer slot for a tree node.
 *
 * @param nodeId Identifier of the node used by the operation.
 * @return Dense buffer slot for the node gray-level statistics.
 */
inline NodeId grayStatsSlotOf(const MorphologicalTree&, NodeId nodeId) noexcept { return nodeId; }

/** @brief Describes the gray-level statistics requested from the attribute pipeline. */
struct GrayLevelStatsRequest {
    /** @brief Indicates whether the gray-level sum was requested. */
    bool level = false;
    /** @brief Indicates whether the mean gray level was requested. */
    bool meanLevel = false;
    /** @brief Indicates whether the gray-level variance was requested. */
    bool varianceLevel = false;
    /** @brief Indicates whether the gray-height attribute was requested. */
    bool grayHeight = false;

    /**
     * @brief Tests whether any requested feature is enabled.
     *
     * @return True when any requested feature is enabled; otherwise false.
     */
    [[nodiscard]] bool any() const noexcept { return level || meanLevel || varianceLevel || grayHeight; }

    /**
     * @brief Tests whether aggregate dependencies holds.
     *
     * @return True when aggregate dependencies; otherwise false.
     */
    [[nodiscard]] bool needsAggregateDependencies() const noexcept { return meanLevel || varianceLevel; }

    /**
     * @brief Builds a request descriptor from the requested attributes.
     *
     * @param requestedAttributes Requested attribute subset.
     * @return Resulting request descriptor from the requested attributes.
     */
    [[nodiscard]] static GrayLevelStatsRequest from(std::span<const Attribute> requestedAttributes) {
        return {.level = containsGrayStatsAttribute(requestedAttributes, LEVEL),
                .meanLevel = containsGrayStatsAttribute(requestedAttributes, MEAN_LEVEL),
                .varianceLevel = containsGrayStatsAttribute(requestedAttributes, VARIANCE_LEVEL),
                .grayHeight = containsGrayStatsAttribute(requestedAttributes, GRAY_HEIGHT)};
    }

  private:
    /**
     * @brief Tests whether gray stats attribute holds.
     *
     * @param requestedAttributes Requested attribute subset.
     * @param attribute Attribute requested by the operation.
     * @return True when gray stats attribute; otherwise false.
     */
    [[nodiscard]] static bool containsGrayStatsAttribute(std::span<const Attribute> requestedAttributes, Attribute attribute) {
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
template <std::floating_point Real, AltitudeValue T>
void computeGrayLevelStatsAttributeKernel(const MorphologicalTree& tree, std::span<const T> altitude, std::span<Real> buffer, const AttributeNames& attrNames,
                                          std::span<const Attribute> requestedAttributes, std::span<const DependencySourceT<Real>> dependencySources) {
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
    const DependencySourceT<Real>* dependencyVol = needsAggregateDependencies ? &dependencies.require(VOLUME) : nullptr;
    const DependencySourceT<Real>* dependencyArea = needsAggregateDependencies ? &dependencies.require(AREA) : nullptr;
    auto indexOfVol = [&](NodeId idx) { return dependencyVol->attrNames->linearIndex(idx, VOLUME); };
    auto indexOfArea = [&](NodeId idx) { return dependencyArea->attrNames->linearIndex(idx, AREA); };

    std::vector<double> sumGrayLevelSquare;
    if (request.varianceLevel) {
        sumGrayLevelSquare.assign(tree.getNumInternalNodeSlots(), 0.0);
    }
    std::vector<Real> subtreeMinAltitude;
    std::vector<Real> subtreeMaxAltitude;
    if (request.grayHeight) {
        subtreeMinAltitude.assign(tree.getNumInternalNodeSlots(), Real{0});
        subtreeMaxAltitude.assign(tree.getNumInternalNodeSlots(), Real{0});
    }

    ::mmcfilters::detail::traversePostOrder(
        tree, tree.getRoot(),
        [&](NodeId nodeId) {
            const NodeId node = detail::grayStatsSlotOf(tree, nodeId);
            const T nodeAltitude = TreeAltitudeAlgorithms::getAltitude(altitude, nodeId);
            const Real nodeAltitudeAsReal = static_cast<Real>(nodeAltitude);
            if (request.varianceLevel) {
                const double nodeAltitudeAsDouble = static_cast<double>(nodeAltitude);
                sumGrayLevelSquare[node] = static_cast<double>(tree.getNumProperParts(nodeId)) * nodeAltitudeAsDouble * nodeAltitudeAsDouble;
            }
            if (request.level)
                buffer[indexOfLevel(node)] = nodeAltitudeAsReal;
            if (request.grayHeight) {
                subtreeMinAltitude[node] = nodeAltitudeAsReal;
                subtreeMaxAltitude[node] = nodeAltitudeAsReal;
            }
        },
        [&](NodeId parentNodeId, NodeId childNodeId) {
            const NodeId parent = detail::grayStatsSlotOf(tree, parentNodeId);
            const NodeId child = detail::grayStatsSlotOf(tree, childNodeId);
            if (request.varianceLevel)
                sumGrayLevelSquare[parent] += sumGrayLevelSquare[child];
            if (request.grayHeight) {
                subtreeMinAltitude[parent] = std::min(subtreeMinAltitude[parent], subtreeMinAltitude[child]);
                subtreeMaxAltitude[parent] = std::max(subtreeMaxAltitude[parent], subtreeMaxAltitude[child]);
            }
        },
        [&](NodeId nodeId) {
            const NodeId node = detail::grayStatsSlotOf(tree, nodeId);
            Real area = needsAggregateDependencies ? dependencyArea->buffer[indexOfArea(node)] : Real{0};
            if (request.meanLevel)
                buffer[indexOfMean(node)] = ::mmcfilters::attributes::numeric::safeDivide(dependencyVol->buffer[indexOfVol(node)], area);
            if (request.varianceLevel) {
                Real meanGrayLevel = ::mmcfilters::attributes::numeric::safeDivide(dependencyVol->buffer[indexOfVol(node)], area);
                double meanGrayLevelSquare = ::mmcfilters::attributes::numeric::safeDivide(sumGrayLevelSquare[node], static_cast<double>(area));
                Real var = static_cast<Real>(meanGrayLevelSquare - (static_cast<double>(meanGrayLevel) * static_cast<double>(meanGrayLevel)));
                buffer[indexOfVariance(node)] = ::mmcfilters::attributes::numeric::clampNonNegative(var);
            }
            if (request.grayHeight) {
                const Real nodeAltitude = static_cast<Real>(TreeAltitudeAlgorithms::getAltitude(altitude, nodeId));
                buffer[indexOfGrayHeight(node)] =
                    std::max(std::abs(nodeAltitude - subtreeMinAltitude[node]), std::abs(subtreeMaxAltitude[node] - nodeAltitude));
            }
        });
}

} // namespace detail

/**
 * @brief Computes grey-level statistics derived from subtree aggregation.
 *
 * @details
 * This computer provides four complementary descriptors:
 * - `LEVEL`: the altitude of the node itself;
 * - `MEAN_LEVEL`: the average grey level over the full subtree support;
 * - `VARIANCE_LEVEL`: the grey-level variance over the full subtree support;
 * - `GRAY_HEIGHT`: the maximum absolute altitude difference between the node
 *   and any node in its subtree.
 *
 * `MEAN_LEVEL` is obtained from the ratio `VOLUME / AREA`, while
 * `VARIANCE_LEVEL` additionally requires the accumulated sum of squared grey
 * levels over the subtree support. `GRAY_HEIGHT` propagates both the minimum
 * and maximum subtree altitudes and measures the farthest endpoint from the
 * current node.
 *
 * This definition is independent of a descriptive tree family and remains
 * meaningful when parent/child altitudes have no global monotonic order. For
 * monotone max-trees and min-trees it reduces to the traditional one-sided
 * grey-height definition.
 */
class GrayLevelStatsComputer {
  public:
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "gray-level-stats";

    /// Stable family id used by the scheduler.
    static constexpr AttributeComputerFamily family = AttributeComputerFamily::GrayLevelStats;

    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Altitude;

    /**
     * @brief Canonical list of grey-level descriptors produced by this computer.
     */
    inline static constexpr std::array<Attribute, 4> producedAttributes{LEVEL, MEAN_LEVEL, VARIANCE_LEVEL, GRAY_HEIGHT};

    /**
     * @brief Computes the requested grey-level statistics.
     *
     * @details
     * Requires a typed altitude span in dense internal-node order. `LEVEL` and
     * `GRAY_HEIGHT` use altitude directly. `MEAN_LEVEL` and `VARIANCE_LEVEL`
     * additionally require `VOLUME` and `AREA` dependencies available through
     * `context.dependencies`.
     *
     * @param context Operation context or diagnostic label.
     */
    template <std::floating_point Real, AltitudeValue T> static void compute(const AltitudeAttributeComputeContext<Real, T>& context) {
        requireAttributeBufferShape(context.tree, context.buffer, context.attrNames);
        detail::computeGrayLevelStatsAttributeKernel(context.tree, context.altitude, context.buffer, context.attrNames, context.requestedAttributes,
                                                     context.dependencySources);
    }

    /**
     * @brief Materializes grey-level statistics for one-pixel unit supports.
     *
     * Unit `LEVEL` and `MEAN_LEVEL` equal the owner-node altitude. Unit
     * `VARIANCE_LEVEL` and `GRAY_HEIGHT` are zero.
     *
     * @param context Operation context or diagnostic label.
     */
    template <std::floating_point Real, AltitudeValue T> static void computeUnitRows(const AltitudeUnitAttributeComputeContext<Real, T>& context) {
        requireUnitAttributeBufferShape(context.tree, context.unitProperParts, context.buffer, context.attrNames);
        TreeAltitudeAlgorithms::validateAltitudeBufferShape(context.tree, context.altitude);

        const detail::GrayLevelStatsRequest request = detail::GrayLevelStatsRequest::from(context.requestedAttributes);
        if (!request.any()) {
            return;
        }

        for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(context.unitProperParts.size()); ++leafIndex) {
            const NodeId properPart = context.unitProperParts[static_cast<size_t>(leafIndex)];
            const Real altitudeValue =
                (request.level || request.meanLevel) ? static_cast<Real>(unitAltitude(context.tree, context.altitude, properPart)) : Real{0};

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
