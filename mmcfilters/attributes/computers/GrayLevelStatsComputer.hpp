#pragma once

#include "AttributeComputerDomain.hpp"
#include "AttributeComputerFamily.hpp"
#include "../detail/AttributeKernelSupport.hpp"
#include "../../trees/detail/TreeTraversalDetail.hpp"
#include "../../trees/detail/CommittedTreeAccess.hpp"
#include "../../trees/TreeAltitudeAlgorithms.hpp"
#include "../../utils/Altitude.hpp"
#include "../../utils/Contract.hpp"

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
namespace kernel {

/**
 * @brief Computes requested gray-level statistics from established altitude and dependency domains.
 * @param context Established tree, altitude span, output layout, and output buffer.
 * @param request Gray-level columns to materialize.
 * @param volumeDependency Optional established volume dependency.
 * @param areaDependency Optional established area dependency.
 */
template <std::floating_point Real, AltitudeValue T>
void computeGrayLevelStats(const AltitudeAttributeComputeContext<Real, T>& context, const GrayLevelStatsRequest& request,
                           const DependencySourceT<Real>* volumeDependency, const DependencySourceT<Real>* areaDependency) {
    if (!request.any()) {
        return;
    }

    const bool needsAggregateDependencies = request.needsAggregateDependencies();
    const int stride = context.attrNames.NUM_ATTRIBUTES;
    const auto offsetOf = [&](Attribute attribute) { return context.attrNames.indexMap.find(attribute)->second; };
    const int meanOffset = request.meanLevel ? offsetOf(MEAN_LEVEL) : 0;
    const int levelOffset = request.level ? offsetOf(LEVEL) : 0;
    const int varianceOffset = request.varianceLevel ? offsetOf(VARIANCE_LEVEL) : 0;
    const int grayHeightOffset = request.grayHeight ? offsetOf(GRAY_HEIGHT) : 0;
    auto indexOfMean = [&](NodeId node) { return static_cast<std::size_t>(node * stride + meanOffset); };
    auto indexOfLevel = [&](NodeId node) { return static_cast<std::size_t>(node * stride + levelOffset); };
    auto indexOfVariance = [&](NodeId node) { return static_cast<std::size_t>(node * stride + varianceOffset); };
    auto indexOfGrayHeight = [&](NodeId node) { return static_cast<std::size_t>(node * stride + grayHeightOffset); };
    const int volumeStride = volumeDependency != nullptr ? volumeDependency->attrNames->NUM_ATTRIBUTES : 0;
    const int volumeOffset = volumeDependency != nullptr ? volumeDependency->attrNames->indexMap.find(VOLUME)->second : 0;
    const int areaStride = areaDependency != nullptr ? areaDependency->attrNames->NUM_ATTRIBUTES : 0;
    const int areaOffset = areaDependency != nullptr ? areaDependency->attrNames->indexMap.find(AREA)->second : 0;
    auto indexOfVolume = [&](NodeId node) { return static_cast<std::size_t>(node * volumeStride + volumeOffset); };
    auto indexOfArea = [&](NodeId node) { return static_cast<std::size_t>(node * areaStride + areaOffset); };

    std::vector<double> sumGrayLevelSquare;
    if (request.varianceLevel) {
        sumGrayLevelSquare.assign(context.tree.getNumInternalNodeSlots(), 0.0);
    }
    std::vector<Real> subtreeMinAltitude;
    std::vector<Real> subtreeMaxAltitude;
    if (request.grayHeight) {
        subtreeMinAltitude.assign(context.tree.getNumInternalNodeSlots(), Real{0});
        subtreeMaxAltitude.assign(context.tree.getNumInternalNodeSlots(), Real{0});
    }

    ::mmcfilters::detail::kernel::traversePostOrder(
        context.tree, context.tree.getRoot(),
        [&](NodeId nodeId) {
            const NodeId node = detail::grayStatsSlotOf(context.tree, nodeId);
            const T nodeAltitude = context.altitude[static_cast<std::size_t>(nodeId)];
            const Real nodeAltitudeAsReal = static_cast<Real>(nodeAltitude);
            if (request.varianceLevel) {
                const double nodeAltitudeAsDouble = static_cast<double>(nodeAltitude);
                sumGrayLevelSquare[static_cast<std::size_t>(node)] =
                    static_cast<double>(::mmcfilters::detail::CommittedTreeAccess::numProperParts(context.tree, nodeId)) * nodeAltitudeAsDouble *
                    nodeAltitudeAsDouble;
            }
            if (request.level)
                context.buffer[indexOfLevel(node)] = nodeAltitudeAsReal;
            if (request.grayHeight) {
                subtreeMinAltitude[node] = nodeAltitudeAsReal;
                subtreeMaxAltitude[node] = nodeAltitudeAsReal;
            }
        },
        [&](NodeId parentNodeId, NodeId childNodeId) {
            const NodeId parent = detail::grayStatsSlotOf(context.tree, parentNodeId);
            const NodeId child = detail::grayStatsSlotOf(context.tree, childNodeId);
            if (request.varianceLevel)
                sumGrayLevelSquare[parent] += sumGrayLevelSquare[child];
            if (request.grayHeight) {
                subtreeMinAltitude[parent] = std::min(subtreeMinAltitude[parent], subtreeMinAltitude[child]);
                subtreeMaxAltitude[parent] = std::max(subtreeMaxAltitude[parent], subtreeMaxAltitude[child]);
            }
        },
        [&](NodeId nodeId) {
            const NodeId node = detail::grayStatsSlotOf(context.tree, nodeId);
            Real area = needsAggregateDependencies ? areaDependency->buffer[indexOfArea(node)] : Real{0};
            if (request.meanLevel)
                context.buffer[indexOfMean(node)] =
                    ::mmcfilters::attributes::numeric::safeDivide(volumeDependency->buffer[indexOfVolume(node)], area);
            if (request.varianceLevel) {
                Real meanGrayLevel = ::mmcfilters::attributes::numeric::safeDivide(volumeDependency->buffer[indexOfVolume(node)], area);
                double meanGrayLevelSquare =
                    ::mmcfilters::attributes::numeric::safeDivide(sumGrayLevelSquare[static_cast<std::size_t>(node)], static_cast<double>(area));
                Real var = static_cast<Real>(meanGrayLevelSquare - (static_cast<double>(meanGrayLevel) * static_cast<double>(meanGrayLevel)));
                context.buffer[indexOfVariance(node)] = ::mmcfilters::attributes::numeric::clampNonNegative(var);
            }
            if (request.grayHeight) {
                const Real nodeAltitude = static_cast<Real>(context.altitude[static_cast<std::size_t>(nodeId)]);
                context.buffer[indexOfGrayHeight(node)] =
                    std::max(std::abs(nodeAltitude - subtreeMinAltitude[node]), std::abs(subtreeMaxAltitude[node] - nodeAltitude));
            }
        });
}

} // namespace kernel

template <std::floating_point Real>
inline const DependencySourceT<Real>* findGrayDependency(std::span<const DependencySourceT<Real>> sources, Attribute attribute) {
    for (const DependencySourceT<Real>& source : sources) {
        if (source.attrNames->contains(attribute)) {
            return &source;
        }
    }
    return nullptr;
}

template <std::floating_point Real, AltitudeValue T>
inline void validateGrayLevelStatsContext(const AltitudeAttributeComputeContext<Real, T>& context, const GrayLevelStatsRequest& request) {
    requireAttributeBufferShape(context.tree, context.buffer, context.attrNames);
    TreeAltitudeAlgorithms::validateAltitudeBufferShape(context.tree, context.altitude);
    const auto requireColumn = [&](bool requested, Attribute attribute, const char* name) {
        if (requested && !context.attrNames.contains(attribute)) {
            throw std::invalid_argument(std::string(name) + " computation requires a matching output column.");
        }
    };
    requireColumn(request.level, LEVEL, "LEVEL");
    requireColumn(request.meanLevel, MEAN_LEVEL, "MEAN_LEVEL");
    requireColumn(request.varianceLevel, VARIANCE_LEVEL, "VARIANCE_LEVEL");
    requireColumn(request.grayHeight, GRAY_HEIGHT, "GRAY_HEIGHT");
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
        const detail::GrayLevelStatsRequest request = detail::GrayLevelStatsRequest::from(context.requestedAttributes);
        MMCFILTERS_CONTRACT_CHECKED_ONLY(detail::validateGrayLevelStatsContext(context, request));
        const DependencySourceT<Real>* volumeDependency = nullptr;
        const DependencySourceT<Real>* areaDependency = nullptr;
        if (request.needsAggregateDependencies()) {
            if constexpr (contract::validationsEnabled) {
                volumeDependency = &context.dependencies.require(VOLUME);
                areaDependency = &context.dependencies.require(AREA);
            } else {
                volumeDependency = detail::findGrayDependency(context.dependencySources, VOLUME);
                areaDependency = detail::findGrayDependency(context.dependencySources, AREA);
            }
        }
        detail::kernel::computeGrayLevelStats(context, request, volumeDependency, areaDependency);
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
                (request.level || request.meanLevel)
                    ? static_cast<Real>(::mmcfilters::detail::kernel::unitAltitude(context.tree, context.altitude, properPart))
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
