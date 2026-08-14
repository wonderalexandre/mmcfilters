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
 * @param nodeId Dense internal node identifier.
 * @return Dense buffer slot for the node gray-level statistics.
 */
inline NodeId grayStatsSlotOf(const MorphologicalTree&, NodeId nodeId) noexcept { return nodeId; }

/** @brief Describes the gray-level statistics requested from the attribute pipeline. */
struct GrayLevelStatsRequest {
    /** @brief Indicates whether the mean gray level was requested. */
    bool meanGrayLevel = false;
    /** @brief Indicates whether the gray-level variance was requested. */
    bool grayLevelVariance = false;
    /** @brief Indicates whether the gray-level-height attribute was requested. */
    bool grayLevelHeight = false;

    /**
     * @brief Tests whether any requested feature is enabled.
     *
     * @return True when any requested feature is enabled; otherwise false.
     */
    [[nodiscard]] bool any() const noexcept { return meanGrayLevel || grayLevelVariance || grayLevelHeight; }

    /**
     * @brief Tests whether aggregate dependencies holds.
     *
     * @return True when aggregate dependencies; otherwise false.
     */
    [[nodiscard]] bool needsAggregateDependencies() const noexcept { return meanGrayLevel || grayLevelVariance; }

    /**
     * @brief Builds a request descriptor from the requested attributes.
     *
     * @param requestedAttributes Requested attribute subset.
     * @return Resulting request descriptor from the requested attributes.
     */
    [[nodiscard]] static GrayLevelStatsRequest from(std::span<const Attribute> requestedAttributes) {
        return {.meanGrayLevel = containsGrayStatsAttribute(requestedAttributes, MeanGrayLevel),
                .grayLevelVariance = containsGrayStatsAttribute(requestedAttributes, GrayLevelVariance),
                .grayLevelHeight = containsGrayStatsAttribute(requestedAttributes, GrayLevelHeight)};
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
 * `GrayLevelHeight` uses only the altitude span. `MeanGrayLevel` and
 * `GrayLevelVariance` require dependencies containing `VOLUME` and `AREA`.
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
    const int meanOffset = request.meanGrayLevel ? offsetOf(MeanGrayLevel) : 0;
    const int varianceOffset = request.grayLevelVariance ? offsetOf(GrayLevelVariance) : 0;
    const int grayHeightOffset = request.grayLevelHeight ? offsetOf(GrayLevelHeight) : 0;
    auto indexOfMean = [&](NodeId node) { return static_cast<std::size_t>(node * stride + meanOffset); };
    auto indexOfVariance = [&](NodeId node) { return static_cast<std::size_t>(node * stride + varianceOffset); };
    auto indexOfGrayHeight = [&](NodeId node) { return static_cast<std::size_t>(node * stride + grayHeightOffset); };
    const int volumeStride = volumeDependency != nullptr ? volumeDependency->attrNames->NUM_ATTRIBUTES : 0;
    const int volumeOffset = volumeDependency != nullptr ? volumeDependency->attrNames->indexMap.find(Volume)->second : 0;
    const int areaStride = areaDependency != nullptr ? areaDependency->attrNames->NUM_ATTRIBUTES : 0;
    const int areaOffset = areaDependency != nullptr ? areaDependency->attrNames->indexMap.find(Area)->second : 0;
    auto indexOfVolume = [&](NodeId node) { return static_cast<std::size_t>(node * volumeStride + volumeOffset); };
    auto indexOfArea = [&](NodeId node) { return static_cast<std::size_t>(node * areaStride + areaOffset); };

    std::vector<double> sumGrayLevelSquare;
    if (request.grayLevelVariance) {
        sumGrayLevelSquare.assign(context.tree.numInternalNodeSlots(), 0.0);
    }
    std::vector<Real> subtreeMinAltitude;
    std::vector<Real> subtreeMaxAltitude;
    if (request.grayLevelHeight) {
        subtreeMinAltitude.assign(context.tree.numInternalNodeSlots(), Real{0});
        subtreeMaxAltitude.assign(context.tree.numInternalNodeSlots(), Real{0});
    }

    ::mmcfilters::detail::kernel::traversePostOrder(
        context.tree, context.tree.root(),
        [&](NodeId nodeId) {
            const NodeId node = detail::grayStatsSlotOf(context.tree, nodeId);
            const T nodeAltitude = context.altitude[static_cast<std::size_t>(nodeId)];
            const Real nodeAltitudeAsReal = static_cast<Real>(nodeAltitude);
            if (request.grayLevelVariance) {
                const double nodeAltitudeAsDouble = static_cast<double>(nodeAltitude);
                sumGrayLevelSquare[static_cast<std::size_t>(node)] =
                    static_cast<double>(::mmcfilters::detail::CommittedTreeAccess::properPartCardinality(context.tree, nodeId)) * nodeAltitudeAsDouble *
                    nodeAltitudeAsDouble;
            }
            if (request.grayLevelHeight) {
                subtreeMinAltitude[node] = nodeAltitudeAsReal;
                subtreeMaxAltitude[node] = nodeAltitudeAsReal;
            }
        },
        [&](NodeId parentNodeId, NodeId childNodeId) {
            const NodeId parent = detail::grayStatsSlotOf(context.tree, parentNodeId);
            const NodeId child = detail::grayStatsSlotOf(context.tree, childNodeId);
            if (request.grayLevelVariance)
                sumGrayLevelSquare[parent] += sumGrayLevelSquare[child];
            if (request.grayLevelHeight) {
                subtreeMinAltitude[parent] = std::min(subtreeMinAltitude[parent], subtreeMinAltitude[child]);
                subtreeMaxAltitude[parent] = std::max(subtreeMaxAltitude[parent], subtreeMaxAltitude[child]);
            }
        },
        [&](NodeId nodeId) {
            const NodeId node = detail::grayStatsSlotOf(context.tree, nodeId);
            Real area = needsAggregateDependencies ? areaDependency->buffer[indexOfArea(node)] : Real{0};
            if (request.meanGrayLevel)
                context.buffer[indexOfMean(node)] =
                    ::mmcfilters::attributes::numeric::safeDivide(volumeDependency->buffer[indexOfVolume(node)], area);
            if (request.grayLevelVariance) {
                Real meanGrayLevel = ::mmcfilters::attributes::numeric::safeDivide(volumeDependency->buffer[indexOfVolume(node)], area);
                double meanGrayLevelSquare =
                    ::mmcfilters::attributes::numeric::safeDivide(sumGrayLevelSquare[static_cast<std::size_t>(node)], static_cast<double>(area));
                Real var = static_cast<Real>(meanGrayLevelSquare - (static_cast<double>(meanGrayLevel) * static_cast<double>(meanGrayLevel)));
                context.buffer[indexOfVariance(node)] = ::mmcfilters::attributes::numeric::clampNonNegative(var);
            }
            if (request.grayLevelHeight) {
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
    TreeAltitudeAlgorithms::validateNodeAltitudeBufferShape(context.tree, context.altitude);
    const auto requireColumn = [&](bool requested, Attribute attribute, const char* name) {
        if (requested && !context.attrNames.contains(attribute)) {
            throw std::invalid_argument(std::string(name) + " computation requires a matching output column.");
        }
    };
    requireColumn(request.meanGrayLevel, MeanGrayLevel, "MEAN_GRAY_LEVEL");
    requireColumn(request.grayLevelVariance, GrayLevelVariance, "GRAY_LEVEL_VARIANCE");
    requireColumn(request.grayLevelHeight, GrayLevelHeight, "GRAY_LEVEL_HEIGHT");
}

} // namespace detail

/**
 * @brief Computes grey-level statistics derived from subtree aggregation.
 *
 * @details
 * This computer provides three complementary descriptors:
 * - `MeanGrayLevel`: the average grey level over the full subtree support;
 * - `GrayLevelVariance`: the grey-level variance over the full subtree support;
 * - `GrayLevelHeight`: the maximum absolute altitude difference between the node
 *   and any node in its subtree.
 *
 * `MeanGrayLevel` is obtained from `VOLUME / AREA`, while
 * `GrayLevelVariance` additionally requires the accumulated sum of squared grey
 * levels over the subtree support. `GrayLevelHeight` propagates both the minimum
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
    inline static constexpr std::array<Attribute, 3> producedAttributes{MeanGrayLevel, GrayLevelVariance, GrayLevelHeight};

    /**
     * @brief Computes the requested grey-level statistics.
     *
     * @details
     * Requires a typed altitude span in dense internal-node order. `GrayLevelHeight`
     * uses altitude directly. `MeanGrayLevel` and `GrayLevelVariance`
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
                volumeDependency = &context.dependencies.require(Volume);
                areaDependency = &context.dependencies.require(Area);
            } else {
                volumeDependency = detail::findGrayDependency(context.dependencySources, Volume);
                areaDependency = detail::findGrayDependency(context.dependencySources, Area);
            }
        }
        detail::kernel::computeGrayLevelStats(context, request, volumeDependency, areaDependency);
    }

    /**
     * @brief Materializes grey-level statistics for one-pixel unit supports.
     *
     * Unit `MeanGrayLevel` equals the smallest-node altitude. Unit
     * `GrayLevelVariance` and `GrayLevelHeight` are zero.
     *
     * @param context Operation context or diagnostic label.
     */
    template <std::floating_point Real, AltitudeValue T> static void computeUnitRows(const AltitudeUnitAttributeComputeContext<Real, T>& context) {
        requireUnitAttributeBufferShape(context.tree, context.unitPixels, context.buffer, context.attrNames);
        TreeAltitudeAlgorithms::validateNodeAltitudeBufferShape(context.tree, context.altitude);

        const detail::GrayLevelStatsRequest request = detail::GrayLevelStatsRequest::from(context.requestedAttributes);
        if (!request.any()) {
            return;
        }

        for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(context.unitPixels.size()); ++leafIndex) {
            const PixelId pixel = context.unitPixels[static_cast<size_t>(leafIndex)];
            const Real altitudeValue = request.meanGrayLevel
                                           ? static_cast<Real>(::mmcfilters::detail::kernel::unitAltitude(context.tree, context.altitude, pixel))
                                           : Real{0};

            if (request.meanGrayLevel) {
                context.buffer[context.attrNames.linearIndex(leafIndex, MeanGrayLevel)] = altitudeValue;
            }
            if (request.grayLevelVariance) {
                context.buffer[context.attrNames.linearIndex(leafIndex, GrayLevelVariance)] = Real{0};
            }
            if (request.grayLevelHeight) {
                context.buffer[context.attrNames.linearIndex(leafIndex, GrayLevelHeight)] = Real{0};
            }
        }
    }
};

} // namespace mmcfilters::attributes::computers
