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

namespace mmcfilters::attributes::computers {

namespace detail {
/**
 * @brief Returns the dense attribute-buffer slot for a tree node.
 *
 * @param nodeId Identifier of the node used by the operation.
 * @return Dense attribute-buffer slot for the node.
 */
inline NodeId volumeSlotOf(const MorphologicalTree&, NodeId nodeId) noexcept { return nodeId; }

/** @brief Describes the volume attributes requested from the attribute pipeline. */
struct VolumeRequest {
    /** @brief Indicates whether the volume attribute was requested. */
    bool volume = false;
    /** @brief Indicates whether relative volume was requested. */
    bool relative = false;

    /**
     * @brief Tests whether any requested feature is enabled.
     *
     * @return True when any requested feature is enabled; otherwise false.
     */
    [[nodiscard]] bool any() const noexcept { return volume || relative; }

    /**
     * @brief Tests whether area dependency holds.
     *
     * @return True when area dependency; otherwise false.
     */
    [[nodiscard]] bool needsAreaDependency() const noexcept { return relative; }

    /**
     * @brief Builds a request descriptor from the requested attributes.
     *
     * @param requestedAttributes Requested attribute subset.
     * @return Resulting request descriptor from the requested attributes.
     */
    [[nodiscard]] static VolumeRequest from(std::span<const Attribute> requestedAttributes) {
        return {.volume = containsVolumeAttribute(requestedAttributes, VOLUME), .relative = containsVolumeAttribute(requestedAttributes, RELATIVE_VOLUME)};
    }

  private:
    /**
     * @brief Tests whether volume attribute holds.
     *
     * @param requestedAttributes Requested attribute subset.
     * @param attribute Attribute requested by the operation.
     * @return True when volume attribute; otherwise false.
     */
    [[nodiscard]] static bool containsVolumeAttribute(std::span<const Attribute> requestedAttributes, Attribute attribute) {
        return std::find(requestedAttributes.begin(), requestedAttributes.end(), attribute) != requestedAttributes.end();
    }
};

/**
 * @brief Typed kernel for `VOLUME` and `RELATIVE_VOLUME`.
 *
 * @param tree Tree whose dense internal node ids index `buffer`.
 * @param altitude Dense altitude span indexed by internal node id.
 * @param buffer Flat internal-node output buffer.
 * @param attrNames Layout containing requested volume attributes.
 * @param requestedAttributes Requested subset of `VOLUME` and
 * `RELATIVE_VOLUME`.
 * @param dependencySources Dependency source `0` must contain `AREA` when
 * `RELATIVE_VOLUME` is requested.
 */
namespace kernel {

/**
 * @brief Computes requested volume descriptors from established altitude and dependency domains.
 * @param context Established tree, altitude span, output layout, and output buffer.
 * @param request Volume columns to materialize.
 * @param areaDependency Optional established area dependency.
 */
template <std::floating_point Real, AltitudeValue T>
void computeVolume(const AltitudeAttributeComputeContext<Real, T>& context, const VolumeRequest& request,
                   const DependencySourceT<Real>* areaDependency) {
    if (!request.any()) {
        return;
    }

    const int stride = context.attrNames.NUM_ATTRIBUTES;
    const auto offsetOf = [&](Attribute attribute) { return context.attrNames.indexMap.find(attribute)->second; };
    const int volumeOffset = request.volume ? offsetOf(VOLUME) : 0;
    const int relativeOffset = request.relative ? offsetOf(RELATIVE_VOLUME) : 0;
    auto indexOfVolume = [&](NodeId node) { return static_cast<std::size_t>(node * stride + volumeOffset); };
    auto indexOfRelative = [&](NodeId node) { return static_cast<std::size_t>(node * stride + relativeOffset); };
    const int areaStride = areaDependency != nullptr ? areaDependency->attrNames->NUM_ATTRIBUTES : 0;
    const int areaOffset = areaDependency != nullptr ? areaDependency->attrNames->indexMap.find(AREA)->second : 0;
    auto indexOfArea = [&](NodeId node) { return static_cast<std::size_t>(node * areaStride + areaOffset); };

    ::mmcfilters::detail::kernel::traversePostOrder(
        context.tree, context.tree.getRoot(),
        [&](NodeId nodeId) {
            const NodeId node = detail::volumeSlotOf(context.tree, nodeId);
            const T nodeAltitude = context.altitude[static_cast<std::size_t>(nodeId)];
            if (request.volume)
                context.buffer[indexOfVolume(node)] =
                    static_cast<Real>(::mmcfilters::detail::CommittedTreeAccess::numProperParts(context.tree, nodeId)) * static_cast<Real>(nodeAltitude);
            if (request.relative)
                context.buffer[indexOfRelative(node)] = Real{0};
        },
        [&](NodeId parentNodeId, NodeId childNodeId) {
            const NodeId parent = detail::volumeSlotOf(context.tree, parentNodeId);
            const NodeId child = detail::volumeSlotOf(context.tree, childNodeId);
            if (request.volume)
                context.buffer[indexOfVolume(parent)] += context.buffer[indexOfVolume(child)];
            if (request.relative)
                context.buffer[indexOfRelative(parent)] +=
                    context.buffer[indexOfRelative(child)] +
                    static_cast<Real>(areaDependency->buffer[indexOfArea(child)] *
                                      std::abs(static_cast<Real>(context.altitude[static_cast<std::size_t>(childNodeId)]) -
                                               static_cast<Real>(context.altitude[static_cast<std::size_t>(parentNodeId)])));
        },
        [&](NodeId nodeId) {
            const NodeId node = detail::volumeSlotOf(context.tree, nodeId);
            if (request.relative)
                context.buffer[indexOfRelative(node)] += areaDependency->buffer[indexOfArea(node)];
        });
}

} // namespace kernel

template <std::floating_point Real>
inline const DependencySourceT<Real>* findDependency(std::span<const DependencySourceT<Real>> sources, Attribute attribute) {
    for (const DependencySourceT<Real>& source : sources) {
        if (source.attrNames->contains(attribute)) {
            return &source;
        }
    }
    return nullptr;
}

template <std::floating_point Real, AltitudeValue T>
inline void validateVolumeContext(const AltitudeAttributeComputeContext<Real, T>& context, const VolumeRequest& request) {
    requireAttributeBufferShape(context.tree, context.buffer, context.attrNames);
    TreeAltitudeAlgorithms::validateAltitudeBufferShape(context.tree, context.altitude);
    if (request.volume && !context.attrNames.contains(VOLUME)) {
        throw std::invalid_argument("VOLUME computation requires a VOLUME column in the output layout.");
    }
    if (request.relative && !context.attrNames.contains(RELATIVE_VOLUME)) {
        throw std::invalid_argument("RELATIVE_VOLUME computation requires a RELATIVE_VOLUME column in the output layout.");
    }
}

} // namespace detail

/**
 * @brief Computes cumulative grey-level volume descriptors on the hierarchy.
 *
 * @details
 * This computer exposes two related attributes:
 * - `VOLUME`: the sum of `altitude * area contribution` over the full subtree
 *   support of the node;
 * - `RELATIVE_VOLUME`: the cumulative contribution of grey-level jumps between
 *   a node and its descendants, weighted by subtree area.
 *
 * `VOLUME` is initialised from the node's own proper parts and then aggregated
 * bottom-up. `RELATIVE_VOLUME` is initialised at zero and accumulates, for
 * each parent/child edge, the absolute altitude difference multiplied by the
 * child subtree area. A final per-node area term is then added so that the
 * descriptor remains expressed in the same unit as a weighted support size.
 *
 * @note `RELATIVE_VOLUME` depends on `AREA` because subtree sizes are needed
 * to weight parent/child altitude differences.
 */
class VolumeComputer {
  public:
    /// Family name used in dependency-plan diagnostics.
    static constexpr std::string_view familyName = "volume";

    /// Stable family id used by the scheduler.
    static constexpr AttributeComputerFamily family = AttributeComputerFamily::Volume;

    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Altitude;

    /**
     * @brief Canonical list of volume descriptors produced by this computer.
     */
    inline static constexpr std::array<Attribute, 2> producedAttributes{VOLUME, RELATIVE_VOLUME};

    /**
     * @brief Computes the requested volume descriptors.
     *
     * @details
     * Requires a typed altitude span in dense internal-node order.
     * `context.requestedAttributes` selects `VOLUME`, `RELATIVE_VOLUME`, or
     * both. `RELATIVE_VOLUME` additionally requires an `AREA` dependency
     * available through `context.dependencies`.
     *
     * @param context Operation context or diagnostic label.
     */
    template <std::floating_point Real, AltitudeValue T> static void compute(const AltitudeAttributeComputeContext<Real, T>& context) {
        const detail::VolumeRequest request = detail::VolumeRequest::from(context.requestedAttributes);
        MMCFILTERS_CONTRACT_CHECKED_ONLY(detail::validateVolumeContext(context, request));
        const DependencySourceT<Real>* areaDependency = nullptr;
        if (request.needsAreaDependency()) {
            if constexpr (contract::validationsEnabled) {
                areaDependency = &context.dependencies.require(AREA);
            } else {
                areaDependency = detail::findDependency(context.dependencySources, AREA);
            }
        }
        detail::kernel::computeVolume(context, request, areaDependency);
    }

    /**
     * @brief Materializes volume descriptors for one-pixel unit supports.
     *
     * Unit `VOLUME` is the owner-node altitude of the proper part. Unit
     * `RELATIVE_VOLUME` is `1`, matching the scalar support-size unit used by
     * the subtree computation.
     *
     * @param context Operation context or diagnostic label.
     */
    template <std::floating_point Real, AltitudeValue T> static void computeUnitRows(const AltitudeUnitAttributeComputeContext<Real, T>& context) {
        requireUnitAttributeBufferShape(context.tree, context.unitProperParts, context.buffer, context.attrNames);
        TreeAltitudeAlgorithms::validateAltitudeBufferShape(context.tree, context.altitude);

        const detail::VolumeRequest request = detail::VolumeRequest::from(context.requestedAttributes);
        if (!request.any()) {
            return;
        }

        for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(context.unitProperParts.size()); ++leafIndex) {
            const NodeId properPart = context.unitProperParts[static_cast<size_t>(leafIndex)];
            if (request.volume) {
                context.buffer[context.attrNames.linearIndex(leafIndex, VOLUME)] =
                    static_cast<Real>(::mmcfilters::detail::kernel::unitAltitude(context.tree, context.altitude, properPart));
            }
            if (request.relative) {
                context.buffer[context.attrNames.linearIndex(leafIndex, RELATIVE_VOLUME)] = Real{1};
            }
        }
    }
};

} // namespace mmcfilters::attributes::computers
