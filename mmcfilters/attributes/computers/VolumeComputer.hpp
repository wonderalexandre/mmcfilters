#pragma once

#include "AttributeComputerDomain.hpp"
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

namespace mmcfilters::attributes::computers {

namespace detail {
inline NodeId volumeSlotOf(const MorphologicalTree&, NodeId nodeId) noexcept {
    return nodeId;
}

struct VolumeRequest {
    bool volume = false;
    bool relative = false;

    [[nodiscard]] bool any() const noexcept {
        return volume || relative;
    }

    [[nodiscard]] bool needsAreaDependency() const noexcept {
        return relative;
    }

    [[nodiscard]] static VolumeRequest from(std::span<const Attribute> requestedAttributes) {
        return {
            .volume = containsVolumeAttribute(requestedAttributes, VOLUME),
            .relative = containsVolumeAttribute(requestedAttributes, RELATIVE_VOLUME)};
    }

private:
    [[nodiscard]] static bool containsVolumeAttribute(
        std::span<const Attribute> requestedAttributes,
        Attribute attribute)
    {
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
template<std::floating_point Real, AltitudeValue T>
void computeVolumeAttributeKernel(
    const MorphologicalTree& tree,
    std::span<const T> altitude,
    std::span<Real> buffer,
    const AttributeNames& attrNames,
    std::span<const Attribute> requestedAttributes,
    std::span<const DependencySourceT<Real>> dependencySources)
{
    TreeAltitudeAlgorithms::validateAltitudeBufferShape(tree, altitude);

    const VolumeRequest request = VolumeRequest::from(requestedAttributes);
    if (!request.any()) {
        return;
    }

    auto indexOfVol = [&](NodeId idx) { return attrNames.linearIndex(idx, VOLUME); };
    auto indexOfRel = [&](NodeId idx) { return attrNames.linearIndex(idx, RELATIVE_VOLUME); };
    const DependencyResolver<Real> dependencies{dependencySources};
    const DependencySourceT<Real>* dependencyArea = request.needsAreaDependency()
        ? &dependencies.require(AREA)
        : nullptr;
    auto indexOfArea = [&](NodeId idx) { return dependencyArea->attrNames->linearIndex(idx, AREA); };

    ::mmcfilters::detail::traversePostOrder(
        tree,
        tree.getRoot(),
        [&](NodeId nodeId) {
            const NodeId node = detail::volumeSlotOf(tree, nodeId);
            const T nodeAltitude = TreeAltitudeAlgorithms::getAltitude(altitude, nodeId);
            if (request.volume)
                buffer[indexOfVol(node)] =
                    static_cast<Real>(tree.getNumProperParts(nodeId)) * static_cast<Real>(nodeAltitude);
            if (request.relative)
                buffer[indexOfRel(node)] = Real{0};
        },
        [&](NodeId parentNodeId, NodeId childNodeId) {
            const NodeId parent = detail::volumeSlotOf(tree, parentNodeId);
            const NodeId child = detail::volumeSlotOf(tree, childNodeId);
            if (request.volume)
                buffer[indexOfVol(parent)] += buffer[indexOfVol(child)];
            if (request.relative)
                buffer[indexOfRel(parent)] +=
                    buffer[indexOfRel(child)] +
                    static_cast<Real>(
                        dependencyArea->buffer[indexOfArea(child)] *
                        std::abs(
                            static_cast<Real>(TreeAltitudeAlgorithms::getAltitude(altitude, childNodeId)) -
                            static_cast<Real>(TreeAltitudeAlgorithms::getAltitude(altitude, parentNodeId))));
        },
        [&](NodeId nodeId) {
            const NodeId node = detail::volumeSlotOf(tree, nodeId);
            if (request.relative)
                buffer[indexOfRel(node)] += dependencyArea->buffer[indexOfArea(node)];
        });
}

}

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

    /// Execution domain required by the computer.
    static constexpr AttributeComputerDomain domain = AttributeComputerDomain::Altitude;

    /**
     * @brief Canonical list of volume descriptors produced by this computer.
     */
    inline static constexpr std::array<Attribute, 2> producedAttributes{
        VOLUME,
        RELATIVE_VOLUME};

    /**
     * @brief Computes the requested volume descriptors.
     *
     * @details
     * Requires a typed altitude span in dense internal-node order.
     * `context.requestedAttributes` selects `VOLUME`, `RELATIVE_VOLUME`, or
     * both. `RELATIVE_VOLUME` additionally requires an `AREA` dependency
     * available through `context.dependencies`.
     */
    template <std::floating_point Real, AltitudeValue T>
    static void compute(const AltitudeAttributeComputeContext<Real, T>& context)
    {
        requireAttributeBufferShape(context.tree, context.buffer, context.attrNames);
        detail::computeVolumeAttributeKernel(
            context.tree,
            context.altitude,
            context.buffer,
            context.attrNames,
            context.requestedAttributes,
            context.dependencySources);
    }

    /**
     * @brief Materializes volume descriptors for one-pixel unit supports.
     *
     * Unit `VOLUME` is the owner-node altitude of the proper part. Unit
     * `RELATIVE_VOLUME` is `1`, matching the scalar support-size unit used by
     * the subtree computation.
     */
    template <std::floating_point Real, AltitudeValue T>
    static void computeUnitRows(const AltitudeUnitAttributeComputeContext<Real, T>& context)
    {
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
                    static_cast<Real>(unitAltitude(context.tree, context.altitude, properPart));
            }
            if (request.relative) {
                context.buffer[context.attrNames.linearIndex(leafIndex, RELATIVE_VOLUME)] = Real{1};
            }
        }
    }

};

} // namespace mmcfilters::attributes::computers
