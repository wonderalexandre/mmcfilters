#pragma once

#include "../AttributeComputer.hpp"
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
inline NodeId volumeSlotOf(const MorphologicalTree&, NodeId nodeId) noexcept {
    return nodeId;
}

inline constexpr std::array<Attribute, 2> VOLUME_ATTRIBUTES{
    VOLUME,
    RELATIVE_VOLUME};

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
template<AltitudeValue T>
void computeVolumeAttributes(
    const MorphologicalTree& tree,
    std::span<const T> altitude,
    std::span<float> buffer,
    const AttributeNames& attrNames,
    std::span<const Attribute> requestedAttributes,
    std::span<const DependencySource> dependencySources)
{
    TreeAltitudeAlgorithms::validateAltitudeBufferShape(tree, altitude);

    const VolumeRequest request = VolumeRequest::from(requestedAttributes);
    if (!request.any()) {
        return;
    }

    auto indexOfVol = [&](NodeId idx) { return attrNames.linearIndex(idx, VOLUME); };
    auto indexOfRel = [&](NodeId idx) { return attrNames.linearIndex(idx, RELATIVE_VOLUME); };
    const DependencySource* dependencyArea = request.needsAreaDependency()
        ? &requireDependencySourceForAttribute(dependencySources, 0, AREA)
        : nullptr;
    auto indexOfArea = [&](NodeId idx) { return dependencyArea->attrNames->linearIndex(idx, AREA); };

    ::mmcfilters::detail::traversePostOrder(
        tree,
        tree.getRoot(),
        [&](NodeId nodeId) {
            const NodeId node = detail::volumeSlotOf(tree, nodeId);
            const T nodeAltitude = TreeAltitudeAlgorithms::getAltitude(altitude, nodeId);
            if (request.volume)
                buffer[indexOfVol(node)] = static_cast<float>(tree.getNumProperParts(nodeId) * nodeAltitude);
            if (request.relative)
                buffer[indexOfRel(node)] = 0.0f;
        },
        [&](NodeId parentNodeId, NodeId childNodeId) {
            const NodeId parent = detail::volumeSlotOf(tree, parentNodeId);
            const NodeId child = detail::volumeSlotOf(tree, childNodeId);
            if (request.volume)
                buffer[indexOfVol(parent)] += buffer[indexOfVol(child)];
            if (request.relative)
                buffer[indexOfRel(parent)] +=
                    buffer[indexOfRel(child)] +
                    static_cast<float>(
                        dependencyArea->buffer[indexOfArea(child)] *
                        std::abs(
                            static_cast<float>(TreeAltitudeAlgorithms::getAltitude(altitude, childNodeId)) -
                            static_cast<float>(TreeAltitudeAlgorithms::getAltitude(altitude, parentNodeId))));
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
class VolumeComputer : public AttributeComputer {
public:
    using AttributeComputer::compute;
    using AttributeComputer::computeUnitAttributes;

    /**
     * @brief Returns the volume descriptors produced by this computer.
     */
    [[nodiscard]] std::vector<Attribute> attributes() const override {
        return {detail::VOLUME_ATTRIBUTES.begin(), detail::VOLUME_ATTRIBUTES.end()};
    }

    /**
     * @brief Computes the requested volume descriptors.
     *
     * Requires an altitude span. `RELATIVE_VOLUME` also requires one dependency
     * source containing `AREA`.
     */
    void compute(
        const MorphologicalTree& tree,
        AttributeAltitudeView altitude,
        std::span<float> buffer,
        const AttributeNames& attrNames,
        std::span<const Attribute> requestedAttributes,
        std::span<const DependencySource> dependencySources) const override
    {
        requireAttributeBufferShape(tree, buffer, attrNames);
        detail::computeVolumeAttributes(
            tree,
            TreeAltitudeAlgorithms::requireAltitudeSpan(altitude),
            buffer,
            attrNames,
            requestedAttributes,
            dependencySources);
    }

    /**
     * @brief Materializes volume descriptors for one-pixel unit supports.
     *
     * Unit `VOLUME` is the owner-node altitude of the proper part. Unit
     * `RELATIVE_VOLUME` is `1`, matching the scalar support-size unit used by
     * the subtree computation.
     */
    void computeUnitAttributes(
        const MorphologicalTree& tree,
        AttributeAltitudeView altitude,
        std::span<const NodeId> unitProperParts,
        std::span<float> buffer,
        const AttributeNames& attrNames,
        std::span<const Attribute> requestedAttributes) const override
    {
        requireUnitAttributeBufferShape(tree, unitProperParts, buffer, attrNames);

        const detail::VolumeRequest request = detail::VolumeRequest::from(requestedAttributes);
        if (!request.any()) {
            return;
        }

        for (NodeId leafIndex = 0; leafIndex < static_cast<NodeId>(unitProperParts.size()); ++leafIndex) {
            const NodeId properPart = unitProperParts[static_cast<size_t>(leafIndex)];
            if (request.volume) {
                buffer[attrNames.linearIndex(leafIndex, VOLUME)] =
                    static_cast<float>(unitAltitude(tree, altitude, properPart));
            }
            if (request.relative) {
                buffer[attrNames.linearIndex(leafIndex, RELATIVE_VOLUME)] = 1.0f;
            }
        }
    }
};

} // namespace mmcfilters::attributes::computers
