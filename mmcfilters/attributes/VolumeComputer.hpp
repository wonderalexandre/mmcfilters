#pragma once

#include "AttributeComputer.hpp"
#include "AttributeComputedIncrementally.hpp"
#include "../trees/WeightedMorphologicalTree.hpp"
namespace mmcfilters {

namespace detail {
inline NodeId volumeSlotOf(const MorphologicalTree&, NodeId nodeId) noexcept {
    return nodeId;
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
    /**
     * @brief Returns the volume descriptors produced by this computer.
     */
    std::vector<Attribute> attributes() const override { return {VOLUME, RELATIVE_VOLUME}; }

    /**
     * @brief Declares the dependencies required by `RELATIVE_VOLUME`.
     */
    std::vector<AttributeOrGroup> requiredAttributes() const override { return {AREA}; }

    /**
     * @brief Computes the requested volume descriptors.
     */
    void compute(
        const MorphologicalTree& tree,
        const AltitudeBuffer* altitude,
        std::span<float> buffer,
        const AttributeNames& attrNames,
        std::span<const Attribute> requestedAttributes,
        std::span<const DependencySource> dependencySources) const override
    {
        if (PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing VOLUME" << std::endl;
        auto indexOfVol = [&](NodeId idx) { return attrNames.linearIndex(idx, VOLUME); };
        auto indexOfRel = [&](NodeId idx) { return attrNames.linearIndex(idx, RELATIVE_VOLUME); };
        const auto& dependencyArea = dependencySources[0];
        auto indexOfArea = [&](NodeId idx) { return dependencyArea.attrNames->linearIndex(idx, AREA); };

        bool computeVolume = std::find(requestedAttributes.begin(), requestedAttributes.end(), VOLUME) != requestedAttributes.end();
        bool computeRelative = std::find(requestedAttributes.begin(), requestedAttributes.end(), RELATIVE_VOLUME) != requestedAttributes.end();

        AttributeComputedIncrementally::traversePostOrder(
            tree,
            tree.getRoot(),
            [&](NodeId nodeId) {
                const NodeId node = detail::volumeSlotOf(tree, nodeId);
                const AltitudeType nodeAltitude = WeightedMorphologicalTree::getAltitude(altitude, nodeId);
                if (computeVolume)
                    buffer[indexOfVol(node)] = static_cast<float>(tree.getNumProperParts(nodeId) * nodeAltitude);
                if (computeRelative)
                    buffer[indexOfRel(node)] = 0.0f;
            },
            [&](NodeId parentNodeId, NodeId childNodeId) {
                const NodeId parent = detail::volumeSlotOf(tree, parentNodeId);
                const NodeId child = detail::volumeSlotOf(tree, childNodeId);
                if (computeVolume)
                    buffer[indexOfVol(parent)] += buffer[indexOfVol(child)];
                if (computeRelative)
                    buffer[indexOfRel(parent)] +=
                        buffer[indexOfRel(child)] +
                        static_cast<float>(
                            dependencyArea.buffer[indexOfArea(child)] *
                            std::abs(
                                static_cast<float>(WeightedMorphologicalTree::getAltitude(altitude, childNodeId)) -
                                static_cast<float>(WeightedMorphologicalTree::getAltitude(altitude, parentNodeId))));
            },
            [&](NodeId nodeId) {
                const NodeId node = detail::volumeSlotOf(tree, nodeId);
                if (computeRelative)
                    buffer[indexOfRel(node)] += dependencyArea.buffer[indexOfArea(node)];
            });
    }
};

} // namespace mmcfilters
