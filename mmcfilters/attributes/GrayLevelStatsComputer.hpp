#pragma once

#include "AttributeComputer.hpp"
#include "AttributeComputedIncrementally.hpp"
#include "../trees/TreeAltitudeOps.hpp"
namespace mmcfilters {

namespace detail {
inline NodeId grayStatsSlotOf(const MorphologicalTree&, NodeId nodeId) noexcept {
    return nodeId;
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
class GrayLevelStatsComputer : public AttributeComputer {
public:
    /**
     * @brief Returns the grey-level descriptors naturally produced together.
     */
    std::vector<Attribute> attributes() const override {
        return {LEVEL, MEAN_LEVEL, VARIANCE_LEVEL, GRAY_HEIGHT};
    }

    /**
     * @brief Declares the dependencies required by the derived statistics.
     */
    std::vector<AttributeOrGroup> requiredAttributes() const override {
        return {VOLUME, AREA};
    }

    /**
     * @brief Computes the requested grey-level statistics.
     */
    void compute(
        const MorphologicalTree& tree,
        const AltitudeBuffer* altitude,
        std::span<float> buffer,
        const AttributeNames& attrNames,
        std::span<const Attribute> requestedAttributes,
        std::span<const DependencySource> dependencySources) const override
    {
        if (PRINT_LOG) std::cout << "\n==== AttributeComputer: Computing GrayLevelStatsComputer " << std::endl;

        auto indexOfMean = [&](NodeId idx) { return attrNames.linearIndex(idx, MEAN_LEVEL); };
        auto indexOfLevel = [&](NodeId idx) { return attrNames.linearIndex(idx, LEVEL); };
        auto indexOfVariance = [&](NodeId idx) { return attrNames.linearIndex(idx, VARIANCE_LEVEL); };
        auto indexOfGrayHeight = [&](NodeId idx) { return attrNames.linearIndex(idx, GRAY_HEIGHT); };

        bool computeMeanLevel = std::find(requestedAttributes.begin(), requestedAttributes.end(), MEAN_LEVEL) != requestedAttributes.end();
        bool computeVarianceLevel = std::find(requestedAttributes.begin(), requestedAttributes.end(), VARIANCE_LEVEL) != requestedAttributes.end();
        bool computeLevel = std::find(requestedAttributes.begin(), requestedAttributes.end(), LEVEL) != requestedAttributes.end();
        bool computeGrayHeight = std::find(requestedAttributes.begin(), requestedAttributes.end(), GRAY_HEIGHT) != requestedAttributes.end();

        const auto& dependencyVol = dependencySources[0];
        const auto& dependencyArea = dependencySources[1];
        auto indexOfVol = [&](NodeId idx) { return dependencyVol.attrNames->linearIndex(idx, VOLUME); };
        auto indexOfArea = [&](NodeId idx) { return dependencyArea.attrNames->linearIndex(idx, AREA); };

        std::vector<long> sumGrayLevelSquare;
        if (computeVarianceLevel) {
            sumGrayLevelSquare.assign(tree.getNumInternalNodeSlots(), 0L);
        }

        AttributeComputedIncrementally::traversePostOrder(
            tree,
            tree.getRoot(),
            [&](NodeId nodeId) {
                const NodeId node = detail::grayStatsSlotOf(tree, nodeId);
                const AltitudeType nodeAltitude = tree_altitude_ops::getAltitude(altitude, nodeId);
                if (computeVarianceLevel)
                    sumGrayLevelSquare[node] = static_cast<long>(tree.getNumProperParts(nodeId) * std::pow(nodeAltitude, 2));
                if (computeLevel)
                    buffer[indexOfLevel(node)] = static_cast<float>(nodeAltitude);
                if (computeGrayHeight)
                    buffer[indexOfGrayHeight(node)] = static_cast<float>(nodeAltitude);
            },
            [&](NodeId parentNodeId, NodeId childNodeId) {
                const NodeId parent = detail::grayStatsSlotOf(tree, parentNodeId);
                const NodeId child = detail::grayStatsSlotOf(tree, childNodeId);
                if (computeVarianceLevel)
                    sumGrayLevelSquare[parent] += sumGrayLevelSquare[child];
                if (computeGrayHeight) {
                    float childValue = buffer[indexOfGrayHeight(child)];
                    float& parentValue = buffer[indexOfGrayHeight(parent)];
                    if (tree.isMaxtree())
                        parentValue = std::max(parentValue, childValue);
                    else
                        parentValue = std::min(parentValue, childValue);
                }
            },
            [&](NodeId nodeId) {
                const NodeId node = detail::grayStatsSlotOf(tree, nodeId);
                float area = dependencyArea.buffer[indexOfArea(node)];
                if (computeMeanLevel)
                    buffer[indexOfMean(node)] = dependencyVol.buffer[indexOfVol(node)] / area;
                if (computeVarianceLevel) {
                    float meanGrayLevel = dependencyVol.buffer[indexOfVol(node)] / area;
                    double meanGrayLevelSquare = sumGrayLevelSquare[node] / area;
                    float var = static_cast<float>(meanGrayLevelSquare - (meanGrayLevel * meanGrayLevel));
                    buffer[indexOfVariance(node)] = var > 0.0f ? var : 0.0f;
                }
            });

        if (computeGrayHeight) {
            for (NodeId nodeId : tree.getPostOrderNodes()) {
                const NodeId node = detail::grayStatsSlotOf(tree, nodeId);
                if (tree.isLeaf(nodeId))
                    buffer[indexOfGrayHeight(node)] = 0.0f;
                else
                    buffer[indexOfGrayHeight(node)] = std::abs(static_cast<float>(tree_altitude_ops::getAltitude(altitude, nodeId)) - buffer[indexOfGrayHeight(node)]) + 1.0f;
            }
        }
    }
};

} // namespace mmcfilters
