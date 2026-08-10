#pragma once

#include "../MorphologicalTree.hpp"
#include "CommittedTreeAccess.hpp"
#include "HierarchyCapabilityValidation.hpp"
#include "../../utils/Altitude.hpp"
#include "../../utils/Common.hpp"
#include "../../utils/Contract.hpp"

#include <cmath>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace mmcfilters::detail {

/**
 * @brief Utilities for altitude-radius neighbourhoods on a morphological tree.
 *
 * @details
 * These helpers operate on topology plus an explicit dense altitude buffer.
 * They are shared by delta-attribute materialization and MSER stability, so the
 * implementation lives beside tree/altitude algorithms instead of under the
 * attribute-computer directory.
 */

/**
 * @brief Computes subtree areas in the dense internal-node slot space.
 *
 * @details
 * The area is the number of direct proper parts owned by a node plus the areas
 * accumulated from its descendants. Delta attributes use this only to select a
 * representative descendant when several nodes map to the same ascendant.
 *
 * @param tree Tree topology used by the operation.
 * @return Values produced by the operation.
 */
inline std::vector<int32_t> computeAttributeDeltaAreasIncrementally(const MorphologicalTree& tree) {
    std::vector<int32_t> area(static_cast<size_t>(tree.getNumInternalNodeSlots()), 0);
    for (NodeId nodeId : tree.getPostOrderNodes()) {
        area[static_cast<size_t>(nodeId)] += static_cast<int32_t>(CommittedTreeAccess::numProperParts(tree, nodeId));
        const NodeId parentNodeId = CommittedTreeAccess::nodeParent(tree, nodeId);
        if (parentNodeId != InvalidNode && parentNodeId != nodeId) {
            area[static_cast<size_t>(parentNodeId)] += area[static_cast<size_t>(nodeId)];
        }
    }
    return area;
}

/**
 * @brief Keeps the largest-area descendant assigned to one delta ascendant.
 *
 * @details
 * Ties are resolved by the smallest `NodeId` to make the mapping deterministic.
 * The `descendants` buffer is indexed by dense internal `NodeId`.
 *
 * @param tree Tree topology used by the operation.
 * @param areaByNode Node identifier represented by `areaByNode`.
 * @param descendants Descendant-node data represented by `descendants`.
 * @param ascendantNodeId Node identifier represented by `ascendantNodeId`.
 * @param candidateNodeId Node identifier represented by `candidateNodeId`.
 */
inline void updateLargestAreaAttributeDeltaDescendant(const MorphologicalTree& tree, const std::vector<int32_t>& areaByNode, std::vector<NodeId>& descendants,
                                                      NodeId ascendantNodeId, NodeId candidateNodeId) {
    if (!tree.isNode(ascendantNodeId) || !tree.isNode(candidateNodeId)) {
        return;
    }
    NodeId& currentNodeId = descendants[static_cast<size_t>(ascendantNodeId)];
    if (ascendantNodeId == candidateNodeId) {
        return;
    }
    if (currentNodeId == InvalidNode || areaByNode[static_cast<size_t>(candidateNodeId)] > areaByNode[static_cast<size_t>(currentNodeId)] ||
        (areaByNode[static_cast<size_t>(candidateNodeId)] == areaByNode[static_cast<size_t>(currentNodeId)] && candidateNodeId < currentNodeId)) {
        currentNodeId = candidateNodeId;
    }
}

/**
 * @brief Validates that altitude values cover the dense internal-node slot domain.
 *
 * @param tree Tree topology used by the operation.
 * @param altitude Altitude data indexed by node identifier.
 */
template <AltitudeValue T> inline void validateAttributeDeltaAltitudeBufferShape(const MorphologicalTree& tree, std::span<const T> altitude) {
    MMCFILTERS_CONTRACT_REQUIRE(altitude.size() == static_cast<size_t>(tree.getNumInternalNodeSlots()),
                                throw std::runtime_error("Altitude buffer size must match the dense internal-node domain."));
}

/**
 * @brief Reads a node altitude from the dense internal-node slot domain.
 *
 * @param altitude Altitude data indexed by node identifier.
 * @param nodeId Identifier of the node used by the operation.
 * @return The requested node altitude from the dense internal-node slot domain.
 */
template <AltitudeValue T> inline T attributeDeltaAltitudeAt(std::span<const T> altitude, NodeId nodeId) {
    MMCFILTERS_CONTRACT_REQUIRE(nodeId >= 0 && static_cast<size_t>(nodeId) < altitude.size(),
                                throw std::invalid_argument("Altitude access requires a valid internal NodeId."));
    return altitude[static_cast<size_t>(nodeId)];
}

/**
 * @brief Rejects invalid altitude deltas before ancestor traversal.
 *
 * @param delta Delta offset or radius used by the operation.
 * @param context Operation context or diagnostic label.
 */
template <AltitudeValue T> inline void validateAltitudeDelta(AltitudeDiff<T> delta, const char* context) {
    if constexpr (std::is_floating_point_v<AltitudeDiff<T>>) {
        MMCFILTERS_CONTRACT_REQUIRE(std::isfinite(delta),
                                    throw std::invalid_argument(std::string(context) + " requires a finite altitude delta."));
    }
    MMCFILTERS_CONTRACT_REQUIRE(delta >= AltitudeDiff<T>{},
                                throw std::invalid_argument(std::string(context) + " requires a non-negative altitude delta."));
}

namespace kernel {

/**
 * @brief Validation-free altitude ascendant search over established domains.
 * @param tree Established tree topology.
 * @param altitude Established altitude span.
 * @param nodeId Established live node.
 * @param delta Non-negative altitude distance.
 * @return First ancestor reaching the distance, with root padding.
 */
template <AltitudeValue T>
inline NodeId findAscendantByAltitudeDelta(const MorphologicalTree& tree, std::span<const T> altitude, NodeId nodeId, AltitudeDiff<T> delta) {
    const AltitudeDiff<T> nodeAltitude = static_cast<AltitudeDiff<T>>(altitude[static_cast<size_t>(nodeId)]);
    NodeId currentNodeId = nodeId;
    while (true) {
        const AltitudeDiff<T> currentAltitude = static_cast<AltitudeDiff<T>>(altitude[static_cast<size_t>(currentNodeId)]);
        if (tree.getAltitudeOrder() == AltitudeOrder::INCREASING_FROM_ROOT) {
            if (nodeAltitude - currentAltitude >= delta) {
                return currentNodeId;
            }
        } else if (currentAltitude - nodeAltitude >= delta) {
            return currentNodeId;
        }
        if (currentNodeId == tree.getRoot()) {
            return currentNodeId;
        }
        currentNodeId = CommittedTreeAccess::nodeParent(tree, currentNodeId);
    }
}

/**
 * @brief Validation-free largest-area update over established node domains.
 * @param areaByNode Established area values.
 * @param descendants Descendant selection updated in place.
 * @param ascendantNodeId Established ascendant slot to update.
 * @param candidateNodeId Established descendant candidate.
 */
inline void updateLargestAreaAttributeDeltaDescendant(const std::vector<int32_t>& areaByNode, std::vector<NodeId>& descendants,
                                                      NodeId ascendantNodeId, NodeId candidateNodeId) {
    NodeId& currentNodeId = descendants[static_cast<size_t>(ascendantNodeId)];
    if (ascendantNodeId == candidateNodeId) {
        return;
    }
    if (currentNodeId == InvalidNode || areaByNode[static_cast<size_t>(candidateNodeId)] > areaByNode[static_cast<size_t>(currentNodeId)] ||
        (areaByNode[static_cast<size_t>(candidateNodeId)] == areaByNode[static_cast<size_t>(currentNodeId)] && candidateNodeId < currentNodeId)) {
        currentNodeId = candidateNodeId;
    }
}

} // namespace kernel

/**
 * @brief Finds the first ancestor at altitude distance `delta` from `nodeId`.
 *
 * @details
 * The polarity follows the component-tree semantics: max-trees climb toward
 * smaller/equal gray levels, while min-trees climb toward larger/equal gray
 * levels. If no ancestor reaches the requested distance before the root, the
 * root is returned, matching the historical delta-attribute behavior.
 *
 * @param tree Tree topology used by the operation.
 * @param altitude Altitude data indexed by node identifier.
 * @param nodeId Identifier of the node used by the operation.
 * @param delta Delta offset or radius used by the operation.
 * @return The located first ancestor at altitude distance delta from nodeId.
 */
template <AltitudeValue T>
inline NodeId findAscendantByAltitudeDelta(const MorphologicalTree& tree, std::span<const T> altitude, NodeId nodeId, AltitudeDiff<T> delta) {
    validateAttributeDeltaAltitudeBufferShape(tree, altitude);
    validateGlobalMonotoneAltitudeOrder(tree, "findAscendantByAltitudeDelta");
    validateAltitudeDelta<T>(delta, "findAscendantByAltitudeDelta");
    MMCFILTERS_CONTRACT_REQUIRE(tree.isAlive(nodeId), throw std::invalid_argument("Node ascendant search requires a live internal NodeId."));
    return kernel::findAscendantByAltitudeDelta(tree, altitude, nodeId, delta);
}

/**
 * @brief Computes delta ascendants and largest-area descendants by altitude distance.
 *
 * @details
 * The returned vectors are indexed by dense internal `NodeId`. `first[node]`
 * stores the altitude-delta ascendant selected for `node`, and `second[asc]`
 * stores the largest-area descendant assigned to that ascendant. This helper is
 * shared by generic delta attributes and MSER stability computation; it is not a
 * weighted-tree structural operation.
 *
 * @param tree Tree topology used by the operation.
 * @param altitude Altitude data indexed by node identifier.
 * @param delta Delta offset or radius used by the operation.
 * @return The computed delta ascendants and largest-area descendants by altitude distance.
 */
template <AltitudeValue T>
inline std::pair<std::vector<NodeId>, std::vector<NodeId>> computeAscendantsAndDescendantsByAltitude(const MorphologicalTree& tree, std::span<const T> altitude,
                                                                                                     AltitudeDiff<T> delta) {
    validateAttributeDeltaAltitudeBufferShape(tree, altitude);
    validateGlobalMonotoneAltitudeOrder(tree, "computeAscendantsAndDescendantsByAltitude");
    validateAltitudeDelta<T>(delta, "computeAscendantsAndDescendantsByAltitude");
    std::vector<NodeId> ascendants(static_cast<size_t>(tree.getNumInternalNodeSlots()), InvalidNode);
    std::vector<NodeId> descendants(static_cast<size_t>(tree.getNumInternalNodeSlots()), InvalidNode);
    const std::vector<int32_t> areaByNode = computeAttributeDeltaAreasIncrementally(tree);

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        const NodeId ascendantNodeId = kernel::findAscendantByAltitudeDelta(tree, altitude, nodeId, delta);
        if (ascendantNodeId == InvalidNode) {
            continue;
        }
        kernel::updateLargestAreaAttributeDeltaDescendant(areaByNode, descendants, ascendantNodeId, nodeId);
        if (descendants[static_cast<size_t>(ascendantNodeId)] != InvalidNode) {
            ascendants[static_cast<size_t>(nodeId)] = ascendantNodeId;
        }
    }

    return {std::move(ascendants), std::move(descendants)};
}

} // namespace mmcfilters::detail
