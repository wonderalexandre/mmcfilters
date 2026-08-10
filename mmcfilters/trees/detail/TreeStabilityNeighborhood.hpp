#pragma once

#include "TreeAltitudeDeltaNeighborhood.hpp"
#include "HierarchyCapabilityValidation.hpp"
#include "../MorphologicalTree.hpp"
#include "CommittedTreeAccess.hpp"
#include "../../utils/Altitude.hpp"
#include "../../utils/Common.hpp"
#include "../../utils/Contract.hpp"

#include <algorithm>
#include <cmath>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace mmcfilters::detail {

/**
 * @brief Dense ancestor/descendant window used by stability criteria.
 *
 * Both buffers are indexed by internal `NodeId` slots. `ascendants[node]` stores
 * the selected context node above `node`; `descendants[node]` stores the selected
 * context node below `node`. Missing context is represented by `InvalidNode`.
 */
struct StabilityNeighborhood {
    /** @brief Stores the ascendants. */
    std::vector<NodeId> ascendants;
    /** @brief Stores the descendants. */
    std::vector<NodeId> descendants;
};

/**
 * @brief Rejects zero or negative topological depth deltas.
 *
 * @param depthDelta Topological depth offset used by the operation.
 * @param context Operation context or diagnostic label.
 */
inline void validatePositiveDepthDelta(int depthDelta, const char* context) {
    MMCFILTERS_CONTRACT_REQUIRE(depthDelta > 0,
                                throw std::invalid_argument(std::string(context) + " requires a positive depth delta."));
}

/**
 * @brief Rejects zero, negative, or non-finite altitude deltas.
 *
 * @param delta Delta offset or radius used by the operation.
 * @param context Operation context or diagnostic label.
 */
template <AltitudeValue T> inline void validatePositiveAltitudeDelta(AltitudeDiff<T> delta, const char* context) {
    if constexpr (std::is_floating_point_v<AltitudeDiff<T>>) {
        MMCFILTERS_CONTRACT_REQUIRE(std::isfinite(delta),
                                    throw std::invalid_argument(std::string(context) + " requires a finite altitude delta."));
    }
    MMCFILTERS_CONTRACT_REQUIRE(delta > AltitudeDiff<T>{},
                                throw std::invalid_argument(std::string(context) + " requires a positive altitude delta."));
}

/**
 * @brief Deterministically keeps the largest-area candidate, then smallest node id.
 *
 * @param areaByNode Node identifier represented by `areaByNode`.
 * @param bestNode Node identifier represented by `bestNode`.
 * @param candidateNode Node identifier represented by `candidateNode`.
 */
inline void updateBestAreaCandidate(const std::vector<int32_t>& areaByNode, NodeId& bestNode, NodeId candidateNode) {
    if (candidateNode == InvalidNode) {
        return;
    }
    if (bestNode == InvalidNode || areaByNode[static_cast<std::size_t>(candidateNode)] > areaByNode[static_cast<std::size_t>(bestNode)] ||
        (areaByNode[static_cast<std::size_t>(candidateNode)] == areaByNode[static_cast<std::size_t>(bestNode)] && candidateNode < bestNode)) {
        bestNode = candidateNode;
    }
}

namespace kernel {

/**
 * @brief Validation-free altitude ascendant search over an established tree and altitude domain.
 * @param tree Established tree topology.
 * @param altitude Established altitude span.
 * @param nodeId Established live node.
 * @param delta Positive altitude distance.
 * @return First strict ancestor reaching the distance, or InvalidNode.
 */
template <AltitudeValue T>
inline NodeId findStrictAscendantByAltitudeDelta(const MorphologicalTree& tree, std::span<const T> altitude, NodeId nodeId,
                                                 AltitudeDiff<T> delta) {
    const AltitudeDiff<T> nodeAltitude = static_cast<AltitudeDiff<T>>(altitude[static_cast<std::size_t>(nodeId)]);
    NodeId currentNodeId = nodeId;
    const NodeId root = tree.getRoot();
    while (currentNodeId != root) {
        currentNodeId = CommittedTreeAccess::nodeParent(tree, currentNodeId);
        if (currentNodeId == InvalidNode) {
            return InvalidNode;
        }

        const AltitudeDiff<T> currentAltitude = static_cast<AltitudeDiff<T>>(altitude[static_cast<std::size_t>(currentNodeId)]);
        if (tree.getAltitudeOrder() == AltitudeOrder::INCREASING_FROM_ROOT) {
            if (nodeAltitude - currentAltitude >= delta) {
                return currentNodeId;
            }
        } else if (currentAltitude - nodeAltitude >= delta) {
            return currentNodeId;
        }
    }
    return InvalidNode;
}

/**
 * @brief Validation-free depth ascendant search over an established live-node domain.
 * @param tree Established tree topology.
 * @param nodeId Established live node.
 * @param depthDelta Positive edge distance.
 * @return Ancestor at the requested distance, or InvalidNode.
 */
inline NodeId findAscendantByDepthDelta(const MorphologicalTree& tree, NodeId nodeId, int depthDelta) {
    NodeId currentNodeId = nodeId;
    for (int depth = 0; depth < depthDelta; ++depth) {
        if (currentNodeId == tree.getRoot()) {
            return InvalidNode;
        }
        currentNodeId = CommittedTreeAccess::nodeParent(tree, currentNodeId);
        if (currentNodeId == InvalidNode) {
            return InvalidNode;
        }
    }
    return currentNodeId;
}

/**
 * @brief Validation-free depth-descendant computation over an established tree domain.
 * @param tree Established tree topology.
 * @param depthDelta Positive edge distance.
 * @param areaByNode Established area values used for deterministic selection.
 * @return Largest-area descendant at the requested depth for every node.
 */
inline std::vector<NodeId> computeDepthDescendants(const MorphologicalTree& tree, int depthDelta, const std::vector<int32_t>& areaByNode) {
    std::vector<NodeId> previous(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), InvalidNode);
    const NodeId numSlots = tree.getNumInternalNodeSlots();
    for (NodeId nodeId = 0; nodeId < numSlots; ++nodeId) {
        if (CommittedTreeAccess::isAlive(tree, nodeId)) {
            previous[static_cast<std::size_t>(nodeId)] = nodeId;
        }
    }

    std::vector<NodeId> current(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), InvalidNode);
    for (int depth = 1; depth <= depthDelta; ++depth) {
        std::fill(current.begin(), current.end(), InvalidNode);
        for (NodeId nodeId = 0; nodeId < numSlots; ++nodeId) {
            if (!CommittedTreeAccess::isAlive(tree, nodeId)) {
                continue;
            }
            NodeId bestNode = InvalidNode;
            for (NodeId childId : CommittedTreeAccess::children(tree, nodeId)) {
                updateBestAreaCandidate(areaByNode, bestNode, previous[static_cast<std::size_t>(childId)]);
            }
            current[static_cast<std::size_t>(nodeId)] = bestNode;
        }
        previous.swap(current);
    }
    return previous;
}

} // namespace kernel

/**
 * @brief Finds the first ancestor that reaches a positive altitude distance.
 *
 * Unlike the generic delta-attribute helper, this returns `InvalidNode` when the
 * root is reached before the requested distance is available. MSER stability
 * needs a complete window instead of padding missing context with the root.
 *
 * @param tree Tree topology used by the operation.
 * @param altitude Altitude data indexed by node identifier.
 * @param nodeId Identifier of the node used by the operation.
 * @param delta Delta offset or radius used by the operation.
 * @return The located first ancestor that reaches a positive altitude distance.
 */
template <AltitudeValue T>
inline NodeId findStrictAscendantByAltitudeDelta(const MorphologicalTree& tree, std::span<const T> altitude, NodeId nodeId, AltitudeDiff<T> delta) {
    validateAttributeDeltaAltitudeBufferShape(tree, altitude);
    validateGlobalMonotoneAltitudeOrder(tree, "findStrictAscendantByAltitudeDelta");
    validatePositiveAltitudeDelta<T>(delta, "findStrictAscendantByAltitudeDelta");
    MMCFILTERS_CONTRACT_REQUIRE(tree.isAlive(nodeId),
                                throw std::invalid_argument("Strict altitude ascendant search requires a live internal NodeId."));
    return kernel::findStrictAscendantByAltitudeDelta(tree, altitude, nodeId, delta);
}

/**
 * @brief Computes strict altitude-delta windows for classical MSER on component trees.
 *
 * `ascendants[node]` is the first ancestor whose altitude is at least `delta`
 * away from `node` in the component-tree polarity. `descendants[asc]` is the
 * largest-area node whose first such ancestor is `asc`, with ties resolved by
 * the smallest `NodeId`.
 *
 * @param tree Tree topology used by the operation.
 * @param altitude Altitude data indexed by node identifier.
 * @param delta Delta offset or radius used by the operation.
 * @return The computed strict altitude-delta windows for classical MSER on component trees.
 */
template <AltitudeValue T>
inline StabilityNeighborhood computeAltitudeStabilityNeighborhood(const MorphologicalTree& tree, std::span<const T> altitude, AltitudeDiff<T> delta) {
    validateAttributeDeltaAltitudeBufferShape(tree, altitude);
    validateGlobalMonotoneAltitudeOrder(tree, "computeAltitudeStabilityNeighborhood");
    validatePositiveAltitudeDelta<T>(delta, "computeAltitudeStabilityNeighborhood");

    StabilityNeighborhood neighborhood{std::vector<NodeId>(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), InvalidNode),
                                       std::vector<NodeId>(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), InvalidNode)};
    const std::vector<int32_t> areaByNode = computeAttributeDeltaAreasIncrementally(tree);

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        const NodeId ascendantNodeId = kernel::findStrictAscendantByAltitudeDelta(tree, altitude, nodeId, delta);
        if (ascendantNodeId == InvalidNode) {
            continue;
        }
        neighborhood.ascendants[static_cast<std::size_t>(nodeId)] = ascendantNodeId;
        updateBestAreaCandidate(areaByNode, neighborhood.descendants[static_cast<std::size_t>(ascendantNodeId)], nodeId);
    }

    return neighborhood;
}

/**
 * @brief Finds the ancestor exactly `depthDelta` parent links above `nodeId`.
 *
 * @param tree Tree topology used by the operation.
 * @param nodeId Identifier of the node used by the operation.
 * @param depthDelta Topological depth offset used by the operation.
 * @return The located ancestor exactly depthDelta parent links above nodeId.
 */
inline NodeId findAscendantByDepthDelta(const MorphologicalTree& tree, NodeId nodeId, int depthDelta) {
    validatePositiveDepthDelta(depthDelta, "findAscendantByDepthDelta");
    MMCFILTERS_CONTRACT_REQUIRE(tree.isAlive(nodeId), throw std::invalid_argument("Depth ascendant search requires a live internal NodeId."));
    return kernel::findAscendantByDepthDelta(tree, nodeId, depthDelta);
}

/**
 * @brief Computes descendant representatives exactly `depthDelta` child links below each node.
 *
 * When several descendants exist at the requested depth, the representative is
 * the largest-area descendant, with ties resolved by the smallest `NodeId`.
 *
 * @param tree Tree topology used by the operation.
 * @param depthDelta Topological depth offset used by the operation.
 * @return The computed descendant representatives exactly depthDelta child links below each node.
 */
inline std::vector<NodeId> computeDepthDescendants(const MorphologicalTree& tree, int depthDelta) {
    validatePositiveDepthDelta(depthDelta, "computeDepthDescendants");
    const std::vector<int32_t> areaByNode = computeAttributeDeltaAreasIncrementally(tree);
    return kernel::computeDepthDescendants(tree, depthDelta, areaByNode);
}

/**
 * @brief Computes edge-count windows for self-dual/topological stability.
 *
 * `depthDelta` means exactly that many tree edges: an ancestor is selected by
 * climbing `depthDelta` parent links, and a descendant is selected among nodes
 * exactly `depthDelta` child links below the center. Altitude is not read.
 *
 * @param tree Tree topology used by the operation.
 * @param depthDelta Topological depth offset used by the operation.
 * @return The computed edge-count windows for self-dual/topological stability.
 */
inline StabilityNeighborhood computeDepthStabilityNeighborhood(const MorphologicalTree& tree, int depthDelta) {
    validatePositiveDepthDelta(depthDelta, "computeDepthStabilityNeighborhood");

    const std::vector<int32_t> areaByNode = computeAttributeDeltaAreasIncrementally(tree);

    StabilityNeighborhood neighborhood{std::vector<NodeId>(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), InvalidNode),
                                       kernel::computeDepthDescendants(tree, depthDelta, areaByNode)};

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        neighborhood.ascendants[static_cast<std::size_t>(nodeId)] = kernel::findAscendantByDepthDelta(tree, nodeId, depthDelta);
    }

    return neighborhood;
}

} // namespace mmcfilters::detail
