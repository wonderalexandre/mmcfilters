#pragma once

#include "TreeAltitudeDeltaNeighborhood.hpp"
#include "TreeKindValidation.hpp"
#include "../MorphologicalTree.hpp"
#include "../../utils/Altitude.hpp"
#include "../../utils/Common.hpp"

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
    std::vector<NodeId> ascendants;
    std::vector<NodeId> descendants;
};

/**
 * @brief Rejects zero or negative topological depth deltas.
 */
inline void validatePositiveDepthDelta(int depthDelta, const char* context) {
    if (depthDelta <= 0) {
        throw std::invalid_argument(std::string(context) + " requires a positive depth delta.");
    }
}

/**
 * @brief Rejects zero, negative, or non-finite altitude deltas.
 */
template<AltitudeValue T>
inline void validatePositiveAltitudeDelta(AltitudeDiff<T> delta, const char* context) {
    if constexpr (std::is_floating_point_v<AltitudeDiff<T>>) {
        if (!std::isfinite(delta)) {
            throw std::invalid_argument(std::string(context) + " requires a finite altitude delta.");
        }
    }
    if (delta <= AltitudeDiff<T>{}) {
        throw std::invalid_argument(std::string(context) + " requires a positive altitude delta.");
    }
}

/**
 * @brief Deterministically keeps the largest-area candidate, then smallest node id.
 */
inline void updateBestAreaCandidate(
    const std::vector<int32_t>& areaByNode,
    NodeId& bestNode,
    NodeId candidateNode) {
    if (candidateNode == InvalidNode) {
        return;
    }
    if (bestNode == InvalidNode ||
        areaByNode[static_cast<std::size_t>(candidateNode)] > areaByNode[static_cast<std::size_t>(bestNode)] ||
        (areaByNode[static_cast<std::size_t>(candidateNode)] == areaByNode[static_cast<std::size_t>(bestNode)] &&
         candidateNode < bestNode)) {
        bestNode = candidateNode;
    }
}

/**
 * @brief Finds the first ancestor that reaches a positive altitude distance.
 *
 * Unlike the generic delta-attribute helper, this returns `InvalidNode` when the
 * root is reached before the requested distance is available. MSER stability
 * needs a complete window instead of padding missing context with the root.
 */
template<AltitudeValue T>
inline NodeId findStrictAscendantByAltitudeDelta(
    const MorphologicalTree& tree,
    std::span<const T> altitude,
    NodeId nodeId,
    AltitudeDiff<T> delta) {
    validateAttributeDeltaAltitudeBufferShape(tree, altitude);
    validateComponentTreeKind(tree, "findStrictAscendantByAltitudeDelta");
    validatePositiveAltitudeDelta<T>(delta, "findStrictAscendantByAltitudeDelta");
    if (!tree.isAlive(nodeId)) {
        throw std::invalid_argument("Strict altitude ascendant search requires a live internal NodeId.");
    }

    const AltitudeDiff<T> nodeAltitude = static_cast<AltitudeDiff<T>>(attributeDeltaAltitudeAt(altitude, nodeId));
    NodeId currentNodeId = nodeId;
    while (!tree.isRoot(currentNodeId)) {
        currentNodeId = tree.getNodeParent(currentNodeId);
        if (currentNodeId == InvalidNode || !tree.isAlive(currentNodeId)) {
            return InvalidNode;
        }

        const AltitudeDiff<T> currentAltitude = static_cast<AltitudeDiff<T>>(attributeDeltaAltitudeAt(altitude, currentNodeId));
        if (tree.getTreeType() == MorphologicalTreeKind::MAX_TREE) {
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
 * @brief Computes strict altitude-delta windows for classical MSER on component trees.
 *
 * `ascendants[node]` is the first ancestor whose altitude is at least `delta`
 * away from `node` in the component-tree polarity. `descendants[asc]` is the
 * largest-area node whose first such ancestor is `asc`, with ties resolved by
 * the smallest `NodeId`.
 */
template<AltitudeValue T>
inline StabilityNeighborhood computeAltitudeStabilityNeighborhood(
    const MorphologicalTree& tree,
    std::span<const T> altitude,
    AltitudeDiff<T> delta) {
    validateAttributeDeltaAltitudeBufferShape(tree, altitude);
    validateComponentTreeKind(tree, "computeAltitudeStabilityNeighborhood");
    validatePositiveAltitudeDelta<T>(delta, "computeAltitudeStabilityNeighborhood");

    StabilityNeighborhood neighborhood{
        std::vector<NodeId>(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), InvalidNode),
        std::vector<NodeId>(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), InvalidNode)};
    const std::vector<int32_t> areaByNode = computeAttributeDeltaAreasIncrementally(tree);

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        const NodeId ascendantNodeId =
            findStrictAscendantByAltitudeDelta(tree, altitude, nodeId, delta);
        if (ascendantNodeId == InvalidNode) {
            continue;
        }
        neighborhood.ascendants[static_cast<std::size_t>(nodeId)] = ascendantNodeId;
        updateBestAreaCandidate(
            areaByNode,
            neighborhood.descendants[static_cast<std::size_t>(ascendantNodeId)],
            nodeId);
    }

    return neighborhood;
}

/**
 * @brief Finds the ancestor exactly `depthDelta` parent links above `nodeId`.
 */
inline NodeId findAscendantByDepthDelta(
    const MorphologicalTree& tree,
    NodeId nodeId,
    int depthDelta) {
    validatePositiveDepthDelta(depthDelta, "findAscendantByDepthDelta");
    if (!tree.isAlive(nodeId)) {
        throw std::invalid_argument("Depth ascendant search requires a live internal NodeId.");
    }

    NodeId currentNodeId = nodeId;
    for (int d = 0; d < depthDelta; ++d) {
        if (tree.isRoot(currentNodeId)) {
            return InvalidNode;
        }
        currentNodeId = tree.getNodeParent(currentNodeId);
        if (currentNodeId == InvalidNode || !tree.isAlive(currentNodeId)) {
            return InvalidNode;
        }
    }
    return currentNodeId;
}

/**
 * @brief Computes descendant representatives exactly `depthDelta` child links below each node.
 *
 * When several descendants exist at the requested depth, the representative is
 * the largest-area descendant, with ties resolved by the smallest `NodeId`.
 */
inline std::vector<NodeId> computeDepthDescendants(
    const MorphologicalTree& tree,
    int depthDelta) {
    validatePositiveDepthDelta(depthDelta, "computeDepthDescendants");
    const std::vector<int32_t> areaByNode = computeAttributeDeltaAreasIncrementally(tree);
    std::vector<NodeId> previous(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), InvalidNode);
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        previous[static_cast<std::size_t>(nodeId)] = nodeId;
    }

    std::vector<NodeId> current(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), InvalidNode);
    for (int depth = 1; depth <= depthDelta; ++depth) {
        std::fill(current.begin(), current.end(), InvalidNode);
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            NodeId bestNode = InvalidNode;
            for (NodeId childId : tree.getChildren(nodeId)) {
                updateBestAreaCandidate(
                    areaByNode,
                    bestNode,
                    previous[static_cast<std::size_t>(childId)]);
            }
            current[static_cast<std::size_t>(nodeId)] = bestNode;
        }
        previous.swap(current);
    }

    return previous;
}

/**
 * @brief Computes edge-count windows for self-dual/topological stability.
 *
 * `depthDelta` means exactly that many tree edges: an ancestor is selected by
 * climbing `depthDelta` parent links, and a descendant is selected among nodes
 * exactly `depthDelta` child links below the center. Altitude is not read.
 */
inline StabilityNeighborhood computeDepthStabilityNeighborhood(
    const MorphologicalTree& tree,
    int depthDelta) {
    validatePositiveDepthDelta(depthDelta, "computeDepthStabilityNeighborhood");

    StabilityNeighborhood neighborhood{
        std::vector<NodeId>(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), InvalidNode),
        computeDepthDescendants(tree, depthDelta)};

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        neighborhood.ascendants[static_cast<std::size_t>(nodeId)] =
            findAscendantByDepthDelta(tree, nodeId, depthDelta);
    }

    return neighborhood;
}

} // namespace mmcfilters::detail
