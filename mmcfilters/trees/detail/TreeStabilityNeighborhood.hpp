#pragma once

#include "TreeAttributeSamplingNeighborhood.hpp"
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
 * @brief Dense ancestor/descendant window used by stability measures.
 *
 * Both buffers are indexed by internal `NodeId` slots. `ancestors[node]` stores
 * the selected context node above `node`; `descendants[node]` stores the selected
 * context node below `node`. Missing context is represented by `InvalidNode`.
 */
struct StabilityNeighborhood {
    /** @brief Dense node identifier of the ancestors. */
    std::vector<NodeId> ancestors;
    /** @brief Dense node identifier of the descendants. */
    std::vector<NodeId> descendants;
};

/**
 * @brief Rejects zero or negative topological depth window radii.
 *
 * @param depthWindowRadius Topological depth-window radius.
 * @param context Operation context or diagnostic label.
 */
inline void validatePositiveDepthWindowRadius(int depthWindowRadius, const char* context) {
    MMCFILTERS_CONTRACT_REQUIRE(depthWindowRadius > 0,
                                throw std::invalid_argument(std::string(context) + " requires a positive depth window radius."));
}

/**
 * @brief Rejects zero, negative, or non-finite altitude window radii.
 *
 * @param altitudeWindowRadius Positive altitude-window radius.
 * @param context Operation context or diagnostic label.
 */
template <AltitudeValue T> inline void validatePositiveAltitudeWindowRadius(AltitudeDifference<T> altitudeWindowRadius, const char* context) {
    if constexpr (std::is_floating_point_v<AltitudeDifference<T>>) {
        MMCFILTERS_CONTRACT_REQUIRE(std::isfinite(altitudeWindowRadius),
                                    throw std::invalid_argument(std::string(context) + " requires a finite altitude window radius."));
    }
    MMCFILTERS_CONTRACT_REQUIRE(altitudeWindowRadius > AltitudeDifference<T>{},
                                throw std::invalid_argument(std::string(context) + " requires a positive altitude window radius."));
}

namespace kernel {

/**
 * @brief Validation-free altitude ancestor search over an established tree and altitude domain.
 * @param tree Established tree topology.
 * @param altitude Established altitude span.
 * @param nodeId Established live node.
 * @param altitudeWindowRadius Positive altitude distance.
 * @return First strict ancestor reaching the distance, or InvalidNode.
 */
template <AltitudeValue T>
inline NodeId findStrictAncestorByAltitudeWindowRadius(const MorphologicalTree& tree, std::span<const T> altitude, NodeId nodeId,
                                                       AltitudeDifference<T> altitudeWindowRadius) {
    const AltitudeDifference<T> nodeAltitude = static_cast<AltitudeDifference<T>>(altitude[static_cast<std::size_t>(nodeId)]);
    NodeId currentNodeId = nodeId;
    const NodeId root = tree.root();
    while (currentNodeId != root) {
        currentNodeId = CommittedTreeAccess::nodeParent(tree, currentNodeId);
        if (currentNodeId == InvalidNode) {
            return InvalidNode;
        }

        const AltitudeDifference<T> currentAltitude = static_cast<AltitudeDifference<T>>(altitude[static_cast<std::size_t>(currentNodeId)]);
        if (tree.nodeAltitudeOrder() == NodeAltitudeOrder::Increasing) {
            if (nodeAltitude - currentAltitude >= altitudeWindowRadius) {
                return currentNodeId;
            }
        } else if (currentAltitude - nodeAltitude >= altitudeWindowRadius) {
            return currentNodeId;
        }
    }
    return InvalidNode;
}

/**
 * @brief Validation-free depth ancestor search over an established live-node domain.
 * @param tree Established tree topology.
 * @param nodeId Established live node.
 * @param depthWindowRadius Positive edge distance.
 * @return Ancestor at the requested distance, or InvalidNode.
 */
inline NodeId findAncestorByDepthWindowRadius(const MorphologicalTree& tree, NodeId nodeId, int depthWindowRadius) {
    NodeId currentNodeId = nodeId;
    for (int depth = 0; depth < depthWindowRadius; ++depth) {
        if (currentNodeId == tree.root()) {
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
 * @param depthWindowRadius Positive edge distance.
 * @param supportMetadata Cached support cardinalities and spatial tie-break keys.
 * @return Greatest-support descendant at the requested depth for every node.
 */
inline std::vector<NodeId> computeDepthDescendants(const MorphologicalTree& tree, int depthWindowRadius,
                                                   const NodeSupportSamplingMetadata& supportMetadata) {
    std::vector<NodeId> previous(static_cast<std::size_t>(tree.numInternalNodeSlots()), InvalidNode);
    const int numSlots = tree.numInternalNodeSlots();
    for (NodeId nodeId = 0; nodeId < numSlots; ++nodeId) {
        if (CommittedTreeAccess::isAlive(tree, nodeId)) {
            previous[static_cast<std::size_t>(nodeId)] = nodeId;
        }
    }

    std::vector<NodeId> current(static_cast<std::size_t>(tree.numInternalNodeSlots()), InvalidNode);
    for (int depth = 1; depth <= depthWindowRadius; ++depth) {
        std::fill(current.begin(), current.end(), InvalidNode);
        for (NodeId nodeId = 0; nodeId < numSlots; ++nodeId) {
            if (!CommittedTreeAccess::isAlive(tree, nodeId)) {
                continue;
            }
            NodeId bestNode = InvalidNode;
            for (NodeId childId : CommittedTreeAccess::children(tree, nodeId)) {
                LargestSupportDescendantPolicy::update(supportMetadata, bestNode, previous[static_cast<std::size_t>(childId)]);
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
 * Unlike padded attribute sampling, this returns `InvalidNode` when the
 * root is reached before the requested distance is available. MSER stability
 * needs a complete window instead of padding missing context with the root.
 *
 * @param tree Tree topology.
 * @param altitude Altitude data indexed by node identifier.
 * @param nodeId Dense internal node identifier.
 * @param altitudeWindowRadius Positive altitude-window radius.
 * @return The located first ancestor that reaches a positive altitude distance.
 */
template <AltitudeValue T>
inline NodeId findStrictAncestorByAltitudeWindowRadius(const MorphologicalTree& tree, std::span<const T> altitude, NodeId nodeId,
                                                       AltitudeDifference<T> altitudeWindowRadius) {
    validateNodeAttributeSamplingAltitudeBufferShape(tree, altitude);
    validateGlobalMonotoneAltitudeOrder(tree, "findStrictAncestorByAltitudeWindowRadius");
    validatePositiveAltitudeWindowRadius<T>(altitudeWindowRadius, "findStrictAncestorByAltitudeWindowRadius");
    MMCFILTERS_CONTRACT_REQUIRE(tree.isAlive(nodeId),
                                throw std::invalid_argument("Strict altitude ancestor search requires a live internal NodeId."));
    return kernel::findStrictAncestorByAltitudeWindowRadius(tree, altitude, nodeId, altitudeWindowRadius);
}

/**
 * @brief Computes strict altitude-distance windows for classical MSER on component trees.
 *
 * `ancestors[node]` is the first ancestor whose altitude is at least the window
 * radius away from `node` in the component-tree polarity. The descendant of an
 * ancestor is the greatest-support node whose first such ancestor is that node.
 * Equal support cardinalities are resolved by the smallest row-major pixel in
 * the candidate support, independently of `NodeId`.
 *
 * @param tree Tree topology.
 * @param altitude Altitude data indexed by node identifier.
 * @param altitudeWindowRadius Positive altitude-window radius.
 * @return The computed strict altitude-distance windows for classical MSER on component trees.
 */
template <AltitudeValue T>
inline StabilityNeighborhood computeAltitudeStabilityNeighborhood(const MorphologicalTree& tree, std::span<const T> altitude,
                                                                   AltitudeDifference<T> altitudeWindowRadius) {
    validateNodeAttributeSamplingAltitudeBufferShape(tree, altitude);
    validateGlobalMonotoneAltitudeOrder(tree, "computeAltitudeStabilityNeighborhood");
    validatePositiveAltitudeWindowRadius<T>(altitudeWindowRadius, "computeAltitudeStabilityNeighborhood");

    StabilityNeighborhood neighborhood{std::vector<NodeId>(static_cast<std::size_t>(tree.numInternalNodeSlots()), InvalidNode),
                                       std::vector<NodeId>(static_cast<std::size_t>(tree.numInternalNodeSlots()), InvalidNode)};
    const NodeSupportSamplingMetadata supportMetadata = computeNodeSupportSamplingMetadata(tree);

    for (NodeId nodeId : tree.aliveNodeIds()) {
        const NodeId ancestorNodeId = kernel::findStrictAncestorByAltitudeWindowRadius(tree, altitude, nodeId, altitudeWindowRadius);
        if (ancestorNodeId == InvalidNode) {
            continue;
        }
        neighborhood.ancestors[static_cast<std::size_t>(nodeId)] = ancestorNodeId;
        LargestSupportDescendantPolicy::update(supportMetadata, neighborhood.descendants[static_cast<std::size_t>(ancestorNodeId)], nodeId);
    }

    return neighborhood;
}

/**
 * @brief Finds the ancestor exactly `depthWindowRadius` parent links above `nodeId`.
 *
 * @param tree Tree topology.
 * @param nodeId Dense internal node identifier.
 * @param depthWindowRadius Topological depth-window radius.
 * @return The located ancestor at the requested depth-window radius.
 */
inline NodeId findAncestorByDepthWindowRadius(const MorphologicalTree& tree, NodeId nodeId, int depthWindowRadius) {
    validatePositiveDepthWindowRadius(depthWindowRadius, "findAncestorByDepthWindowRadius");
    MMCFILTERS_CONTRACT_REQUIRE(tree.isAlive(nodeId), throw std::invalid_argument("Depth ancestor search requires a live internal NodeId."));
    return kernel::findAncestorByDepthWindowRadius(tree, nodeId, depthWindowRadius);
}

/**
 * @brief Computes descendant representatives at the requested depth-window radius.
 *
 * When several descendants exist at the requested depth, the representative is
 * the greatest-support descendant. Equal support cardinalities are resolved by
 * the smallest row-major pixel in the candidate support, never by `NodeId`.
 *
 * @param tree Tree topology.
 * @param depthWindowRadius Topological depth-window radius.
 * @return The computed descendant representatives at the requested radius.
 */
inline std::vector<NodeId> computeDepthDescendants(const MorphologicalTree& tree, int depthWindowRadius) {
    validatePositiveDepthWindowRadius(depthWindowRadius, "computeDepthDescendants");
    const NodeSupportSamplingMetadata supportMetadata = computeNodeSupportSamplingMetadata(tree);
    return kernel::computeDepthDescendants(tree, depthWindowRadius, supportMetadata);
}

/**
 * @brief Computes edge-count windows for self-dual/topological stability.
 *
 * `depthWindowRadius` means exactly that many tree edges: an ancestor is
 * selected by climbing that many parent links, and a descendant is selected
 * among nodes exactly that many child links below the center. Altitude is not read.
 *
 * @param tree Tree topology.
 * @param depthWindowRadius Topological depth-window radius.
 * @return The computed edge-count windows for self-dual/topological stability.
 */
inline StabilityNeighborhood computeDepthStabilityNeighborhood(const MorphologicalTree& tree, int depthWindowRadius) {
    validatePositiveDepthWindowRadius(depthWindowRadius, "computeDepthStabilityNeighborhood");

    const NodeSupportSamplingMetadata supportMetadata = computeNodeSupportSamplingMetadata(tree);

    StabilityNeighborhood neighborhood{std::vector<NodeId>(static_cast<std::size_t>(tree.numInternalNodeSlots()), InvalidNode),
                                       kernel::computeDepthDescendants(tree, depthWindowRadius, supportMetadata)};

    for (NodeId nodeId : tree.aliveNodeIds()) {
        neighborhood.ancestors[static_cast<std::size_t>(nodeId)] =
            kernel::findAncestorByDepthWindowRadius(tree, nodeId, depthWindowRadius);
    }

    return neighborhood;
}

} // namespace mmcfilters::detail
