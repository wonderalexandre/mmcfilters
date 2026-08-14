#pragma once

#include "../MorphologicalTree.hpp"
#include "CommittedTreeAccess.hpp"
#include "HierarchyCapabilityValidation.hpp"
#include "../../attributes/AttributeTypes.hpp"
#include "../../utils/Altitude.hpp"
#include "../../utils/Common.hpp"
#include "../../utils/Contract.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace mmcfilters::detail {

/**
 * @brief Dense ancestor and representative-descendant samples for one altitude distance.
 */
struct NodeAttributeSamplingNeighborhood {
    /// First ancestors that reach the requested altitude distance.
    std::vector<NodeId> ancestors;

    /// Descendants selected by the requested sampling policy.
    std::vector<NodeId> representativeDescendants;
};

/**
 * @brief Support metadata used by representative-descendant selection.
 */
struct NodeSupportSamplingMetadata {
    /// Number of pixels in every live node support.
    std::span<const std::int32_t> cardinalityByNode;

    /// Smallest row-major pixel identifier in every live node support.
    std::span<const PixelId> smallestPixelByNode;
};

/**
 * @brief Returns cached node-support cardinalities and spatial tie-break keys.
 */
inline NodeSupportSamplingMetadata computeNodeSupportSamplingMetadata(const MorphologicalTree& tree) {
    return {CommittedTreeAccess::nodeSupportCardinalities(tree), CommittedTreeAccess::smallestNodeSupportPixels(tree)};
}

/**
 * @brief Selects the greatest-support descendant, breaking ties spatially.
 */
struct LargestSupportDescendantPolicy {
    /**
     * @brief Updates one representative with a valid descendant candidate.
     * @param metadata Node-support comparison metadata.
     * @param representative Current representative, updated when appropriate.
     * @param candidate Candidate descendant.
     */
    static void update(const NodeSupportSamplingMetadata& metadata, NodeId& representative, NodeId candidate) {
        if (candidate == InvalidNode) {
            return;
        }

        const std::size_t candidateIndex = static_cast<std::size_t>(candidate);
        if (metadata.smallestPixelByNode[candidateIndex] == std::numeric_limits<PixelId>::max()) {
            throw std::logic_error("Node-attribute sampling requires every candidate to have a non-empty support.");
        }

        if (representative == InvalidNode) {
            representative = candidate;
            return;
        }

        const std::size_t representativeIndex = static_cast<std::size_t>(representative);
        const std::int32_t candidateCardinality = metadata.cardinalityByNode[candidateIndex];
        const std::int32_t representativeCardinality = metadata.cardinalityByNode[representativeIndex];
        if (candidateCardinality > representativeCardinality ||
            (candidateCardinality == representativeCardinality &&
             metadata.smallestPixelByNode[candidateIndex] < metadata.smallestPixelByNode[representativeIndex])) {
            representative = candidate;
        }
    }
};

/**
 * @brief Validates the dense altitude buffer used by node-attribute sampling.
 */
template <AltitudeValue T>
inline void validateNodeAttributeSamplingAltitudeBufferShape(const MorphologicalTree& tree, std::span<const T> altitude) {
    MMCFILTERS_CONTRACT_REQUIRE(altitude.size() == static_cast<std::size_t>(tree.numInternalNodeSlots()),
                                throw std::runtime_error("Altitude buffer size must match the dense internal-node domain."));
}

/**
 * @brief Rejects a zero, negative, or non-finite altitude step.
 */
template <AltitudeValue T> inline void validatePositiveAltitudeStep(AltitudeDifference<T> altitudeStep, const char* context) {
    if constexpr (std::is_floating_point_v<AltitudeDifference<T>>) {
        MMCFILTERS_CONTRACT_REQUIRE(std::isfinite(altitudeStep),
                                    throw std::invalid_argument(std::string(context) + " requires a finite altitude step."));
    }
    MMCFILTERS_CONTRACT_REQUIRE(altitudeStep > AltitudeDifference<T>{},
                                throw std::invalid_argument(std::string(context) + " requires a positive altitude step."));
}

/** @brief Rejects unsupported representative-descendant sampling policies. */
inline void validateNodeAttributeSamplingPolicy(NodeAttributeSamplingPolicy samplingPolicy) {
    switch (samplingPolicy) {
    case NodeAttributeSamplingPolicy::LargestSupportDescendant:
        return;
    default:
        throw std::invalid_argument("Unknown node-attribute sampling policy.");
    }
}

namespace kernel {

/**
 * @brief Finds the first strict ancestor that reaches an established altitude distance.
 * @tparam T Node-altitude type.
 * @param tree Established tree topology.
 * @param altitude Dense node-altitude buffer.
 * @param nodeId Starting node.
 * @param altitudeDistance Required absolute altitude distance.
 * @return First qualifying ancestor, or `InvalidNode` when none exists.
 */
template <AltitudeValue T>
inline NodeId findAncestorByAltitudeDistance(const MorphologicalTree& tree, std::span<const T> altitude, NodeId nodeId,
                                             AltitudeDifference<T> altitudeDistance) {
    const AltitudeDifference<T> nodeAltitude = static_cast<AltitudeDifference<T>>(altitude[static_cast<std::size_t>(nodeId)]);
    NodeId currentNodeId = nodeId;
    while (currentNodeId != tree.root()) {
        currentNodeId = CommittedTreeAccess::nodeParent(tree, currentNodeId);
        if (currentNodeId == InvalidNode) {
            return InvalidNode;
        }

        const AltitudeDifference<T> currentAltitude =
            static_cast<AltitudeDifference<T>>(altitude[static_cast<std::size_t>(currentNodeId)]);
        if (tree.nodeAltitudeOrder() == NodeAltitudeOrder::Increasing) {
            if (nodeAltitude - currentAltitude >= altitudeDistance) {
                return currentNodeId;
            }
        } else if (currentAltitude - nodeAltitude >= altitudeDistance) {
            return currentNodeId;
        }
    }
    return InvalidNode;
}

/**
 * @brief Fills reusable ancestor and representative-descendant buffers for one altitude distance.
 * @param tree Established tree topology.
 * @param aliveNodes Established live-node range.
 * @param altitude Dense node-altitude buffer.
 * @param altitudeDistance Required absolute altitude distance.
 * @param metadata Reusable support metadata, independent of altitude distance.
 * @param neighborhood Reusable dense buffers to reset and fill.
 */
template <AltitudeValue T>
inline void fillNodeAttributeSamplingNeighborhood(const MorphologicalTree& tree, const MorphologicalTree::AliveNodeRange& aliveNodes,
                                                  std::span<const T> altitude, AltitudeDifference<T> altitudeDistance,
                                                  const NodeSupportSamplingMetadata& metadata,
                                                  NodeAttributeSamplingNeighborhood& neighborhood) {
    std::fill(neighborhood.ancestors.begin(), neighborhood.ancestors.end(), InvalidNode);
    std::fill(neighborhood.representativeDescendants.begin(), neighborhood.representativeDescendants.end(), InvalidNode);

    for (NodeId nodeId : aliveNodes) {
        const NodeId ancestorNodeId = findAncestorByAltitudeDistance(tree, altitude, nodeId, altitudeDistance);
        if (ancestorNodeId == InvalidNode) {
            continue;
        }

        neighborhood.ancestors[static_cast<std::size_t>(nodeId)] = ancestorNodeId;
        NodeId& representative = neighborhood.representativeDescendants[static_cast<std::size_t>(ancestorNodeId)];
        LargestSupportDescendantPolicy::update(metadata, representative, nodeId);
    }
}

} // namespace kernel

/**
 * @brief Finds the first ancestor at the requested positive altitude distance.
 *
 * @return The first qualifying ancestor, or `InvalidNode` when the root is
 * reached before the requested distance is available.
 */
template <AltitudeValue T>
inline NodeId findAncestorByAltitudeDistance(const MorphologicalTree& tree, std::span<const T> altitude, NodeId nodeId,
                                             AltitudeDifference<T> altitudeDistance) {
    validateNodeAttributeSamplingAltitudeBufferShape(tree, altitude);
    validateGlobalMonotoneAltitudeOrder(tree, "findAncestorByAltitudeDistance");
    validatePositiveAltitudeStep<T>(altitudeDistance, "findAncestorByAltitudeDistance");
    MMCFILTERS_CONTRACT_REQUIRE(tree.isAlive(nodeId),
                                throw std::invalid_argument("Altitude-based ancestor search requires a live internal NodeId."));
    return kernel::findAncestorByAltitudeDistance(tree, altitude, nodeId, altitudeDistance);
}

/**
 * @brief Computes ancestor and representative-descendant samples for one altitude distance.
 */
template <AltitudeValue T>
inline NodeAttributeSamplingNeighborhood
computeNodeAttributeSamplingNeighborhood(const MorphologicalTree& tree, std::span<const T> altitude,
                                         AltitudeDifference<T> altitudeDistance, NodeAttributeSamplingPolicy samplingPolicy) {
    validateNodeAttributeSamplingAltitudeBufferShape(tree, altitude);
    validateGlobalMonotoneAltitudeOrder(tree, "computeNodeAttributeSamplingNeighborhood");
    validatePositiveAltitudeStep<T>(altitudeDistance, "computeNodeAttributeSamplingNeighborhood");
    validateNodeAttributeSamplingPolicy(samplingPolicy);

    NodeAttributeSamplingNeighborhood neighborhood{
        std::vector<NodeId>(static_cast<std::size_t>(tree.numInternalNodeSlots()), InvalidNode),
        std::vector<NodeId>(static_cast<std::size_t>(tree.numInternalNodeSlots()), InvalidNode),
    };
    const NodeSupportSamplingMetadata metadata = computeNodeSupportSamplingMetadata(tree);
    const auto aliveNodes = tree.aliveNodeIds();

    kernel::fillNodeAttributeSamplingNeighborhood(tree, aliveNodes, altitude, altitudeDistance, metadata, neighborhood);

    return neighborhood;
}

} // namespace mmcfilters::detail
