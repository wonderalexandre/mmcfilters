#pragma once

#include "../MorphologicalTree.hpp"
#include "TreeKernelAccess.hpp"
#include "../../utils/Altitude.hpp"
#include "../../utils/Common.hpp"

#include <algorithm>
#include <cstdint>
#include <span>
#include <vector>

namespace mmcfilters::detail {

/**
 * @brief Computes reusable ancestor and representative-descendant neighborhoods.
 *
 * Stability and node-attribute sampling share these dense node-slot buffers.
 * Altitude neighborhoods use the first ancestor at the requested altitude
 * distance. Depth neighborhoods use exactly the requested number of tree edges.
 * Descendants maximize support cardinality, with ties resolved by the smallest
 * row-major support pixel. Missing neighbors are represented by InvalidNode.
 *
 * The caller supplies a committed rooted tree with nonempty node supports.
 * Altitude computations additionally require matching globally monotone finite
 * altitudes and a positive finite distance; depth computations require a
 * positive radius. The tree and its support metadata must remain unchanged
 * while this object is used. Recreate it after topology or proper-part edits.
 * Consumers may move out the result buffers before discarding the object.
 */
class TreeNeighborhood {
  private:
    /** @brief Borrowed committed topology used for all computations. */
    const MorphologicalTree& tree_;
    /** @brief Cached support cardinalities indexed by node slot. */
    std::span<const std::int32_t> cardinalityByNode_;
    /** @brief Cached minimum support pixels used for spatial tie-breaking. */
    std::span<const PixelId> smallestPixelByNode_;
    /** @brief Indicates whether a previous computation populated the buffers. */
    bool populated_ = false;

    /** @brief Prepares the existing buffers for another neighborhood computation. */
    void prepare() {
        if (populated_) {
            std::fill(ancestors.begin(), ancestors.end(), InvalidNode);
            std::fill(descendants.begin(), descendants.end(), InvalidNode);
        }
        populated_ = true;
    }

    /**
     * @brief Updates the representative with the greatest-support available candidate.
     * @param representative Current representative, updated in place.
     * @param candidate Live descendant, or InvalidNode when unavailable.
     */
    void updateRepresentative(NodeId& representative, NodeId candidate) const {
        if (candidate == InvalidNode) {
            return;
        }
        if (representative == InvalidNode) {
            representative = candidate;
            return;
        }

        const auto candidateIndex = static_cast<std::size_t>(candidate);
        const auto representativeIndex = static_cast<std::size_t>(representative);
        const auto candidateCardinality = cardinalityByNode_[candidateIndex];
        const auto representativeCardinality = cardinalityByNode_[representativeIndex];
        if (candidateCardinality > representativeCardinality ||
            (candidateCardinality == representativeCardinality &&
             smallestPixelByNode_[candidateIndex] < smallestPixelByNode_[representativeIndex])) {
            representative = candidate;
        }
    }

    /**
     * @brief Finds the first strict ancestor reaching an established altitude distance.
     * @tparam T Node-altitude type.
     * @param altitude Matching dense altitude buffer.
     * @param nodeId Starting live node.
     * @param distance Positive finite altitude distance.
     * @return First qualifying ancestor, or InvalidNode when none exists.
     */
    template <AltitudeValue T>
    NodeId findAncestorByAltitude(std::span<const T> altitude, NodeId nodeId, AltitudeDifference<T> distance) const {
        const auto nodeAltitude = static_cast<AltitudeDifference<T>>(altitude[static_cast<std::size_t>(nodeId)]);
        NodeId current = nodeId;
        const NodeId root = tree_.root();
        while (current != root) {
            current = TreeKernelAccess::nodeParent(tree_, current);
            const auto currentAltitude = static_cast<AltitudeDifference<T>>(altitude[static_cast<std::size_t>(current)]);
            if (tree_.nodeAltitudeOrder() == NodeAltitudeOrder::Increasing) {
                if (nodeAltitude - currentAltitude >= distance) {
                    return current;
                }
            } else if (currentAltitude - nodeAltitude >= distance) {
                return current;
            }
        }
        return InvalidNode;
    }

    /**
     * @brief Finds the ancestor at an established edge distance.
     * @param nodeId Starting live node.
     * @param radius Positive number of parent links.
     * @return Ancestor at the requested distance, or InvalidNode when unavailable.
     */
    NodeId findAncestorByDepth(NodeId nodeId, int radius) const {
        NodeId current = nodeId;
        for (int depth = 0; depth < radius; ++depth) {
            if (current == tree_.root()) {
                return InvalidNode;
            }
            current = TreeKernelAccess::nodeParent(tree_, current);
        }
        return current;
    }

  public:
    /** @brief Selected ancestors indexed by dense internal node slot. */
    std::vector<NodeId> ancestors;
    /** @brief Greatest-support representative descendants indexed by node slot. */
    std::vector<NodeId> descendants;

    /**
     * @brief Allocates result buffers and borrows current support metadata once.
     * @param tree Committed topology that outlives this object.
     */
    explicit TreeNeighborhood(const MorphologicalTree& tree)
        : tree_(tree), cardinalityByNode_(TreeKernelAccess::nodeSupportCardinalities(tree)),
          smallestPixelByNode_(TreeKernelAccess::smallestNodeSupportPixels(tree)),
          ancestors(static_cast<std::size_t>(tree.numInternalNodeSlots()), InvalidNode),
          descendants(static_cast<std::size_t>(tree.numInternalNodeSlots()), InvalidNode) {}

    /**
     * @brief Computes an altitude neighborhood using the existing result buffers.
     *
     * A node becomes a descendant candidate of its first qualifying ancestor.
     * Missing context remains InvalidNode; sampling policies handle it later.
     * @tparam T Node-altitude type.
     * @param altitude Matching globally monotone finite altitudes.
     * @param distance Positive finite altitude distance.
     */
    template <AltitudeValue T> void computeByAltitude(std::span<const T> altitude, AltitudeDifference<T> distance) {
        prepare();
        for (NodeId nodeId : tree_.aliveNodeIds()) {
            const NodeId ancestor = findAncestorByAltitude(altitude, nodeId, distance);
            if (ancestor == InvalidNode) {
                continue;
            }
            ancestors[static_cast<std::size_t>(nodeId)] = ancestor;
            updateRepresentative(descendants[static_cast<std::size_t>(ancestor)], nodeId);
        }
    }

    /**
     * @brief Computes an edge-count neighborhood without reading altitudes.
     *
     * Descendants are selected among nodes exactly radius edges below each
     * center. Two descendant buffers alternate over depth steps.
     * @param radius Positive number of tree edges.
     */
    void computeByDepth(int radius) {
        prepare();
        const int numSlots = tree_.numInternalNodeSlots();
        for (NodeId nodeId = 0; nodeId < numSlots; ++nodeId) {
            if (TreeKernelAccess::isAlive(tree_, nodeId)) {
                descendants[static_cast<std::size_t>(nodeId)] = nodeId;
            }
        }

        std::vector<NodeId> current(static_cast<std::size_t>(numSlots), InvalidNode);
        for (int depth = 1; depth <= radius; ++depth) {
            std::fill(current.begin(), current.end(), InvalidNode);
            for (NodeId nodeId = 0; nodeId < numSlots; ++nodeId) {
                if (!TreeKernelAccess::isAlive(tree_, nodeId)) {
                    continue;
                }
                NodeId best = InvalidNode;
                for (NodeId child : TreeKernelAccess::children(tree_, nodeId)) {
                    updateRepresentative(best, descendants[static_cast<std::size_t>(child)]);
                }
                current[static_cast<std::size_t>(nodeId)] = best;
            }
            descendants.swap(current);
        }

        for (NodeId nodeId : tree_.aliveNodeIds()) {
            ancestors[static_cast<std::size_t>(nodeId)] = findAncestorByDepth(nodeId, radius);
        }
    }
};

} // namespace mmcfilters::detail
