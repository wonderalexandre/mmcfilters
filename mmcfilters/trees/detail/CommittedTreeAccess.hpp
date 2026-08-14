#pragma once

#include "../MorphologicalTree.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>

namespace mmcfilters::detail {

/**
 * @brief Direct reads over a committed tree inside established kernels.
 *
 * Public tree accessors remain defensive.  This facade is deliberately kept
 * in `detail` and may only be used when a dominating public boundary or a
 * construction invariant already established the node/pixel domain.
 */
class CommittedTreeAccess {
  private:
    /**
     * @brief Materializes reusable support metadata when committed-tree versions changed.
     * @param tree Committed tree whose support metadata is required.
     */
    static void ensureNodeSupportMetadata(const MorphologicalTree& tree) {
        auto& cache = tree.nodeSupportMetadataCache_;
        if (cache.nodeStructureVersion == tree.nodeStructureVersion_ && cache.topologyVersion == tree.topologyVersion_ &&
            cache.properPartVersion == tree.properPartVersion_) {
            return;
        }

        constexpr PixelId MissingPixel = std::numeric_limits<PixelId>::max();
        const std::size_t numSlots = static_cast<std::size_t>(tree.numInternalNodeSlots());
        cache.cardinalityByNode.assign(numSlots, 0);
        cache.smallestPixelByNode.assign(numSlots, MissingPixel);

        for (PixelId pixel = 0; pixel < static_cast<PixelId>(tree.smallestNodeMap_.size()); ++pixel) {
            const NodeId smallestNodeId = tree.smallestNodeMap_[static_cast<std::size_t>(pixel)];
            const std::size_t smallestNodeIndex = static_cast<std::size_t>(smallestNodeId);
            if (cache.smallestPixelByNode[smallestNodeIndex] == MissingPixel) {
                cache.smallestPixelByNode[smallestNodeIndex] = pixel;
            }
        }

        for (NodeId nodeId : tree.postOrder()) {
            const std::size_t nodeIndex = static_cast<std::size_t>(nodeId);
            cache.cardinalityByNode[nodeIndex] += static_cast<std::int32_t>(tree.properPartCardinalityByNode_[nodeIndex]);

            const NodeId parentNodeId = tree.nodeParent_[nodeIndex];
            if (parentNodeId != InvalidNode && parentNodeId != nodeId) {
                const std::size_t parentIndex = static_cast<std::size_t>(parentNodeId);
                cache.cardinalityByNode[parentIndex] += cache.cardinalityByNode[nodeIndex];
                cache.smallestPixelByNode[parentIndex] =
                    std::min(cache.smallestPixelByNode[parentIndex], cache.smallestPixelByNode[nodeIndex]);
            }
        }

        cache.nodeStructureVersion = tree.nodeStructureVersion_;
        cache.topologyVersion = tree.topologyVersion_;
        cache.properPartVersion = tree.properPartVersion_;
    }

  public:
    /** @brief Returns the established grid domain. @param tree Committed tree. @return Stored grid domain. */
    [[nodiscard]] static const GridDomain2D& gridDomain2D(const MorphologicalTree& tree) noexcept { return *tree.gridDomain2D_; }

    /**
     * @brief Returns cached support cardinalities for the current committed tree.
     * @param tree Committed tree whose support cardinalities are required.
     * @return Dense cardinalities indexed by internal node slot.
     */
    [[nodiscard]] static std::span<const std::int32_t> nodeSupportCardinalities(const MorphologicalTree& tree) {
        ensureNodeSupportMetadata(tree);
        return tree.nodeSupportMetadataCache_.cardinalityByNode;
    }

    /**
     * @brief Returns cached smallest row-major support pixels for the current committed tree.
     * @param tree Committed tree whose spatial support keys are required.
     * @return Dense smallest-pixel keys indexed by internal node slot.
     */
    [[nodiscard]] static std::span<const PixelId> smallestNodeSupportPixels(const MorphologicalTree& tree) {
        ensureNodeSupportMetadata(tree);
        return tree.nodeSupportMetadataCache_.smallestPixelByNode;
    }

    /**
     * @brief Reads the smallest node of an established pixel.
     * @param tree Committed tree.
     * @param pixel Established pixel id.
     * @return Smallest node id.
     */
    [[nodiscard]] static NodeId smallestNodeMap(const MorphologicalTree& tree, PixelId pixel) noexcept {
        return tree.smallestNodeMap_[static_cast<std::size_t>(pixel)];
    }

    /** @brief Reads the parent of an established node. @param tree Committed tree. @param nodeId Established node id. @return Parent node id. */
    [[nodiscard]] static NodeId nodeParent(const MorphologicalTree& tree, NodeId nodeId) noexcept { return tree.nodeParent_[static_cast<std::size_t>(nodeId)]; }

    /** @brief Iterates children of an established node. @param tree Committed tree. @param nodeId Established node id. @return Child range. */
    [[nodiscard]] static MorphologicalTree::ChildrenRange children(const MorphologicalTree& tree, NodeId nodeId) {
        return MorphologicalTree::ChildrenRange(&tree, tree.firstChild_[static_cast<std::size_t>(nodeId)], tree.topologyVersion_);
    }

    /** @brief Iterates proper parts of an established node. @param tree Committed tree. @param nodeId Established node id. @return Proper-part range. */
    [[nodiscard]] static MorphologicalTree::ProperPartRange properParts(const MorphologicalTree& tree, NodeId nodeId) {
        return MorphologicalTree::ProperPartRange(&tree, tree.properHead_[static_cast<std::size_t>(nodeId)], tree.properPartVersion_);
    }

    /** @brief Counts children of an established node. @param tree Committed tree. @param nodeId Established node id. @return Child count. */
    [[nodiscard]] static int numChildren(const MorphologicalTree& tree, NodeId nodeId) noexcept {
        return tree.numChildrenByNode_[static_cast<std::size_t>(nodeId)];
    }

    /** @brief Counts proper parts of an established node. @param tree Committed tree. @param nodeId Established node id. @return Proper-part count. */
    [[nodiscard]] static int properPartCardinality(const MorphologicalTree& tree, NodeId nodeId) noexcept {
        return tree.properPartCardinalityByNode_[static_cast<std::size_t>(nodeId)];
    }

    /** @brief Tests liveness of an established node slot. @param tree Committed tree. @param nodeId Established node id. @return True for a live node. */
    [[nodiscard]] static bool isAlive(const MorphologicalTree& tree, NodeId nodeId) noexcept { return tree.alive_[static_cast<std::size_t>(nodeId)] != 0; }

    /** @brief Maps a dense node to its preserved Higra id. @param tree Committed tree. @param nodeId Established node id. @return Higra id or InvalidNode. */
    [[nodiscard]] static NodeId higraNodeId(const MorphologicalTree& tree, NodeId nodeId) noexcept {
        return tree.preservedExternalNodeIdOffset_ ? *tree.preservedExternalNodeIdOffset_ + nodeId : InvalidNode;
    }

    /** @brief Iterates an established subtree. @param tree Committed tree. @param root Established subtree root. @return Subtree node range. */
    [[nodiscard]] static MorphologicalTree::SubtreeNodeRange subtree(const MorphologicalTree& tree, NodeId root) {
        return MorphologicalTree::SubtreeNodeRange(&tree, root, tree.topologyVersion_);
    }

    /**
     * @brief Finds the lowest common ancestor of two established nodes.
     * @param tree Committed tree.
     * @param lhs First established node.
     * @param rhs Second established node.
     * @return Lowest common ancestor.
     */
    [[nodiscard]] static NodeId lowestCommonAncestor(const MorphologicalTree& tree, NodeId lhs, NodeId rhs) {
        if (lhs == rhs) {
            return lhs;
        }
        tree.ensureDfsIntervalCache();
        const auto isAncestorCached = [&](NodeId ancestor, NodeId node) {
            return tree.dfsIntervalCache_.entryIndex[static_cast<std::size_t>(ancestor)] <=
                       tree.dfsIntervalCache_.entryIndex[static_cast<std::size_t>(node)] &&
                   tree.dfsIntervalCache_.exitIndex[static_cast<std::size_t>(ancestor)] >=
                       tree.dfsIntervalCache_.exitIndex[static_cast<std::size_t>(node)];
        };
        if (isAncestorCached(lhs, rhs)) {
            return lhs;
        }
        if (isAncestorCached(rhs, lhs)) {
            return rhs;
        }
        return tree.ensureLcaCache().findLowestCommonAncestor(lhs, rhs);
    }

    /**
     * @brief Tests the inclusive ancestor relation between established nodes.
     * @param tree Committed tree.
     * @param ancestor Candidate ancestor.
     * @param node Candidate descendant.
     * @return True when ancestor equals or contains node.
     */
    [[nodiscard]] static bool isAncestor(const MorphologicalTree& tree, NodeId ancestor, NodeId node) {
        tree.ensureDfsIntervalCache();
        return tree.dfsIntervalCache_.entryIndex[static_cast<std::size_t>(ancestor)] <=
                   tree.dfsIntervalCache_.entryIndex[static_cast<std::size_t>(node)] &&
               tree.dfsIntervalCache_.exitIndex[static_cast<std::size_t>(ancestor)] >=
                   tree.dfsIntervalCache_.exitIndex[static_cast<std::size_t>(node)];
    }
};

} // namespace mmcfilters::detail
