#pragma once

#include "../MorphologicalTree.hpp"

#include <cstddef>

namespace mmcfilters::detail {

/**
 * @brief Direct reads over a committed tree inside established kernels.
 *
 * Public tree accessors remain defensive.  This facade is deliberately kept
 * in `detail` and may only be used when a dominating public boundary or a
 * construction invariant already established the node/proper-part domain.
 */
class CommittedTreeAccess {
  public:
    /** @brief Returns the established grid domain. @param tree Committed tree. @return Stored grid domain. */
    [[nodiscard]] static const GridDomain2D& gridDomain2D(const MorphologicalTree& tree) noexcept { return *tree.gridDomain2D_; }

    /**
     * @brief Reads the owner of an established proper-part id.
     * @param tree Committed tree.
     * @param properPart Established proper-part id.
     * @return Owning node id.
     */
    [[nodiscard]] static NodeId properPartOwner(const MorphologicalTree& tree, NodeId properPart) noexcept {
        return tree.properPartOwner_[static_cast<std::size_t>(properPart)];
    }

    /** @brief Reads the parent of an established node. @param tree Committed tree. @param nodeId Established node id. @return Parent node id. */
    [[nodiscard]] static NodeId nodeParent(const MorphologicalTree& tree, NodeId nodeId) noexcept { return tree.nodeParent_[static_cast<std::size_t>(nodeId)]; }

    /** @brief Iterates children of an established node. @param tree Committed tree. @param nodeId Established node id. @return Child range. */
    [[nodiscard]] static MorphologicalTree::ChildrenRange children(const MorphologicalTree& tree, NodeId nodeId) {
        return MorphologicalTree::ChildrenRange(&tree, tree.firstChild_[static_cast<std::size_t>(nodeId)], tree.topologyVersion_);
    }

    /** @brief Iterates proper parts of an established node. @param tree Committed tree. @param nodeId Established node id. @return Proper-part range. */
    [[nodiscard]] static MorphologicalTree::ProperPartsRange properParts(const MorphologicalTree& tree, NodeId nodeId) {
        return MorphologicalTree::ProperPartsRange(&tree, tree.properHead_[static_cast<std::size_t>(nodeId)], tree.properPartVersion_);
    }

    /** @brief Counts children of an established node. @param tree Committed tree. @param nodeId Established node id. @return Child count. */
    [[nodiscard]] static int numChildren(const MorphologicalTree& tree, NodeId nodeId) noexcept {
        return tree.numChildrenByNode_[static_cast<std::size_t>(nodeId)];
    }

    /** @brief Counts proper parts of an established node. @param tree Committed tree. @param nodeId Established node id. @return Proper-part count. */
    [[nodiscard]] static int numProperParts(const MorphologicalTree& tree, NodeId nodeId) noexcept {
        return tree.numProperPartsByNode_[static_cast<std::size_t>(nodeId)];
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
        tree.ensurePrePostOrderCache();
        const auto isAncestorCached = [&](NodeId ancestor, NodeId node) {
            return tree.prePostOrderCache_.timePreOrder[static_cast<std::size_t>(ancestor)] <=
                       tree.prePostOrderCache_.timePreOrder[static_cast<std::size_t>(node)] &&
                   tree.prePostOrderCache_.timePostOrder[static_cast<std::size_t>(ancestor)] >=
                       tree.prePostOrderCache_.timePostOrder[static_cast<std::size_t>(node)];
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
        tree.ensurePrePostOrderCache();
        return tree.prePostOrderCache_.timePreOrder[static_cast<std::size_t>(ancestor)] <=
                   tree.prePostOrderCache_.timePreOrder[static_cast<std::size_t>(node)] &&
               tree.prePostOrderCache_.timePostOrder[static_cast<std::size_t>(ancestor)] >=
                   tree.prePostOrderCache_.timePostOrder[static_cast<std::size_t>(node)];
    }
};

} // namespace mmcfilters::detail
