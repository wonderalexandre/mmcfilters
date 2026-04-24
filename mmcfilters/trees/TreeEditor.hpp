#pragma once

#include "MorphologicalTree.hpp"

#include <stdexcept>

namespace mmcfilters {

/**
 * @brief Thin edit-session facade for multi-step topology updates.
 *
 * `TreeEditor` intentionally does not implement transactions or rollback. Its
 * role is to make explicit that the caller is performing a staged edit where
 * the tree may become temporarily disconnected, and to centralize the final
 * strong validation in `commit()`.
 */
class TreeEditor {
private:
    MorphologicalTree& tree_;

public:
    explicit TreeEditor(MorphologicalTree& tree) noexcept : tree_(tree) {}

    /**
     * @brief Creates a live detached node in the topological hierarchy.
     */
    NodeId createDetachedNode() {
        return tree_.createDetachedNode();
    }

    /**
     * @brief Detaches one non-root node from the connected rooted component.
     */
    void detach(NodeId nodeId) {
        if (!tree_.isAlive(nodeId)) {
            throw std::invalid_argument("TreeEditor::detach requires a live node.");
        }
        if (tree_.isRoot(nodeId)) {
            throw std::invalid_argument("TreeEditor::detach cannot detach the connected root.");
        }
        tree_.detachNode(nodeId);
    }

    /**
     * @brief Reparents one live non-root node under another live node.
     */
    void reparent(NodeId nodeId, NodeId newParentId) {
        if (!tree_.isAlive(nodeId) || !tree_.isAlive(newParentId)) {
            throw std::invalid_argument("TreeEditor::reparent requires live node ids.");
        }
        if (tree_.isRoot(nodeId)) {
            throw std::invalid_argument("TreeEditor::reparent cannot move the connected root.");
        }
        if (nodeId == newParentId) {
            throw std::invalid_argument("TreeEditor::reparent requires distinct node ids.");
        }
        tree_.moveNode(nodeId, newParentId);
    }

    /**
     * @brief Attaches one detached node back under the connected rooted tree.
     */
    void attach(NodeId parentId, NodeId detachedNodeId) {
        if (!tree_.isAlive(parentId) || !tree_.isAlive(detachedNodeId)) {
            throw std::invalid_argument("TreeEditor::attach requires live node ids.");
        }
        if (parentId == detachedNodeId) {
            throw std::invalid_argument("TreeEditor::attach requires distinct node ids.");
        }
        if (tree_.getNodeParent(detachedNodeId) != detachedNodeId) {
            throw std::invalid_argument("TreeEditor::attach expects a detached self-parented node.");
        }
        tree_.attachNode(parentId, detachedNodeId);
    }

    /**
     * @brief Transfers every direct child of `sourceId` under `parentId`.
     */
    void moveChildren(NodeId parentId, NodeId sourceId) {
        if (!tree_.isAlive(parentId) || !tree_.isAlive(sourceId)) {
            throw std::invalid_argument("TreeEditor::moveChildren requires live node ids.");
        }
        if (parentId == sourceId) {
            throw std::invalid_argument("TreeEditor::moveChildren requires distinct node ids.");
        }
        tree_.moveChildren(parentId, sourceId);
    }

    /**
     * @brief Transfers one direct proper part from `sourceNodeId` to `targetNodeId`.
     */
    void moveProperPart(NodeId targetNodeId, NodeId sourceNodeId, NodeId properPartId) {
        if (!tree_.isAlive(targetNodeId) || !tree_.isAlive(sourceNodeId)) {
            throw std::invalid_argument("TreeEditor::moveProperPart requires live node ids.");
        }
        if (targetNodeId == sourceNodeId) {
            throw std::invalid_argument("TreeEditor::moveProperPart requires distinct source and target nodes.");
        }
        if (!tree_.isProperPart(properPartId)) {
            throw std::invalid_argument("TreeEditor::moveProperPart requires a valid proper-part id.");
        }
        tree_.moveProperPart(targetNodeId, sourceNodeId, properPartId);
    }

    /**
     * @brief Transfers every direct proper part from `sourceNodeId` to `targetNodeId`.
     */
    void moveProperParts(NodeId targetNodeId, NodeId sourceNodeId) {
        if (!tree_.isAlive(targetNodeId) || !tree_.isAlive(sourceNodeId)) {
            throw std::invalid_argument("TreeEditor::moveProperParts requires live node ids.");
        }
        if (targetNodeId == sourceNodeId) {
            throw std::invalid_argument("TreeEditor::moveProperParts requires distinct source and target nodes.");
        }
        tree_.moveProperParts(targetNodeId, sourceNodeId);
    }

    /**
     * @brief Promotes `nodeId` to become the connected root.
     */
    void setRoot(NodeId nodeId) {
        if (!tree_.isAlive(nodeId)) {
            throw std::invalid_argument("TreeEditor::setRoot requires a live node.");
        }
        tree_.setRoot(nodeId);
    }

    /**
     * @brief Returns whether the staged edit still has detached alive nodes.
     */
    bool hasDetachedAliveNodes() const noexcept {
        return tree_.hasDetachedAliveNodes();
    }

    /**
     * @brief Finalizes the edit by validating that the tree is connected again.
     */
    void commit() const {
        tree_.validateConnectedRootedTree();
    }
};

} // namespace mmcfilters
