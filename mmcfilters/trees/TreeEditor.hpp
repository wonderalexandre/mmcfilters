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
 * strong validation in `validateAndCommit()` / `commit()`.
 *
 * Callers cannot construct `TreeEditor` directly. Use
 * `MorphologicalTree::edit()` so the edit boundary remains explicit at the
 * call site.
 */
class TreeEditor {
    friend class MorphologicalTree;

private:
    MorphologicalTree* tree_ = nullptr;
    bool active_ = false;

    /**
     * @brief Temporarily exposes the tree to committed-safe mutators.
     *
     * Public mutators such as `pruneNode()` and `mergeNodeIntoParent()` reject
     * open edit sessions. The editor uses this guard only around those complete
     * local operations, then restores the staged-edit flag even if they throw.
     */
    class EditSessionPause {
    private:
        MorphologicalTree& tree_;
        bool wasEditing_ = false;

    public:
        explicit EditSessionPause(MorphologicalTree& tree) noexcept
            : tree_(tree), wasEditing_(tree.editSessionOpen_) {
            tree_.editSessionOpen_ = false;
        }

        ~EditSessionPause() noexcept {
            tree_.editSessionOpen_ = wasEditing_;
        }
    };

    /**
     * @brief Opens an edit session on construction.
     */
    explicit TreeEditor(MorphologicalTree& tree) : tree_(&tree), active_(true) {
        tree_->beginEditSession();
    }

    /**
     * @brief Returns the edited tree or rejects use after commit/move-out.
     */
    MorphologicalTree& tree() const {
        if (!active_ || tree_ == nullptr) {
            throw std::logic_error("TreeEditor operation requires an active edit session.");
        }
        return *tree_;
    }

public:
    TreeEditor(const TreeEditor&) = delete;
    TreeEditor& operator=(const TreeEditor&) = delete;

    /**
     * @brief Transfers the open edit-session handle without closing it.
     */
    TreeEditor(TreeEditor&& other) noexcept
        : tree_(other.tree_), active_(other.active_) {
        other.tree_ = nullptr;
        other.active_ = false;
    }

    TreeEditor& operator=(TreeEditor&&) = delete;

    /**
     * @brief Leaves an unfinished edit session open.
     *
     * There is intentionally no rollback in the destructor. Destroying an active
     * editor without `commit()`, `validateAndCommit()`, or `commitUnchecked()`
     * keeps the tree in editing mode so committed-tree APIs continue to reject it.
     */
    ~TreeEditor() = default;

    /**
     * @brief Creates a live detached node in the topological hierarchy.
     *
     * The caller must attach or release the node before a checked commit can
     * succeed.
     */
    [[nodiscard("Discarding a detached node id makes the staged edit impossible to complete safely")]] NodeId createDetachedNode() {
        MorphologicalTree& t = tree();
        return t.createDetachedNode();
    }

    /**
     * @brief Detaches one non-root node from the connected rooted component.
     */
    void detach(NodeId nodeId) {
        MorphologicalTree& t = tree();
        if (!t.isAlive(nodeId)) {
            throw std::invalid_argument("TreeEditor::detach requires a live node.");
        }
        if (t.isRoot(nodeId)) {
            throw std::invalid_argument("TreeEditor::detach cannot detach the connected root.");
        }
        t.detachNode(nodeId);
    }

    /**
     * @brief Reparents one live non-root node under another live node.
     *
     * This is a staged structural edit: cycle freedom is enforced at commit
     * time rather than by an ancestry walk on every call.
     */
    void reparent(NodeId nodeId, NodeId newParentId) {
        MorphologicalTree& t = tree();
        if (!t.isAlive(nodeId) || !t.isAlive(newParentId)) {
            throw std::invalid_argument("TreeEditor::reparent requires live node ids.");
        }
        if (t.isRoot(nodeId)) {
            throw std::invalid_argument("TreeEditor::reparent cannot move the connected root.");
        }
        if (nodeId == newParentId) {
            throw std::invalid_argument("TreeEditor::reparent requires distinct node ids.");
        }
        t.moveNode(nodeId, newParentId);
    }

    /**
     * @brief Attaches one detached node back under the connected rooted tree.
     *
     * The node must be self-parented. Full connectivity and cycle validation
     * remains the responsibility of `commit()` / `validateAndCommit()`.
     */
    void attach(NodeId parentId, NodeId detachedNodeId) {
        MorphologicalTree& t = tree();
        if (!t.isAlive(parentId) || !t.isAlive(detachedNodeId)) {
            throw std::invalid_argument("TreeEditor::attach requires live node ids.");
        }
        if (parentId == detachedNodeId) {
            throw std::invalid_argument("TreeEditor::attach requires distinct node ids.");
        }
        if (t.getNodeParent(detachedNodeId) != detachedNodeId) {
            throw std::invalid_argument("TreeEditor::attach expects a detached self-parented node.");
        }
        t.attachNode(parentId, detachedNodeId);
    }

    /**
     * @brief Transfers every direct child of `sourceId` under `parentId`.
     *
     * This keeps the operation local and does not run an ancestry check; callers
     * that move children inside a subtree must repair or reject cycles before
     * committing.
     */
    void moveChildren(NodeId parentId, NodeId sourceId) {
        MorphologicalTree& t = tree();
        if (!t.isAlive(parentId) || !t.isAlive(sourceId)) {
            throw std::invalid_argument("TreeEditor::moveChildren requires live node ids.");
        }
        if (parentId == sourceId) {
            throw std::invalid_argument("TreeEditor::moveChildren requires distinct node ids.");
        }
        t.moveChildren(parentId, sourceId);
    }

    /**
     * @brief Transfers one direct proper part from `sourceNodeId` to `targetNodeId`.
     *
     * Moving one proper part is proportional to the linked-list update only; it
     * does not rebuild the full ownership index.
     */
    void moveProperPart(NodeId targetNodeId, NodeId sourceNodeId, NodeId properPartId) {
        MorphologicalTree& t = tree();
        if (!t.isAlive(targetNodeId) || !t.isAlive(sourceNodeId)) {
            throw std::invalid_argument("TreeEditor::moveProperPart requires live node ids.");
        }
        if (targetNodeId == sourceNodeId) {
            throw std::invalid_argument("TreeEditor::moveProperPart requires distinct source and target nodes.");
        }
        if (!t.isProperPart(properPartId)) {
            throw std::invalid_argument("TreeEditor::moveProperPart requires a valid proper-part id.");
        }
        t.moveProperPart(targetNodeId, sourceNodeId, properPartId);
    }

    /**
     * @brief Transfers every direct proper part from `sourceNodeId` to `targetNodeId`.
     *
     * The underlying splice is linear only in the number of moved direct proper
     * parts.
     */
    void moveProperParts(NodeId targetNodeId, NodeId sourceNodeId) {
        MorphologicalTree& t = tree();
        if (!t.isAlive(targetNodeId) || !t.isAlive(sourceNodeId)) {
            throw std::invalid_argument("TreeEditor::moveProperParts requires live node ids.");
        }
        if (targetNodeId == sourceNodeId) {
            throw std::invalid_argument("TreeEditor::moveProperParts requires distinct source and target nodes.");
        }
        t.moveProperParts(targetNodeId, sourceNodeId);
    }

    /**
     * @brief Detaches a direct child from its parent and optionally releases an empty detached slot.
     */
    void removeChild(NodeId parentNodeId, NodeId childId, bool releaseNodeFlag) {
        MorphologicalTree& t = tree();
        if (!t.isAlive(parentNodeId) || !t.isAlive(childId)) {
            throw std::invalid_argument("TreeEditor::removeChild requires live node ids.");
        }
        if (!t.hasChild(parentNodeId, childId)) {
            throw std::invalid_argument("TreeEditor::removeChild requires a direct parent-child relation.");
        }
        t.removeChild(parentNodeId, childId, releaseNodeFlag);
    }

    /**
     * @brief Releases an empty detached non-root node slot.
     */
    void releaseNode(NodeId nodeId) {
        MorphologicalTree& t = tree();
        if (!t.isAlive(nodeId)) {
            throw std::invalid_argument("TreeEditor::releaseNode requires a live node.");
        }
        if (t.isRoot(nodeId)) {
            throw std::invalid_argument("TreeEditor::releaseNode cannot release the connected root.");
        }
        if (t.getNodeParent(nodeId) != nodeId) {
            throw std::invalid_argument("TreeEditor::releaseNode expects a detached self-parented node.");
        }
        t.releaseNode(nodeId);
    }

    /**
     * @brief Promotes `nodeId` to become the connected root.
     *
     * The previous root becomes detached; callers must reconnect or explicitly
     * tolerate that state before a checked commit.
     */
    void setRoot(NodeId nodeId) {
        MorphologicalTree& t = tree();
        if (!t.isAlive(nodeId)) {
            throw std::invalid_argument("TreeEditor::setRoot requires a live node.");
        }
        t.setRoot(nodeId);
    }

    /**
     * @brief Applies the committed-safe subtree prune inside the staged edit.
     */
    void pruneNode(NodeId nodeId) {
        MorphologicalTree& t = tree();
        EditSessionPause pause(t);
        t.pruneNode(nodeId);
    }

    /**
     * @brief Applies the committed-safe parent merge inside the staged edit.
     */
    void mergeNodeIntoParent(NodeId nodeId) {
        MorphologicalTree& t = tree();
        EditSessionPause pause(t);
        t.mergeNodeIntoParent(nodeId);
    }

    /**
     * @brief Returns whether the staged edit still has detached alive nodes.
     */
    [[nodiscard]] bool hasDetachedAliveNodes() const noexcept {
        return active_ && tree_ != nullptr && tree_->hasDetachedAliveNodes();
    }

    /**
     * @brief Runs strong validation without closing the session.
     */
    [[nodiscard]] TreeValidationResult validate() const noexcept {
        if (!active_ || tree_ == nullptr) {
            return {false, "TreeEditor validation requires an active edit session."};
        }
        return tree_->validateConnectedRootedTreeResult();
    }

    /**
     * @brief Validates and closes the edit session on success.
     */
    [[nodiscard("Inspect the validation result or use commit() for exception-based failure handling")]] TreeValidationResult validateAndCommit() noexcept {
        TreeValidationResult result = validate();
        if (!result.ok) {
            return result;
        }
        tree_->endEditSession();
        active_ = false;
        return result;
    }

    /**
     * @brief Finalizes the edit by validating that the tree is connected again.
     *
     * This validation is linear in the current dense internal-node domain plus
     * the direct proper-part domain.
     */
    void commit() {
        TreeValidationResult result = validateAndCommit();
        if (!result.ok) {
            throw std::runtime_error(result.message);
        }
    }

    /**
     * @brief Finalizes a staged edit without running the linear validation.
     *
     * This is intended for internal algorithms that preserve the tree
     * invariants by construction and may commit inside a hot loop. Public and
     * test-facing edit sessions should keep using `commit()`.
     */
    void commitUnchecked() noexcept {
        if (active_ && tree_ != nullptr) {
#if defined(MMCFILTERS_ENABLE_ASSERTS)
            const TreeValidationResult result = validate();
            assert(result.ok && "TreeEditor::commitUnchecked requires a valid tree when MMCFILTERS_ENABLE_ASSERTS is enabled.");
#endif
            tree_->endEditSession();
            active_ = false;
        }
    }
};

inline TreeEditor MorphologicalTree::edit() {
    return TreeEditor(*this);
}

} // namespace mmcfilters
