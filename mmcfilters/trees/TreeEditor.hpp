#pragma once

#include "MorphologicalTree.hpp"

#include <cstdint>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>

namespace mmcfilters {

template <AltitudeValue T> class ValuedMorphologicalTreeEditor;

/**
 * @brief Thin edit-session facade for multi-step topology updates.
 *
 * Its role is to make explicit that the caller is performing a staged edit
 * where the tree may become temporarily disconnected, and to centralize the
 * final strong validation in `validateAndCommit()` / `commit()`. Public
 * editors keep an undo journal proportional to the mutation delta; the
 * internal by-construction editor remains journal-free.
 *
 * Callers cannot construct `TreeEditor` directly. Use
 * `MorphologicalTree::edit()` so the edit boundary remains explicit at the
 * call site.
 */
class TreeEditor {
    friend class MorphologicalTree;
    template <AltitudeValue T> friend class ValuedMorphologicalTreeEditor;

  public:
    /**
     * @brief Move-only evidence that the current edit revision satisfies the
     * generic topology invariants.
     */
    class IncrementalProof {
        friend class TreeEditor;

        /** @brief Editor. */
        const TreeEditor* editor_ = nullptr;
        /** @brief References the tree used by the component. */
        const MorphologicalTree* tree_ = nullptr;
        /** @brief Mutation version used to detect stale derived state. */
        std::size_t mutationVersion_ = 0;
        /** @brief Validation mode. */
        TreeEditValidationMode validationMode_ = TreeEditValidationMode::Complete;

        /**
         * @brief Constructs `IncrementalProof` from the supplied inputs.
         *
         * @param editor Tree editor associated with the proof.
         * @param validationMode Validation strategy applied to the edit.
         */
        IncrementalProof(const TreeEditor& editor, TreeEditValidationMode validationMode) noexcept
            : editor_(&editor), tree_(editor.tree_), mutationVersion_(editor.tree_ != nullptr ? editor.tree_->getMutationVersion() : 0),
              validationMode_(validationMode) {}

      public:
        /**
         * @brief Disables copy construction.
         */
        IncrementalProof(const IncrementalProof&) = delete;
        /**
         * @brief Disables copy assignment.
         */
        IncrementalProof& operator=(const IncrementalProof&) = delete;

        /**
         * @brief Transfers proof ownership and invalidates `other`.
         *
         * @param other Object to compare with or transfer from.
         */
        IncrementalProof(IncrementalProof&& other) noexcept
            : editor_(other.editor_), tree_(other.tree_), mutationVersion_(other.mutationVersion_), validationMode_(other.validationMode_) {
            other.editor_ = nullptr;
            other.tree_ = nullptr;
        }

        /**
         * @brief Disables move assignment.
         */
        IncrementalProof& operator=(IncrementalProof&&) = delete;

        /**
         * @brief Returns whether this proof was established by complete validation.
         *
         * @return Whether this proof was established by complete validation.
         */
        [[nodiscard]] bool usedCompleteValidation() const noexcept { return validationMode_ == TreeEditValidationMode::Complete; }
    };

  private:
    /**
     * @brief Flat open-addressing set of identifiers touched by the edit.
     *
     * Storage is proportional to the mutation delta. Repeated proper-part
     * moves between the same nodes remain O(1) insertions and do not grow the
     * ledger with duplicate ids.
     */
    template <class Id, Id InvalidId> class DeltaIdSet {
        /** @brief Identifiers buffer. */
        std::vector<Id> ids_;
        /** @brief Table buffer. */
        std::vector<Id> table_;

        /**
         * @brief Hashes an identifier for sparse-set storage.
         *
         * @param id Identifier.
         * @return Hash value for the identifier.
         */
        [[nodiscard]] static std::size_t hash(Id id) noexcept {
            return static_cast<std::size_t>(static_cast<std::uint32_t>(id) * std::uint32_t{2654435761u});
        }

        /**
         * @brief Rebuilds the identifier set with the requested capacity.
         *
         * @param capacity Requested storage capacity.
         */
        void rebuild(std::size_t capacity) {
            std::vector<Id> newTable(capacity, InvalidId);
            const std::size_t mask = capacity - 1;
            for (Id id : ids_) {
                std::size_t slot = hash(id) & mask;
                while (newTable[slot] != InvalidId) {
                    slot = (slot + 1) & mask;
                }
                newTable[slot] = id;
            }
            table_ = std::move(newTable);
        }

      public:
        /**
         * @brief Tests whether contains holds.
         *
         * @param id Identifier.
         * @return True when contains; otherwise false.
         */
        [[nodiscard]] bool contains(Id id) const noexcept {
            if (id == InvalidId || table_.empty()) {
                return false;
            }
            const std::size_t mask = table_.size() - 1;
            std::size_t slot = hash(id) & mask;
            while (table_[slot] != InvalidId) {
                if (table_[slot] == id) {
                    return true;
                }
                slot = (slot + 1) & mask;
            }
            return false;
        }

        /**
         * @brief Inserts an identifier into the set.
         *
         * @param id Identifier.
         * @return True when the documented condition holds; otherwise false.
         */
        [[nodiscard]] bool insert(Id id) {
            if (id == InvalidId) {
                return false;
            }
            if (table_.empty()) {
                ids_.reserve(16);
                rebuild(32);
            } else if ((ids_.size() + 1) * 10 > table_.size() * 7) {
                rebuild(table_.size() * 2);
            }

            const std::size_t mask = table_.size() - 1;
            std::size_t slot = hash(id) & mask;
            while (table_[slot] != InvalidId) {
                if (table_[slot] == id) {
                    return false;
                }
                slot = (slot + 1) & mask;
            }
            ids_.push_back(id);
            table_[slot] = id;
            return true;
        }

        [[nodiscard]] const std::vector<Id>&
        /**
         * @brief Returns the stored entries.
         *
         * @return The stored entries.
         */
        entries() const noexcept {
            return ids_;
        }
    };

    /** @brief Deduplicated node identifiers captured by one edit delta. */
    using DeltaNodeSet = DeltaIdSet<NodeId, InvalidNode>;

    /** @brief Deduplicated pixel identifiers captured by one edit delta. */
    using DeltaPixelSet = DeltaIdSet<PixelId, InvalidPixel>;

    /** @brief Captures one node slot so a recoverable tree edit can be rolled back. */
    struct NodeRollbackState {
        /** @brief Dense node identifier held by this record. */
        NodeId node = InvalidNode;
        /** @brief Dense node identifier of the parent. */
        NodeId parent = InvalidNode;
        /** @brief Dense node identifier of the first child. */
        NodeId firstChild = InvalidNode;
        /** @brief Dense node identifier of the next sibling. */
        NodeId nextSibling = InvalidNode;
        /** @brief Dense node identifier of the previous sibling. */
        NodeId prevSibling = InvalidNode;
        /** @brief Dense node identifier of the last child. */
        NodeId lastChild = InvalidNode;
        /** @brief Number of children. */
        int numChildren = 0;
        /** @brief Alive. */
        std::uint8_t alive = 0;
        /** @brief Pixel identifier of the proper head. */
        PixelId properHead = InvalidPixel;
        /** @brief Pixel identifier of the proper tail. */
        PixelId properTail = InvalidPixel;
        /** @brief Number of proper parts. */
        int properPartCardinality = 0;
    };

    /** @brief Captures one proper-part slot so a recoverable tree edit can be rolled back. */
    struct PixelRollbackState {
        /** @brief Identifies the pixel slot whose links are restored. */
        PixelId pixel = InvalidPixel;
        /** @brief Identifies the pixel's smallest node. */
        NodeId smallestNodeId = InvalidNode;
        /** @brief Identifies the next pixel in the same proper part. */
        PixelId next = InvalidPixel;
        /** @brief Identifies the previous pixel in the same proper part. */
        PixelId previous = InvalidPixel;
    };

    /** @brief Enumerates the supported free list mutation values. */
    enum class FreeListMutation { Popped, Pushed };

    /** @brief Records one free-list mutation for reverse replay during rollback. */
    struct FreeListRollbackState {
        /** @brief Mutation. */
        FreeListMutation mutation = FreeListMutation::Popped;
        /** @brief Dense node identifier held by this record. */
        NodeId node = InvalidNode;
    };

    /**
     * @brief Copy-on-first-write journal for recoverable public edits.
     *
     * Only directly affected node/proper-part slots are retained. Dense storage
     * growth is undone by restoring the original vector sizes, while free-list
     * pushes and pops are replayed in reverse.
     */
    struct RollbackJournal {
        /** @brief Original node slots. */
        std::size_t originalNodeSlots = 0;
        /** @brief Dense node identifier of the root. */
        NodeId root = InvalidNode;
        /** @brief Number of nodes. */
        int numNodes = 0;
        /** @brief Dense node identifier of the preserved external node identifier offset. */
        std::optional<NodeId> preservedExternalNodeIdOffset;
        /** @brief Node structure version used to detect stale derived state. */
        std::size_t nodeStructureVersion = 0;
        /** @brief Topology version used to detect stale derived state. */
        std::size_t topologyVersion = 0;
        /** @brief Proper part version used to detect stale derived state. */
        std::size_t properPartVersion = 0;
        /** @brief Mutation version used to detect stale derived state. */
        std::size_t mutationVersion = 0;
        /** @brief Captured nodes. */
        DeltaNodeSet capturedNodes;
        /** @brief Records pixel slots captured by the rollback journal. */
        DeltaPixelSet capturedPixels;
        /** @brief Nodes buffer. */
        std::vector<NodeRollbackState> nodes;
        /** @brief Captured pixel-link states buffer. */
        std::vector<PixelRollbackState> pixels;
        /** @brief Free list mutations buffer. */
        std::vector<FreeListRollbackState> freeListMutations;
    };

    /** @brief References the tree used by the component. */
    MorphologicalTree* tree_ = nullptr;
    /** @brief Indicates whether this editor owns the active edit session. */
    bool active_ = false;
    /** @brief Indicates whether the active edit can be rolled back. */
    bool recoverable_ = false;
    /** @brief Rollback journal. */
    std::unique_ptr<RollbackJournal> rollbackJournal_;
    /** @brief Touched nodes. */
    DeltaNodeSet touchedNodes_;
    /** @brief Detached node balance. */
    int detachedNodeBalance_ = 0;
    /** @brief Unsupported leaf balance. */
    int unsupportedLeafBalance_ = 0;
    /** @brief Indicates whether incremental validation remains supported. */
    bool incrementalValidationSupported_ = true;
    /** @brief Indicates whether construction directly established the required invariants. */
    bool invariantsEstablishedByConstruction_ = false;

    /**
     * @brief Temporarily exposes the tree to committed-safe mutators.
     *
     * Public mutators such as `pruneNode()` and `mergeNodeIntoParent()` reject
     * open edit sessions. The editor uses this guard only around those complete
     * local operations, then restores the staged-edit flag even if they throw.
     */
    class EditSessionPause {
      private:
        /** @brief References the tree used by the component. */
        MorphologicalTree& tree_;
        /** @brief Indicates whether an edit session was active before the pause. */
        bool wasEditing_ = false;

      public:
        /**
         * @brief Constructs `EditSessionPause` from the supplied inputs.
         *
         * @param tree Tree topology.
         */
        explicit EditSessionPause(MorphologicalTree& tree) noexcept : tree_(tree), wasEditing_(tree.editSessionOpen_) { tree_.editSessionOpen_ = false; }

        /**
         * @brief Destroys `EditSessionPause`.
         */
        ~EditSessionPause() noexcept { tree_.editSessionOpen_ = wasEditing_; }
    };

    /**
     * @brief Opens an edit session on construction.
     *
     * @param tree Tree topology.
     * @param invariantsEstablishedByConstruction Whether construction already established every edit invariant.
     */
    explicit TreeEditor(MorphologicalTree& tree, bool invariantsEstablishedByConstruction = false)
        : tree_(&tree), active_(true), recoverable_(!invariantsEstablishedByConstruction),
          invariantsEstablishedByConstruction_(invariantsEstablishedByConstruction) {
        tree_->beginEditSession();
    }

    /**
     * @brief Lazily creates the public rollback journal before first mutation.
     */
    void ensureRollbackJournal() {
        if (!recoverable_ || rollbackJournal_) {
            return;
        }
        auto journal = std::make_unique<RollbackJournal>();
        journal->originalNodeSlots = tree_->nodeParent_.size();
        journal->root = tree_->rootNodeId_;
        journal->numNodes = tree_->numNodes_;
        journal->preservedExternalNodeIdOffset = tree_->preservedExternalNodeIdOffset_;
        journal->nodeStructureVersion = tree_->nodeStructureVersion_;
        journal->topologyVersion = tree_->topologyVersion_;
        journal->properPartVersion = tree_->properPartVersion_;
        journal->mutationVersion = tree_->mutationVersion_;
        rollbackJournal_ = std::move(journal);
    }

    /**
     * @brief Returns the edited tree or rejects use after commit/move-out.
     *
     * @return The edited tree or rejects use after commit/move-out.
     */
    MorphologicalTree& tree() const {
        if (!active_ || tree_ == nullptr) {
            throw std::logic_error("TreeEditor operation requires an active edit session.");
        }
        return *tree_;
    }

    /**
     * @brief Closes an edit session after its caller established the invariants.
     *
     * This boundary is private so ordinary callers cannot publish an
     * unvalidated topology. The topology editor reaches it after validation;
     * the valuedTree editor reaches it only after its corresponding topology
     * and altitude proof has been established.
     *
     * @param validationMode Validation strategy used by the edit session.
     */
    void finishCommit(TreeEditValidationMode validationMode) noexcept {
        if (active_ && tree_ != nullptr) {
            tree_->recordEditCommit(validationMode);
            tree_->endEditSession();
            active_ = false;
            recoverable_ = false;
            rollbackJournal_.reset();
        }
    }

    /**
     * @brief Captures one original node slot before its first mutation.
     *
     * @param node Node identifier.
     */
    void captureNodeForRollback(NodeId node) {
        if (!recoverable_ || node < 0 || static_cast<std::size_t>(node) >= tree_->nodeParent_.size()) {
            return;
        }
        ensureRollbackJournal();
        if (static_cast<std::size_t>(node) >= rollbackJournal_->originalNodeSlots || rollbackJournal_->capturedNodes.contains(node)) {
            return;
        }
        rollbackJournal_->nodes.reserve(rollbackJournal_->nodes.size() + 1);
        if (!rollbackJournal_->capturedNodes.insert(node)) {
            return;
        }
        const std::size_t slot = static_cast<std::size_t>(node);
        rollbackJournal_->nodes.push_back({node, tree_->nodeParent_[slot], tree_->firstChild_[slot], tree_->nextSibling_[slot], tree_->prevSibling_[slot],
                                           tree_->lastChild_[slot], tree_->numChildrenByNode_[slot], tree_->alive_[slot], tree_->properHead_[slot],
                                           tree_->properTail_[slot], tree_->properPartCardinalityByNode_[slot]});
    }

    /**
     * @brief Captures one original proper-part slot before its first mutation.
     *
     * @param pixel Proper-part data.
     */
    void capturePixelForRollback(PixelId pixel) {
        if (!recoverable_ || !tree_->isPixel(pixel)) {
            return;
        }
        ensureRollbackJournal();
        if (rollbackJournal_->capturedPixels.contains(pixel)) {
            return;
        }
        rollbackJournal_->pixels.reserve(rollbackJournal_->pixels.size() + 1);
        if (!rollbackJournal_->capturedPixels.insert(pixel)) {
            return;
        }
        const std::size_t slot = static_cast<std::size_t>(pixel);
        rollbackJournal_->pixels.push_back({pixel, tree_->smallestNodeMap_[slot], tree_->nextProperPart_[slot], tree_->prevProperPart_[slot]});
    }

    /**
     * @brief Captures node link neighborhood.
     *
     * @param node Node identifier.
     */
    void captureNodeLinkNeighborhood(NodeId node) {
        if (node < 0 || static_cast<std::size_t>(node) >= tree_->nodeParent_.size()) {
            return;
        }
        captureNodeForRollback(node);
        captureNodeForRollback(tree_->nodeParent_[static_cast<std::size_t>(node)]);
        captureNodeForRollback(tree_->prevSibling_[static_cast<std::size_t>(node)]);
        captureNodeForRollback(tree_->nextSibling_[static_cast<std::size_t>(node)]);
    }

    /**
     * @brief Captures proper part link neighborhood.
     *
     * @param pixel Proper-part identifier.
     */
    void capturePixelLinkNeighborhood(PixelId pixel) {
        if (!tree_->isPixel(pixel)) {
            return;
        }
        capturePixelForRollback(pixel);
        capturePixelForRollback(tree_->prevProperPart_[static_cast<std::size_t>(pixel)]);
        capturePixelForRollback(tree_->nextProperPart_[static_cast<std::size_t>(pixel)]);
    }

    /**
     * @brief Prepares free list growth.
     *
     * @param count Number of entries involved in the operation.
     */
    void prepareFreeListGrowth(std::size_t count) {
        if (!recoverable_ || count == 0) {
            return;
        }
        ensureRollbackJournal();
        rollbackJournal_->freeListMutations.reserve(rollbackJournal_->freeListMutations.size() + count);
        tree_->freeNodeIds_.reserve(tree_->freeNodeIds_.size() + count);
    }

    /**
     * @brief Records free list growth.
     *
     * @param previousSize Count.
     */
    void recordFreeListGrowth(std::size_t previousSize) noexcept {
        if (!rollbackJournal_) {
            return;
        }
        for (std::size_t i = previousSize; i < tree_->freeNodeIds_.size(); ++i) {
            rollbackJournal_->freeListMutations.push_back({FreeListMutation::Pushed, tree_->freeNodeIds_[i]});
        }
    }

    /**
     * @brief Restores the copy-on-first-write journal without allocation.
     */
    void restoreRollbackJournal() noexcept {
        if (!active_ || tree_ == nullptr || !recoverable_) {
            return;
        }
        if (!rollbackJournal_) {
            tree_->endEditSession();
            recoverable_ = false;
            active_ = false;
            return;
        }

        for (auto it = rollbackJournal_->freeListMutations.rbegin(); it != rollbackJournal_->freeListMutations.rend(); ++it) {
            if (it->mutation == FreeListMutation::Pushed) {
                assert(!tree_->freeNodeIds_.empty());
                assert(tree_->freeNodeIds_.back() == it->node);
                tree_->freeNodeIds_.pop_back();
            } else {
                tree_->freeNodeIds_.push_back(it->node);
            }
        }

        const std::size_t originalSlots = rollbackJournal_->originalNodeSlots;
        tree_->nodeParent_.resize(originalSlots);
        tree_->firstChild_.resize(originalSlots);
        tree_->nextSibling_.resize(originalSlots);
        tree_->prevSibling_.resize(originalSlots);
        tree_->lastChild_.resize(originalSlots);
        tree_->numChildrenByNode_.resize(originalSlots);
        tree_->alive_.resize(originalSlots);
        tree_->properHead_.resize(originalSlots);
        tree_->properTail_.resize(originalSlots);
        tree_->properPartCardinalityByNode_.resize(originalSlots);

        for (const NodeRollbackState& state : rollbackJournal_->nodes) {
            const std::size_t slot = static_cast<std::size_t>(state.node);
            tree_->nodeParent_[slot] = state.parent;
            tree_->firstChild_[slot] = state.firstChild;
            tree_->nextSibling_[slot] = state.nextSibling;
            tree_->prevSibling_[slot] = state.prevSibling;
            tree_->lastChild_[slot] = state.lastChild;
            tree_->numChildrenByNode_[slot] = state.numChildren;
            tree_->alive_[slot] = state.alive;
            tree_->properHead_[slot] = state.properHead;
            tree_->properTail_[slot] = state.properTail;
            tree_->properPartCardinalityByNode_[slot] = state.properPartCardinality;
        }

        for (const PixelRollbackState& state : rollbackJournal_->pixels) {
            const std::size_t slot = static_cast<std::size_t>(state.pixel);
            tree_->smallestNodeMap_[slot] = state.smallestNodeId;
            tree_->nextProperPart_[slot] = state.next;
            tree_->prevProperPart_[slot] = state.previous;
        }

        tree_->rootNodeId_ = rollbackJournal_->root;
        tree_->numNodes_ = rollbackJournal_->numNodes;
        tree_->preservedExternalNodeIdOffset_ = rollbackJournal_->preservedExternalNodeIdOffset;
        tree_->nodeStructureVersion_ = rollbackJournal_->nodeStructureVersion;
        tree_->topologyVersion_ = rollbackJournal_->topologyVersion;
        tree_->properPartVersion_ = rollbackJournal_->properPartVersion;
        tree_->mutationVersion_ = rollbackJournal_->mutationVersion;
        tree_->invalidateDfsIntervalCache();
        tree_->invalidateLcaCache();
        tree_->endEditSession();
        rollbackJournal_.reset();
        recoverable_ = false;
        active_ = false;
    }

    /**
     * @brief Tests whether detached holds.
     *
     * @param tree Tree topology.
     * @param node Node identifier.
     * @return True when detached; otherwise false.
     */
    [[nodiscard]] static bool isDetached(const MorphologicalTree& tree, NodeId node) noexcept {
        return tree.isAlive(node) && node != tree.rootNodeId_ && tree.nodeParent_[static_cast<std::size_t>(node)] == node;
    }

    /**
     * @brief Records detached transition.
     *
     * @param wasDetached Records whether the node was detached during the edit.
     * @param isNowDetached Flag controlling is now detached.
     */
    void recordDetachedTransition(bool wasDetached, bool isNowDetached) noexcept {
        if (invariantsEstablishedByConstruction_) {
            return;
        }
        detachedNodeBalance_ += static_cast<int>(isNowDetached) - static_cast<int>(wasDetached);
    }

    /**
     * @brief Tests whether unsupported leaf holds.
     *
     * @param tree Tree topology.
     * @param node Node identifier.
     * @return True when unsupported leaf; otherwise false.
     */
    [[nodiscard]] static bool isUnsupportedLeaf(const MorphologicalTree& tree, NodeId node) noexcept {
        return tree.isAlive(node) && tree.numChildrenByNode_[static_cast<std::size_t>(node)] == 0 &&
               tree.properPartCardinalityByNode_[static_cast<std::size_t>(node)] == 0;
    }

    /**
     * @brief Records unsupported leaf transition.
     *
     * @param wasUnsupported Records whether the node lost all supporting descendants.
     * @param isNowUnsupported Flag controlling is now unsupported.
     */
    void recordUnsupportedLeafTransition(bool wasUnsupported, bool isNowUnsupported) noexcept {
        if (invariantsEstablishedByConstruction_) {
            return;
        }
        unsupportedLeafBalance_ += static_cast<int>(isNowUnsupported) - static_cast<int>(wasUnsupported);
    }

    /**
     * @brief Marks touch.
     *
     * @param node Node identifier.
     */
    void touch(NodeId node) {
        if (invariantsEstablishedByConstruction_) {
            return;
        }
        try {
            static_cast<void>(touchedNodes_.insert(node));
        } catch (...) {
            incrementalValidationSupported_ = false;
            throw;
        }
    }

  public:
    /**
     * @brief Disables copy construction.
     */
    TreeEditor(const TreeEditor&) = delete;
    /**
     * @brief Disables copy assignment.
     */
    TreeEditor& operator=(const TreeEditor&) = delete;

    /**
     * @brief Transfers the open edit-session handle without closing it.
     *
     * @param other Object to compare with or transfer from.
     */
    TreeEditor(TreeEditor&& other) noexcept
        : tree_(other.tree_), active_(other.active_), recoverable_(other.recoverable_), rollbackJournal_(std::move(other.rollbackJournal_)),
          touchedNodes_(std::move(other.touchedNodes_)), detachedNodeBalance_(other.detachedNodeBalance_),
          unsupportedLeafBalance_(other.unsupportedLeafBalance_), incrementalValidationSupported_(other.incrementalValidationSupported_),
          invariantsEstablishedByConstruction_(other.invariantsEstablishedByConstruction_) {
        other.tree_ = nullptr;
        other.active_ = false;
        other.recoverable_ = false;
        other.rollbackJournal_.reset();
        other.detachedNodeBalance_ = 0;
        other.unsupportedLeafBalance_ = 0;
        other.incrementalValidationSupported_ = false;
        other.invariantsEstablishedByConstruction_ = false;
    }

    /**
     * @brief Disables move assignment.
     */
    TreeEditor& operator=(TreeEditor&&) = delete;

    /**
     * @brief Rolls back an unfinished recoverable public edit.
     *
     * The internal by-construction editor is deliberately journal-free and
     * therefore still requires its algorithm to publish explicitly.
     */
    ~TreeEditor() { restoreRollbackJournal(); }

    /**
     * @brief Tests whether this editor owns a delta rollback journal.
     *
     * @return True if this editor owns a delta rollback journal; otherwise false.
     */
    [[nodiscard]] bool canRollback() const noexcept { return active_ && recoverable_; }

    /**
     * @brief Aborts a recoverable edit and restores its original state.
     */
    void rollback() {
        if (!active_ || tree_ == nullptr) {
            throw std::logic_error("TreeEditor::rollback requires an active edit session.");
        }
        if (!recoverable_) {
            throw std::logic_error("TreeEditor::rollback is unavailable for the internal journal-free editor.");
        }
        restoreRollbackJournal();
    }

    /**
     * @brief Creates a live detached node in the topological hierarchy.
     *
     * The caller must attach or release the node before a checked commit can
     * succeed.
     *
     * @return The created live detached node in the topological hierarchy.
     */
    [[nodiscard("Discarding a detached node id makes the staged edit impossible to complete safely")]] NodeId createDetachedNode() {
        MorphologicalTree& t = tree();
        if (invariantsEstablishedByConstruction_) {
            return t.createDetachedNode();
        }

        ensureRollbackJournal();
        const bool reusesFreeSlot = !t.freeNodeIds_.empty();
        const NodeId candidate = reusesFreeSlot ? t.freeNodeIds_.back() : static_cast<NodeId>(t.nodeParent_.size());
        captureNodeForRollback(candidate);
        if (reusesFreeSlot) {
            rollbackJournal_->freeListMutations.reserve(rollbackJournal_->freeListMutations.size() + 1);
            rollbackJournal_->freeListMutations.push_back({FreeListMutation::Popped, candidate});
        }
        const NodeId nodeId = t.createDetachedNode();
        recordDetachedTransition(false, isDetached(t, nodeId));
        recordUnsupportedLeafTransition(false, isUnsupportedLeaf(t, nodeId));
        touch(nodeId);
        return nodeId;
    }

    /**
     * @brief Detaches one non-root node from the connected rooted component.
     *
     * @param nodeId Dense internal node identifier.
     */
    void detach(NodeId nodeId) {
        MorphologicalTree& t = tree();
        if (!t.isAlive(nodeId)) {
            throw std::invalid_argument("TreeEditor::detach requires a live node.");
        }
        if (t.isRoot(nodeId)) {
            throw std::invalid_argument("TreeEditor::detach cannot detach the connected root.");
        }
        if (invariantsEstablishedByConstruction_) {
            t.detachNode(nodeId);
            return;
        }
        const NodeId oldParent = t.parent(nodeId);
        captureNodeLinkNeighborhood(nodeId);
        touch(nodeId);
        const bool wasDetached = isDetached(t, nodeId);
        const bool oldParentWasUnsupported = isUnsupportedLeaf(t, oldParent);
        t.detachNode(nodeId);
        recordDetachedTransition(wasDetached, isDetached(t, nodeId));
        recordUnsupportedLeafTransition(oldParentWasUnsupported, isUnsupportedLeaf(t, oldParent));
    }

    /**
     * @brief Reparents one live non-root node under another live node.
     *
     * This is a staged structural edit: cycle freedom is enforced at commit
     * time rather than by an ancestry walk on every call.
     *
     * @param nodeId Dense internal node identifier.
     * @param newParentId Parent-node value.
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
        if (invariantsEstablishedByConstruction_) {
            t.moveNode(nodeId, newParentId);
            return;
        }
        const NodeId oldParent = t.parent(nodeId);
        captureNodeLinkNeighborhood(nodeId);
        captureNodeForRollback(newParentId);
        captureNodeForRollback(t.lastChild_[static_cast<std::size_t>(newParentId)]);
        touch(nodeId);
        const bool wasDetached = isDetached(t, nodeId);
        const bool oldParentWasUnsupported = isUnsupportedLeaf(t, oldParent);
        const bool newParentWasUnsupported = isUnsupportedLeaf(t, newParentId);
        t.moveNode(nodeId, newParentId);
        recordDetachedTransition(wasDetached, isDetached(t, nodeId));
        recordUnsupportedLeafTransition(oldParentWasUnsupported, isUnsupportedLeaf(t, oldParent));
        if (oldParent != newParentId) {
            recordUnsupportedLeafTransition(newParentWasUnsupported, isUnsupportedLeaf(t, newParentId));
        }
    }

    /**
     * @brief Attaches one detached node back under the connected rooted tree.
     *
     * The node must be self-parented. Full connectivity and cycle validation
     * remains the responsibility of `commit()` / `validateAndCommit()`.
     *
     * @param parentId Identifier of the parent node.
     * @param detachedNodeId Node identifier.
     */
    void attach(NodeId parentId, NodeId detachedNodeId) {
        MorphologicalTree& t = tree();
        if (!t.isAlive(parentId) || !t.isAlive(detachedNodeId)) {
            throw std::invalid_argument("TreeEditor::attach requires live node ids.");
        }
        if (parentId == detachedNodeId) {
            throw std::invalid_argument("TreeEditor::attach requires distinct node ids.");
        }
        if (t.parent(detachedNodeId) != detachedNodeId) {
            throw std::invalid_argument("TreeEditor::attach expects a detached self-parented node.");
        }
        if (invariantsEstablishedByConstruction_) {
            t.attachNode(parentId, detachedNodeId);
            return;
        }
        captureNodeLinkNeighborhood(detachedNodeId);
        captureNodeForRollback(parentId);
        captureNodeForRollback(t.lastChild_[static_cast<std::size_t>(parentId)]);
        touch(detachedNodeId);
        const bool wasDetached = isDetached(t, detachedNodeId);
        const bool parentWasUnsupported = isUnsupportedLeaf(t, parentId);
        t.attachNode(parentId, detachedNodeId);
        recordDetachedTransition(wasDetached, isDetached(t, detachedNodeId));
        recordUnsupportedLeafTransition(parentWasUnsupported, isUnsupportedLeaf(t, parentId));
    }

    /**
     * @brief Transfers every direct child of `sourceId` under `parentId`.
     *
     * This keeps the operation local and does not run an ancestry check; callers
     * that move children inside a subtree must repair or reject cycles before
     * committing.
     *
     * @param parentId Identifier of the parent node.
     * @param sourceId Input.
     */
    void moveChildren(NodeId parentId, NodeId sourceId) {
        MorphologicalTree& t = tree();
        if (!t.isAlive(parentId) || !t.isAlive(sourceId)) {
            throw std::invalid_argument("TreeEditor::moveChildren requires live node ids.");
        }
        if (parentId == sourceId) {
            throw std::invalid_argument("TreeEditor::moveChildren requires distinct node ids.");
        }
        if (invariantsEstablishedByConstruction_) {
            t.moveChildren(parentId, sourceId);
            return;
        }
        captureNodeForRollback(parentId);
        captureNodeForRollback(sourceId);
        captureNodeForRollback(t.lastChild_[static_cast<std::size_t>(parentId)]);
        for (NodeId child = t.firstChild_[static_cast<std::size_t>(sourceId)]; child != InvalidNode; child = t.nextSibling_[static_cast<std::size_t>(child)]) {
            captureNodeForRollback(child);
        }
        touch(parentId);
        const bool parentWasUnsupported = isUnsupportedLeaf(t, parentId);
        const bool sourceWasUnsupported = isUnsupportedLeaf(t, sourceId);
        t.moveChildren(parentId, sourceId);
        recordUnsupportedLeafTransition(parentWasUnsupported, isUnsupportedLeaf(t, parentId));
        recordUnsupportedLeafTransition(sourceWasUnsupported, isUnsupportedLeaf(t, sourceId));
    }

    /**
     * @brief Transfers one direct proper part from `sourceNodeId` to `targetNodeId`.
     *
     * Moving one proper part is proportional to the linked-list update only; it
     * does not rebuild the full smallest-node map.
     *
     * @param targetNodeId Node identifier.
     * @param sourceNodeId Node identifier.
     * @param pixel Proper-part identifier.
     */
    void movePixelToProperPart(NodeId targetNodeId, NodeId sourceNodeId, PixelId pixel) {
        MorphologicalTree& t = tree();
        if (!t.isAlive(targetNodeId) || !t.isAlive(sourceNodeId)) {
            throw std::invalid_argument("TreeEditor::movePixelToProperPart requires live node ids.");
        }
        if (targetNodeId == sourceNodeId) {
            throw std::invalid_argument("TreeEditor::movePixelToProperPart requires distinct source and target nodes.");
        }
        if (!t.isPixel(pixel)) {
            throw std::invalid_argument("TreeEditor::movePixelToProperPart requires a valid proper-part id.");
        }
        if (invariantsEstablishedByConstruction_) {
            t.movePixelToProperPart(targetNodeId, sourceNodeId, pixel);
            return;
        }
        captureNodeForRollback(targetNodeId);
        captureNodeForRollback(sourceNodeId);
        capturePixelLinkNeighborhood(pixel);
        capturePixelForRollback(t.properTail_[static_cast<std::size_t>(targetNodeId)]);
        const bool targetWasUnsupported = isUnsupportedLeaf(t, targetNodeId);
        const bool sourceWasUnsupported = isUnsupportedLeaf(t, sourceNodeId);
        t.movePixelToProperPart(targetNodeId, sourceNodeId, pixel);
        recordUnsupportedLeafTransition(targetWasUnsupported, isUnsupportedLeaf(t, targetNodeId));
        recordUnsupportedLeafTransition(sourceWasUnsupported, isUnsupportedLeaf(t, sourceNodeId));
    }

    /**
     * @brief Transfers every direct proper part from `sourceNodeId` to `targetNodeId`.
     *
     * The underlying splice is linear only in the number of moved direct proper
     * parts.
     *
     * @param targetNodeId Node identifier.
     * @param sourceNodeId Node identifier.
     */
    void mergeProperParts(NodeId targetNodeId, NodeId sourceNodeId) {
        MorphologicalTree& t = tree();
        if (!t.isAlive(targetNodeId) || !t.isAlive(sourceNodeId)) {
            throw std::invalid_argument("TreeEditor::mergeProperParts requires live node ids.");
        }
        if (targetNodeId == sourceNodeId) {
            throw std::invalid_argument("TreeEditor::mergeProperParts requires distinct source and target nodes.");
        }
        if (invariantsEstablishedByConstruction_) {
            t.mergeProperParts(targetNodeId, sourceNodeId);
            return;
        }
        captureNodeForRollback(targetNodeId);
        captureNodeForRollback(sourceNodeId);
        capturePixelForRollback(t.properTail_[static_cast<std::size_t>(targetNodeId)]);
        for (PixelId pixel = t.properHead_[static_cast<std::size_t>(sourceNodeId)]; pixel != InvalidPixel;
             pixel = t.nextProperPart_[static_cast<std::size_t>(pixel)]) {
            capturePixelForRollback(pixel);
        }
        const bool targetWasUnsupported = isUnsupportedLeaf(t, targetNodeId);
        const bool sourceWasUnsupported = isUnsupportedLeaf(t, sourceNodeId);
        t.mergeProperParts(targetNodeId, sourceNodeId);
        recordUnsupportedLeafTransition(targetWasUnsupported, isUnsupportedLeaf(t, targetNodeId));
        recordUnsupportedLeafTransition(sourceWasUnsupported, isUnsupportedLeaf(t, sourceNodeId));
    }

    /**
     * @brief Detaches a direct child from its parent and optionally releases an empty detached slot.
     *
     * @param parentNodeId Identifier of the parent node.
     * @param childId Identifier of the child node.
     * @param releaseNodeFlag Flag controlling release node flag.
     */
    void removeChild(NodeId parentNodeId, NodeId childId, bool releaseNodeFlag) {
        MorphologicalTree& t = tree();
        if (!t.isAlive(parentNodeId) || !t.isAlive(childId)) {
            throw std::invalid_argument("TreeEditor::removeChild requires live node ids.");
        }
        if (!t.hasChild(parentNodeId, childId)) {
            throw std::invalid_argument("TreeEditor::removeChild requires a direct parent-child relation.");
        }
        if (invariantsEstablishedByConstruction_) {
            t.removeChild(parentNodeId, childId, releaseNodeFlag);
            return;
        }
        captureNodeLinkNeighborhood(childId);
        captureNodeForRollback(parentNodeId);
        const bool willRelease =
            releaseNodeFlag && t.numChildrenByNode_[static_cast<std::size_t>(childId)] == 0 && t.properPartCardinalityByNode_[static_cast<std::size_t>(childId)] == 0;
        if (willRelease) {
            prepareFreeListGrowth(1);
        }
        touch(childId);
        const bool wasDetached = isDetached(t, childId);
        const bool parentWasUnsupported = isUnsupportedLeaf(t, parentNodeId);
        const bool childWasUnsupported = isUnsupportedLeaf(t, childId);
        const std::size_t freeListSize = t.freeNodeIds_.size();
        t.removeChild(parentNodeId, childId, releaseNodeFlag);
        recordFreeListGrowth(freeListSize);
        recordDetachedTransition(wasDetached, isDetached(t, childId));
        recordUnsupportedLeafTransition(parentWasUnsupported, isUnsupportedLeaf(t, parentNodeId));
        recordUnsupportedLeafTransition(childWasUnsupported, isUnsupportedLeaf(t, childId));
    }

    /**
     * @brief Releases an empty detached non-root node slot.
     *
     * @param nodeId Dense internal node identifier.
     */
    void releaseNode(NodeId nodeId) {
        MorphologicalTree& t = tree();
        if (!t.isAlive(nodeId)) {
            throw std::invalid_argument("TreeEditor::releaseNode requires a live node.");
        }
        if (t.isRoot(nodeId)) {
            throw std::invalid_argument("TreeEditor::releaseNode cannot release the connected root.");
        }
        if (t.parent(nodeId) != nodeId) {
            throw std::invalid_argument("TreeEditor::releaseNode expects a detached self-parented node.");
        }
        if (invariantsEstablishedByConstruction_) {
            t.releaseNode(nodeId);
            return;
        }
        captureNodeForRollback(nodeId);
        const bool willRelease = t.numChildrenByNode_[static_cast<std::size_t>(nodeId)] == 0 && t.properPartCardinalityByNode_[static_cast<std::size_t>(nodeId)] == 0;
        if (willRelease) {
            prepareFreeListGrowth(1);
        }
        const bool wasDetached = isDetached(t, nodeId);
        const bool wasUnsupported = isUnsupportedLeaf(t, nodeId);
        const std::size_t freeListSize = t.freeNodeIds_.size();
        t.releaseNode(nodeId);
        recordFreeListGrowth(freeListSize);
        recordDetachedTransition(wasDetached, isDetached(t, nodeId));
        recordUnsupportedLeafTransition(wasUnsupported, isUnsupportedLeaf(t, nodeId));
    }

    /**
     * @brief Promotes `nodeId` to become the connected root.
     *
     * The previous root becomes detached; callers must reconnect or explicitly
     * tolerate that state before a checked commit.
     *
     * @param nodeId Dense internal node identifier.
     */
    void setRoot(NodeId nodeId) {
        MorphologicalTree& t = tree();
        if (!t.isAlive(nodeId)) {
            throw std::invalid_argument("TreeEditor::setRoot requires a live node.");
        }
        if (invariantsEstablishedByConstruction_) {
            t.setRoot(nodeId);
            return;
        }
        const NodeId oldRoot = t.root();
        const NodeId oldParent = t.parent(nodeId);
        captureNodeForRollback(oldRoot);
        captureNodeLinkNeighborhood(nodeId);
        captureNodeForRollback(oldParent);
        touch(oldRoot);
        touch(oldParent);
        touch(nodeId);
        const bool oldRootWasDetached = isDetached(t, oldRoot);
        const bool nodeWasDetached = isDetached(t, nodeId);
        t.setRoot(nodeId);
        recordDetachedTransition(oldRootWasDetached, isDetached(t, oldRoot));
        if (nodeId != oldRoot) {
            recordDetachedTransition(nodeWasDetached, isDetached(t, nodeId));
        }
    }

    /**
     * @brief Applies the committed-safe subtree prune inside the staged edit.
     *
     * @param nodeId Dense internal node identifier.
     */
    void pruneNode(NodeId nodeId) {
        MorphologicalTree& t = tree();
        if (invariantsEstablishedByConstruction_) {
            EditSessionPause pause(t);
            t.pruneNode(nodeId);
            return;
        }
        incrementalValidationSupported_ = false;

        if (t.isAlive(nodeId) && !t.isRoot(nodeId)) {
            const NodeId parent = t.nodeParent_[static_cast<std::size_t>(nodeId)];
            if (parent != InvalidNode && parent != nodeId) {
                std::vector<NodeId> subtree;
                subtree.push_back(nodeId);
                for (std::size_t i = 0; i < subtree.size(); ++i) {
                    const NodeId current = subtree[i];
                    for (NodeId child = t.firstChild_[static_cast<std::size_t>(current)]; child != InvalidNode;
                         child = t.nextSibling_[static_cast<std::size_t>(child)]) {
                        subtree.push_back(child);
                    }
                }

                captureNodeForRollback(parent);
                captureNodeLinkNeighborhood(nodeId);
                capturePixelForRollback(t.properTail_[static_cast<std::size_t>(parent)]);
                for (NodeId current : subtree) {
                    captureNodeForRollback(current);
                    for (PixelId pixel = t.properHead_[static_cast<std::size_t>(current)]; pixel != InvalidPixel;
                         pixel = t.nextProperPart_[static_cast<std::size_t>(pixel)]) {
                        capturePixelForRollback(pixel);
                    }
                }
                prepareFreeListGrowth(subtree.size());
            }
        }

        const std::size_t freeListSize = t.freeNodeIds_.size();
        EditSessionPause pause(t);
        t.pruneNode(nodeId);
        recordFreeListGrowth(freeListSize);
    }

    /**
     * @brief Applies the committed-safe parent merge inside the staged edit.
     *
     * @param nodeId Dense internal node identifier.
     */
    void mergeNodeIntoParent(NodeId nodeId) {
        MorphologicalTree& t = tree();
        if (invariantsEstablishedByConstruction_) {
            EditSessionPause pause(t);
            t.mergeNodeIntoParent(nodeId);
            return;
        }
        incrementalValidationSupported_ = false;

        if (t.isAlive(nodeId) && !t.isRoot(nodeId)) {
            const NodeId parent = t.nodeParent_[static_cast<std::size_t>(nodeId)];
            if (parent != InvalidNode && parent != nodeId) {
                captureNodeForRollback(parent);
                captureNodeLinkNeighborhood(nodeId);
                capturePixelForRollback(t.properTail_[static_cast<std::size_t>(parent)]);
                for (PixelId pixel = t.properHead_[static_cast<std::size_t>(nodeId)]; pixel != InvalidPixel;
                     pixel = t.nextProperPart_[static_cast<std::size_t>(pixel)]) {
                    capturePixelForRollback(pixel);
                }
                for (NodeId child = t.firstChild_[static_cast<std::size_t>(nodeId)]; child != InvalidNode;
                     child = t.nextSibling_[static_cast<std::size_t>(child)]) {
                    captureNodeForRollback(child);
                }
                prepareFreeListGrowth(1);
            }
        }

        const std::size_t freeListSize = t.freeNodeIds_.size();
        EditSessionPause pause(t);
        t.mergeNodeIntoParent(nodeId);
        recordFreeListGrowth(freeListSize);
    }

    /**
     * @brief Returns whether the staged edit still has detached alive nodes.
     *
     * @return Whether the staged edit still has detached alive nodes.
     */
    [[nodiscard]] bool hasDetachedAliveNodes() const noexcept { return active_ && tree_ != nullptr && tree_->hasDetachedAliveNodes(); }

  private:
    /**
     * @brief Validates only nodes and parent paths affected by supported edit
     * primitives.
     *
     * @param changedParentCyclesExcluded Parent-node value.
     * @return Validation result for only nodes and parent paths affected by supported edit primitives.
     */
    [[nodiscard]] TreeValidationResult validateIncrementalTopology(bool changedParentCyclesExcluded) const noexcept {
        try {
            if (!active_ || tree_ == nullptr) {
                return {false, "Incremental topology validation requires an active edit session."};
            }
            if (detachedNodeBalance_ != 0) {
                return {false, "Incremental topology validation found detached alive nodes."};
            }
            if (unsupportedLeafBalance_ != 0) {
                return {false, "Incremental topology validation found a live node whose subtree support is empty."};
            }

            const NodeId root = tree_->root();
            if (!tree_->isAlive(root) || tree_->parent(root) != root) {
                return {false, "Incremental topology validation requires a live self-parented root."};
            }

            const int numNodes = tree_->numNodes();
            for (NodeId node : touchedNodes_.entries()) {
                if (!tree_->isAlive(node)) {
                    continue;
                }

                if (!changedParentCyclesExcluded) {
                    NodeId cursor = node;
                    int pathLength = 0;
                    while (cursor != root) {
                        if (!tree_->isAlive(cursor)) {
                            return {false, "Incremental topology validation found a parent path outside the alive node domain."};
                        }
                        const NodeId parent = tree_->parent(cursor);
                        if (parent == InvalidNode || parent == cursor || !tree_->isAlive(parent)) {
                            return {false, "Incremental topology validation found a detached or invalid parent path."};
                        }
                        cursor = parent;
                        if (++pathLength > numNodes) {
                            return {false, "Incremental topology validation found a parent cycle."};
                        }
                    }
                }
            }
            return {true, ""};
        } catch (const std::exception& ex) {
            return {false, ex.what()};
        } catch (...) {
            return {false, "Incremental topology validation failed with an unknown error."};
        }
    }

    /**
     * @brief Builds the validation proof for the current incremental edit state.
     *
     * @param changedParentCyclesExcluded Parent-related value.
     * @return Proof describing the invariants established by the edit.
     */
    [[nodiscard]] IncrementalProof proveIncrementalImpl(bool changedParentCyclesExcluded) {
        TreeValidationResult result;
        TreeEditValidationMode validationMode = TreeEditValidationMode::Incremental;
        if (invariantsEstablishedByConstruction_ && incrementalValidationSupported_) {
            result = {true, ""};
        } else if (incrementalValidationSupported_) {
            result = validateIncrementalTopology(changedParentCyclesExcluded);
        } else {
            result = validate();
            validationMode = TreeEditValidationMode::Complete;
        }
        if (!result.ok) {
            throw std::runtime_error(result.message);
        }
#ifndef NDEBUG
        if (validationMode == TreeEditValidationMode::Incremental) {
            const TreeValidationResult oracle = validate();
            if (!oracle.ok) {
                throw std::runtime_error(std::string("Incremental topology proof disagrees with the complete validation oracle: ") + oracle.message);
            }
        }
#endif
        return IncrementalProof(*this, validationMode);
    }

    /**
     * @brief Uses strict monotone changed arcs as generic acyclicity evidence.
     *
     * Only `ValuedMorphologicalTreeEditor` can issue this evidence after checking the
     * declared altitude order around every touched node.
     *
     * @return Proof produced by uses strict monotone changed arcs as generic acyclicity evidence.
     */
    [[nodiscard]] IncrementalProof proveIncrementalWithStrictAltitudeAcyclicity() { return proveIncrementalImpl(true); }

  public:
    /**
     * @brief Produces move-only evidence for the current edit revision.
     *
     * Supported primitives are checked only on the mutation delta and changed
     * parent paths. An edit containing a primitive without incremental support
     * falls back to the complete validator.
     *
     * @return The produced move-only evidence for the current edit revision.
     */
    [[nodiscard]] IncrementalProof proveIncremental() { return proveIncrementalImpl(false); }

    /**
     * @brief Commits the exact edit revision represented by `proof`.
     *
     * @param proof Validation proof consumed by the operation.
     */
    void commit(IncrementalProof&& proof) {
        if (!active_ || tree_ == nullptr || proof.editor_ != this || proof.tree_ != tree_ || proof.mutationVersion_ != tree_->getMutationVersion()) {
            throw std::logic_error("Incremental topology proof is stale or belongs to another edit session.");
        }

        const TreeEditValidationMode validationMode = proof.validationMode_;
        proof.editor_ = nullptr;
        proof.tree_ = nullptr;
        finishCommit(validationMode);
    }

    /**
     * @brief Runs strong validation without closing the session.
     *
     * @return Result produced by running strong validation without closing the session.
     */
    [[nodiscard]] TreeValidationResult validate() const noexcept {
        if (!active_ || tree_ == nullptr) {
            return {false, "TreeEditor validation requires an active edit session."};
        }
        return tree_->validateConnectedRootedTreeResult();
    }

    /**
     * @brief Validates and closes the edit session on success.
     *
     * @return Validation result for and closes the edit session on success.
     */
    [[nodiscard("Inspect the validation result or use commit() for exception-based failure handling")]] TreeValidationResult validateAndCommit() noexcept {
        TreeValidationResult result = validate();
        if (!result.ok) {
            return result;
        }
        finishCommit(TreeEditValidationMode::Complete);
        return result;
    }

    /**
     * @brief Finalizes the edit by validating that the tree is connected again.
     *
     * This validation is linear in the current dense internal-node domain plus
     * the direct pixel domain.
     */
    void commit() {
        TreeValidationResult result = validateAndCommit();
        if (!result.ok) {
            throw std::runtime_error(result.message);
        }
    }
};

inline TreeEditor MorphologicalTree::edit() { return TreeEditor(*this); }

} // namespace mmcfilters
