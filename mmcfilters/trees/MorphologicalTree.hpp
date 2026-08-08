#pragma once

#include "../utils/RegularGridAdjacency2D.hpp"
#include "../utils/Assert.hpp"
#include "../utils/Altitude.hpp"
#include "../utils/Common.hpp"
#include "../dataStructure/FastQueue.hpp"
#include "HierarchySemantics.hpp"
#include "ProperPartDomain.hpp"
#include "detail/MorphologicalTreeConstructionTag.hpp"
#include "detail/NativeHierarchyValidationDetail.hpp"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <list>
#include <memory>
#include <optional>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace mmcfilters {

/**
 * @brief Selects the node-id domain used by attribute buffers exposed to callers.
 */
enum class NodeIdSpace { MORPHOLOGICAL_TREE, HIGRA };

/**
 * @brief Non-throwing validation result returned by edit-session checks.
 *
 * `ok` carries the success flag and `message` carries a diagnostic when
 * validation fails. The explicit boolean conversion mirrors `ok` so callers can
 * write compact checks while still preserving the detailed message.
 */
struct [[nodiscard]] TreeValidationResult {
    /// True when the checked tree topology satisfies the validation contract.
    bool ok = false;

    /// Human-readable diagnostic, especially useful when `ok` is false.
    std::string message;

    /**
     * @brief Allows validation results to be tested directly in boolean contexts.
     */
    explicit operator bool() const noexcept { return ok; }
};

/**
 * @brief Counts committed edit sessions by validation strategy.
 *
 * The counters make performance-sensitive callers testable without embedding
 * timing thresholds in correctness tests.
 */
struct TreeEditValidationStatistics {
    /// Number of commits checked with complete validation.
    std::size_t completeValidationCommits = 0;
    /// Number of commits checked with incremental proof validation.
    std::size_t incrementalValidationCommits = 0;
};

enum class TreeEditValidationMode { Complete, Incremental };

// Forward declaration for the edit-session wrapper.
class TreeEditor;

/**
 * @brief Mutable morphological tree built directly on proper parts and dense node ids.
 *
 * `MorphologicalTree` is the central mutable hierarchy of this project. It can
 * represent a rooted morphological tree of partial partitions and exposes a
 * dense `NodeId` domain together with direct proper-part ownership. The class keeps
 * explicit parent/child links, linked lists of direct proper parts, optional
 * regular 2D proper-part-domain metadata, and a small set of
 * structural caches used by the
 * public traversal and ancestry queries.
 *
 * Data model:
 *
 * - proper parts are indexed by `NodeId` in the range `[0, getNumTotalProperParts())`;
 * - internal nodes are indexed by `NodeId` in the range `[0, getNumInternalNodeSlots())`;
 * - each live node owns zero or more direct proper parts and may have direct children;
 * - the full support of a node is the union of the direct proper parts in its subtree;
 * - every node in a checked committed hierarchy has non-empty full support;
 * - the root and detached nodes point to themselves as parent.
 *
 * Main responsibilities:
 *
 * - expose mutable topology operations used by filters and adjusters;
 * - provide structural traversals over live nodes, children, proper parts, and subtrees;
 * - serve as the topology/ownership core used by weighted wrappers such as
 *   `WeightedMorphologicalTree<std::uint8_t>`.
 *
 * The class intentionally keeps only the canonical structural state and a small
 * number of derived caches. Higher-level attribute computation is delegated to
 * the incremental attribute computers in `mmcfilters/attributes`.
 */
class MorphologicalTree {
  private:
    friend class TreeEditor;
    friend class MorphologicalTreeFactory;

    class LCAEulerRMQ; // Forward declaration for the LCA cache implementation.

    // ========================= Private attributes ========================= //
    /** @brief Stores the root node identifier. */
    NodeId rootNodeId_ = InvalidNode;
    /** @brief Stores the semantics. */
    HierarchySemantics semantics_;
    /** @brief Stores the grid domain2 d. */
    std::optional<GridDomain2D> gridDomain2D_;
    /** @brief Stores the num nodes. */
    int numNodes_ = 0;
    /** @brief Stores the preserved external node identifier offset. */
    std::optional<NodeId> preservedExternalNodeIdOffset_;
    /** @brief Indicates whether a guarded edit session is open. */
    bool editSessionOpen_ = false;
    /** @brief Stores the edit validation statistics. */
    TreeEditValidationStatistics editValidationStatistics_;

    // Proper-part ownership, indexed by proper-part global id [0, getNumTotalProperParts()).
    /** @brief Stores the proper part owner. */
    std::vector<NodeId> properPartOwner_;

    // Parent links, indexed by local node-slot id [0, getNumInternalNodeSlots()).
    /** @brief Stores the node parent. */
    std::vector<NodeId> nodeParent_;

    // Internal hierarchy linked structure, indexed by local node-slot id [0, getNumInternalNodeSlots()).
    /** @brief Stores the first child. */
    std::vector<NodeId> firstChild_;
    /** @brief Stores the next sibling. */
    std::vector<NodeId> nextSibling_;
    /** @brief Stores the prev sibling. */
    std::vector<NodeId> prevSibling_;
    /** @brief Stores the last child. */
    std::vector<NodeId> lastChild_;
    /** @brief Stores the num children by node. */
    std::vector<int> numChildrenByNode_;

    // Free slot management and alive-node iteration, indexed by local node-slot id [0, getNumInternalNodeSlots()).
    /** @brief Stores the alive. */
    std::vector<uint8_t> alive_;
    /** @brief Stores the free node identifiers. */
    std::vector<NodeId> freeNodeIds_;

    // Proper-parts linked lists, indexed by local node-slot id [0, getNumTotalProperParts()).
    /** @brief Stores the proper head. */
    std::vector<NodeId> properHead_;
    /** @brief Stores the proper tail. */
    std::vector<NodeId> properTail_;
    /** @brief Stores the num proper parts by node. */
    std::vector<int> numProperPartsByNode_;
    /** @brief Stores the next proper part. */
    std::vector<NodeId> nextProperPart_;
    /** @brief Stores the prev proper part. */
    std::vector<NodeId> prevProperPart_;

    // Structural caches for traversal-based queries, indexed by local node-slot id [0, getNumInternalNodeSlots()).
    /** @brief Caches preorder and postorder indexes for the current topology version. */
    struct PrePostOrderCache {
        /** @brief Stores the time pre order. */
        std::vector<int> timePreOrder;
        /** @brief Stores the time post order. */
        std::vector<int> timePostOrder;
        /** @brief Indicates whether the cached preorder/postorder data is valid. */
        bool valid = false;

        /**
         * @brief Marks the cached traversal timestamps as stale.
         */
        void invalidate() noexcept { valid = false; }
    };
    /** @brief Stores the pre post order cache. */
    mutable PrePostOrderCache prePostOrderCache_;
    /** @brief Stores the lowest-common-ancestor cache. */
    mutable std::unique_ptr<LCAEulerRMQ> lcaCache_;

    // Version counters for iterator invalidation.
    /** @brief Stores the node structure version. */
    std::size_t nodeStructureVersion_ = 0;
    /** @brief Stores the topology version. */
    std::size_t topologyVersion_ = 0;
    /** @brief Stores the proper part version. */
    std::size_t properPartVersion_ = 0;
    /** @brief Stores the mutation version. */
    std::size_t mutationVersion_ = 0;

    /** @brief Enumerates the supported child splice policy values. */
    enum class ChildSplicePolicy { AppendToTargetTail, ReplaceSourceSlotWhenDirectChild };

    // ========================= Private methods ========================= //
    /**
     * @brief Returns a reusable slot or appends a fresh one to the dense storage.
     *
     * @return A reusable slot or appends a fresh one to the dense storage.
     */
    inline NodeId allocateSlot() {
        if (!freeNodeIds_.empty()) {
            const NodeId slotId = freeNodeIds_.back();
            freeNodeIds_.pop_back();

            nodeParent_[slotId] = InvalidNode;
            firstChild_[slotId] = InvalidNode;
            nextSibling_[slotId] = InvalidNode;
            prevSibling_[slotId] = InvalidNode;
            lastChild_[slotId] = InvalidNode;
            numChildrenByNode_[slotId] = 0;
            alive_[slotId] = 1;
            properHead_[slotId] = InvalidNode;
            properTail_[slotId] = InvalidNode;
            numProperPartsByNode_[slotId] = 0;
            return slotId;
        }

        const NodeId slotId = static_cast<NodeId>(nodeParent_.size());
        nodeParent_.push_back(InvalidNode);
        firstChild_.push_back(InvalidNode);
        nextSibling_.push_back(InvalidNode);
        prevSibling_.push_back(InvalidNode);
        lastChild_.push_back(InvalidNode);
        numChildrenByNode_.push_back(0);
        alive_.push_back(1);
        properHead_.push_back(InvalidNode);
        properTail_.push_back(InvalidNode);
        numProperPartsByNode_.push_back(0);
        return slotId;
    }

    /**
     * @brief Resets the proper-part linked-list storage to the requested size.
     *
     * @param numProperParts Proper-part data represented by `numProperParts`.
     */
    inline void initializeProperPartStorage(size_t numProperParts) {
        nextProperPart_.assign(numProperParts, InvalidNode);
        prevProperPart_.assign(numProperParts, InvalidNode);
    }

    /**
     * @brief Drops the preserved imported Higra id-domain marker after topology changes.
     */
    inline void invalidateHigraNodeIdSpace() noexcept { preservedExternalNodeIdOffset_.reset(); }

    /**
     * @brief Preserves an affine external-id mapping for Higra interoperability.
     *
     * Concrete import adapters own the external layout policy. The generic
     * topology only retains the offset needed to map its dense internal ids
     * while the structure remains unchanged.
     *
     * @param internalNodeOffset Offset of the internal-node segment in the imported domain.
     */
    inline void preserveExternalNodeIdOffset(NodeId internalNodeOffset) {
        if (internalNodeOffset < 0) {
            throw std::invalid_argument("An external internal-node offset must be non-negative.");
        }
        preservedExternalNodeIdOffset_ = internalNodeOffset;
    }

    /**
     * @brief Opens a staged edit session and rejects nested sessions.
     */
    inline void beginEditSession() {
        if (editSessionOpen_) {
            throw std::logic_error("A MorphologicalTree edit session is already open.");
        }
        editSessionOpen_ = true;
    }

    /**
     * @brief Closes the current staged edit session flag.
     */
    inline void endEditSession() noexcept { editSessionOpen_ = false; }

    /**
     * @brief Records edit commit.
     *
     * @param mode Mode selected for the operation.
     */
    inline void recordEditCommit(TreeEditValidationMode mode) noexcept {
        if (mode == TreeEditValidationMode::Incremental) {
            ++editValidationStatistics_.incrementalValidationCommits;
        } else {
            ++editValidationStatistics_.completeValidationCommits;
        }
    }

    /**
     * @brief Clears the current topology and keeps only an empty proper-part domain.
     *
     * @param numProperParts Proper-part data represented by `numProperParts`.
     */
    inline void initializeEmptyStorage(size_t numProperParts) {
        rootNodeId_ = InvalidNode;
        numNodes_ = 0;
        properPartOwner_.assign(numProperParts, InvalidNode);

        nodeParent_.clear();
        firstChild_.clear();
        nextSibling_.clear();
        prevSibling_.clear();
        lastChild_.clear();
        numChildrenByNode_.clear();
        alive_.clear();
        freeNodeIds_.clear();
        properHead_.clear();
        properTail_.clear();
        numProperPartsByNode_.clear();
        initializeProperPartStorage(numProperParts);
        invalidateHigraNodeIdSpace();

        prePostOrderCache_.timePreOrder.clear();
        prePostOrderCache_.timePostOrder.clear();
        invalidatePrePostOrderCache();
        invalidateAllIterators();
    }

    /**
     * @brief Clears one slot and returns it to the free-slot pool.
     *
     * @param slotId Dense internal-node slot identifier.
     */
    inline void releaseSlotStorage(NodeId slotId) {
        // Grow the free-list before changing the slot. If allocation fails, the
        // live topology is therefore left untouched instead of terminating from
        // inside a falsely-noexcept operation or publishing a half-released slot.
        freeNodeIds_.push_back(slotId);
        nodeParent_[slotId] = InvalidNode;
        firstChild_[slotId] = InvalidNode;
        nextSibling_[slotId] = InvalidNode;
        prevSibling_[slotId] = InvalidNode;
        lastChild_[slotId] = InvalidNode;
        numChildrenByNode_[slotId] = 0;
        alive_[slotId] = 0;
        properHead_[slotId] = InvalidNode;
        properTail_[slotId] = InvalidNode;
        numProperPartsByNode_[slotId] = 0;
    }

    /**
     * @brief Tests whether a dense slot currently belongs to the free list.
     *
     * @param slotId Dense internal-node slot identifier.
     * @return True if a dense slot currently belongs to the free list; otherwise false.
     */
    inline bool isFreeSlot(NodeId slotId) const noexcept { return slotId >= 0 && slotId < static_cast<NodeId>(alive_.size()) && alive_[slotId] == 0; }

    /**
     * @brief Rejects invalid or released internal-node ids before indexed reads.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @param context Operation context or diagnostic label.
     */
    inline void requireAliveNode(NodeId nodeId, const char* context) const {
        if (!isAlive(nodeId)) {
            throw std::invalid_argument(std::string(context) + " requires a live internal NodeId.");
        }
    }

    /**
     * @brief Rejects invalid, released, or root node ids for non-root-only edits.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @param context Operation context or diagnostic label.
     */
    inline void requireAliveNonRootNode(NodeId nodeId, const char* context) const {
        requireAliveNode(nodeId, context);
        if (isRoot(nodeId)) {
            throw std::invalid_argument(std::string(context) + " cannot target the root node.");
        }
    }

    /**
     * @brief Reconstructs direct proper-part linked lists from `properPartOwner_`.
     */
    inline void rebuildProperPartLinksFromOwnership() {
        properHead_.assign(nodeParent_.size(), InvalidNode);
        properTail_.assign(nodeParent_.size(), InvalidNode);
        numProperPartsByNode_.assign(nodeParent_.size(), 0);
        initializeProperPartStorage(properPartOwner_.size());

        for (NodeId properPartId = 0; properPartId < static_cast<NodeId>(properPartOwner_.size()); ++properPartId) {
            const NodeId ownerSlotId = properPartOwner_[properPartId];
            if (ownerSlotId == InvalidNode || ownerSlotId >= static_cast<NodeId>(nodeParent_.size()) || isFreeSlot(ownerSlotId)) {
                continue;
            }

            const NodeId tailProperPartId = properTail_[ownerSlotId];
            if (tailProperPartId == InvalidNode) {
                properHead_[ownerSlotId] = properPartId;
                properTail_[ownerSlotId] = properPartId;
            } else {
                nextProperPart_[tailProperPartId] = properPartId;
                prevProperPart_[properPartId] = tailProperPartId;
                properTail_[ownerSlotId] = properPartId;
            }
            numProperPartsByNode_[ownerSlotId]++;
        }
    }

    /**
     * @brief Removes one direct proper part from the linked list of its current owner.
     *
     * @param ownerSlotId Dense slot of the owning node.
     * @param properPartId Proper-part identifier used by the operation.
     */
    inline void unlinkProperPartFromOwner(NodeId ownerSlotId, NodeId properPartId) noexcept {
        const NodeId prev = prevProperPart_[static_cast<size_t>(properPartId)];
        const NodeId next = nextProperPart_[static_cast<size_t>(properPartId)];

        if (prev == InvalidNode) {
            properHead_[static_cast<size_t>(ownerSlotId)] = next;
        } else {
            nextProperPart_[static_cast<size_t>(prev)] = next;
        }

        if (next == InvalidNode) {
            properTail_[static_cast<size_t>(ownerSlotId)] = prev;
        } else {
            prevProperPart_[static_cast<size_t>(next)] = prev;
        }

        nextProperPart_[static_cast<size_t>(properPartId)] = InvalidNode;
        prevProperPart_[static_cast<size_t>(properPartId)] = InvalidNode;
        --numProperPartsByNode_[static_cast<size_t>(ownerSlotId)];
    }

    /**
     * @brief Appends one detached proper part to the target owner's direct list.
     *
     * @param targetSlotId Destination represented by `targetSlotId`.
     * @param properPartId Proper-part identifier used by the operation.
     */
    inline void appendDetachedProperPart(NodeId targetSlotId, NodeId properPartId) noexcept {
        properPartOwner_[static_cast<size_t>(properPartId)] = targetSlotId;
        const NodeId tail = properTail_[static_cast<size_t>(targetSlotId)];
        if (tail == InvalidNode) {
            properHead_[static_cast<size_t>(targetSlotId)] = properPartId;
            properTail_[static_cast<size_t>(targetSlotId)] = properPartId;
        } else {
            nextProperPart_[static_cast<size_t>(tail)] = properPartId;
            prevProperPart_[static_cast<size_t>(properPartId)] = tail;
            properTail_[static_cast<size_t>(targetSlotId)] = properPartId;
        }
        ++numProperPartsByNode_[static_cast<size_t>(targetSlotId)];
    }

    /**
     * @brief Splices all direct proper parts from `sourceSlotId` to the tail of `targetSlotId`.
     *
     * @param targetSlotId Destination represented by `targetSlotId`.
     * @param sourceSlotId Input represented by `sourceSlotId`.
     */
    inline void spliceProperPartsSlots(NodeId targetSlotId, NodeId sourceSlotId) noexcept {
        const NodeId sourceHead = properHead_[static_cast<size_t>(sourceSlotId)];
        if (sourceHead == InvalidNode) {
            return;
        }

        for (NodeId properPartId = sourceHead; properPartId != InvalidNode; properPartId = nextProperPart_[static_cast<size_t>(properPartId)]) {
            properPartOwner_[static_cast<size_t>(properPartId)] = targetSlotId;
        }

        const NodeId sourceTail = properTail_[static_cast<size_t>(sourceSlotId)];
        const int sourceCount = numProperPartsByNode_[static_cast<size_t>(sourceSlotId)];
        const NodeId targetTail = properTail_[static_cast<size_t>(targetSlotId)];

        if (targetTail == InvalidNode) {
            properHead_[static_cast<size_t>(targetSlotId)] = sourceHead;
            properTail_[static_cast<size_t>(targetSlotId)] = sourceTail;
        } else {
            nextProperPart_[static_cast<size_t>(targetTail)] = sourceHead;
            prevProperPart_[static_cast<size_t>(sourceHead)] = targetTail;
            properTail_[static_cast<size_t>(targetSlotId)] = sourceTail;
        }

        numProperPartsByNode_[static_cast<size_t>(targetSlotId)] += sourceCount;
        properHead_[static_cast<size_t>(sourceSlotId)] = InvalidNode;
        properTail_[static_cast<size_t>(sourceSlotId)] = InvalidNode;
        numProperPartsByNode_[static_cast<size_t>(sourceSlotId)] = 0;
    }

    /**
     * @brief Releases a live detached slot and invalidates dependent caches.
     *
     * @param slotNodeId Node identifier represented by `slotNodeId`.
     */
    inline void releaseSlotNode(NodeId slotNodeId) {
        releaseSlotStorage(slotNodeId);
        numNodes_--;
        invalidatePrePostOrderCache();
        bumpNodeStructureVersion();
    }

    /**
     * @brief Invalidates iterator-version counters after a full rebuild.
     */
    inline void invalidateAllIterators() noexcept {
        nodeStructureVersion_ = 0;
        topologyVersion_ = 0;
        properPartVersion_ = 0;
        ++mutationVersion_;
        lcaCache_.reset();
    }

    /**
     * @brief Drops the cached LCA structure after a topology-changing update.
     */
    inline void invalidateLcaCache() const noexcept { lcaCache_.reset(); }

    /**
     * @brief Marks alive-node iterators as invalid after node-structure changes.
     */
    inline void bumpNodeStructureVersion() noexcept {
        ++nodeStructureVersion_;
        ++mutationVersion_;
        invalidateLcaCache();
        invalidateHigraNodeIdSpace();
    }

    /**
     * @brief Marks topology-derived traversals and LCA queries as stale.
     */
    inline void bumpTopologyVersion() noexcept {
        ++topologyVersion_;
        ++mutationVersion_;
        invalidateLcaCache();
        invalidatePrePostOrderCache();
        invalidateHigraNodeIdSpace();
    }

    /**
     * @brief Marks direct proper-part iterators as stale.
     */
    inline void bumpProperPartVersion() noexcept {
        ++properPartVersion_;
        ++mutationVersion_;
        invalidateHigraNodeIdSpace();
    }

    /**
     * @brief Debug-only fail-fast check for alive-node iterators.
     *
     * @param expectedVersion Mutation version required by the operation.
     */
    inline void checkNodeIteratorVersion([[maybe_unused]] std::size_t expectedVersion) const {
        assert(expectedVersion == nodeStructureVersion_ && "Alive-node iterator invalidated by node-structure mutation.");
    }

    /**
     * @brief Debug-only fail-fast check for topology iterators.
     *
     * @param expectedVersion Mutation version required by the operation.
     */
    inline void checkTopologyIteratorVersion([[maybe_unused]] std::size_t expectedVersion) const {
        assert(expectedVersion == topologyVersion_ && "Topology iterator invalidated by tree-structure mutation.");
    }

    /**
     * @brief Debug-only fail-fast check for direct proper-part iterators.
     *
     * @param expectedVersion Mutation version required by the operation.
     */
    inline void checkProperPartIteratorVersion([[maybe_unused]] std::size_t expectedVersion) const {
        assert(expectedVersion == properPartVersion_ && "Proper-parts iterator invalidated by proper-part mutation.");
    }

    /**
     * @brief Invalidates preorder/postorder timestamps used by ancestry queries.
     */
    inline void invalidatePrePostOrderCache() const noexcept { prePostOrderCache_.invalidate(); }

    /**
     * @brief Recomputes preorder/postorder timestamps for the connected tree.
     */
    inline void recomputePrePostOrderCache() const {
        prePostOrderCache_.timePreOrder.assign(nodeParent_.size(), -1);
        prePostOrderCache_.timePostOrder.assign(nodeParent_.size(), -1);

        if (rootNodeId_ == InvalidNode) {
            prePostOrderCache_.valid = true;
            return;
        }

        int timer = 0;
        computeIncrementalAttributes(
            const_cast<MorphologicalTree*>(this), getRoot(), [&](NodeId nodeId) -> void { prePostOrderCache_.timePreOrder[nodeId] = timer++; },
            [&](NodeId, NodeId) -> void {}, [&](NodeId nodeId) -> void { prePostOrderCache_.timePostOrder[nodeId] = timer++; });

        prePostOrderCache_.valid = true;
    }

    /**
     * @brief Ensures that preorder/postorder timestamps are available.
     */
    inline void ensurePrePostOrderCache() const {
        if (!prePostOrderCache_.valid) {
            recomputePrePostOrderCache();
        }
    }

    /**
     * @brief Lazily builds the Euler-tour LCA cache on first use.
     *
     * @return Reference to the resulting object.
     */
    inline const LCAEulerRMQ& ensureLcaCache() const {
        if (!lcaCache_) {
            lcaCache_ = std::make_unique<LCAEulerRMQ>(this);
        }
        return *lcaCache_;
    }

    /**
     * @brief Generic depth-first accumulation helper used by incremental computations.
     *
     * @param tree Tree topology used by the operation.
     * @param rootNodeId Identifier of the traversal root.
     * @param preProcessing Callback invoked before visiting a node children.
     * @param mergeProcessing Callback invoked after completing one child.
     * @param postProcessing Callback invoked after all children are merged.
     */
    template <class PreProcessing, class MergeProcessing, class PostProcessing>
    static void computeIncrementalAttributes(MorphologicalTree* tree, NodeId rootNodeId, PreProcessing&& preProcessing, MergeProcessing&& mergeProcessing,
                                             PostProcessing&& postProcessing) {
        preProcessing(rootNodeId);
        for (NodeId childNodeId : tree->getChildren(rootNodeId)) {
            computeIncrementalAttributes(tree, childNodeId, preProcessing, mergeProcessing, postProcessing);
            mergeProcessing(rootNodeId, childNodeId);
        }
        postProcessing(rootNodeId);
    }

    // ========================= Internal topology helpers ========================= //
    // These helpers mutate the dense node-indexed topology storage directly.
    /**
     * @brief Appends `childId` to the end of the child list of `parentSlotId`.
     *
     * @param parentSlotId Parent-node value represented by `parentSlotId`.
     * @param childId Identifier of the child node.
     */
    void linkChildSlot(NodeId parentSlotId, NodeId childId) {
        if (parentSlotId < 0 || childId < 0)
            return;

        // If the child already has a parent, detach it before reading sibling links.
        if (nodeParent_[childId] != InvalidNode) {
            unlinkChildSlot(nodeParent_[childId], childId);
        }

        nodeParent_[childId] = parentSlotId;
        prevSibling_[childId] = lastChild_[parentSlotId];
        nextSibling_[childId] = InvalidNode;

        if (firstChild_[parentSlotId] == InvalidNode) {
            firstChild_[parentSlotId] = lastChild_[parentSlotId] = childId;
        } else {
            nextSibling_[lastChild_[parentSlotId]] = childId;
            lastChild_[parentSlotId] = childId;
        }
        ++numChildrenByNode_[parentSlotId];
        bumpTopologyVersion();
    }

    /**
     * @brief Removes one child from its parent list.
     *
     * @param parentSlotId Parent-node value represented by `parentSlotId`.
     * @param childId Identifier of the child node.
     */
    inline void unlinkChildSlot(NodeId parentSlotId, NodeId childId) {
        if (parentSlotId < 0 || childId < 0)
            return;
        if (nodeParent_[childId] != parentSlotId)
            return;

        const NodeId prev = prevSibling_[childId];
        const NodeId next = nextSibling_[childId];

        if (prev == InvalidNode)
            firstChild_[parentSlotId] = next;
        else
            nextSibling_[prev] = next;

        if (next == InvalidNode)
            lastChild_[parentSlotId] = prev;
        else
            prevSibling_[next] = prev;

        if (numChildrenByNode_[parentSlotId] > 0)
            --numChildrenByNode_[parentSlotId];

        nodeParent_[childId] = InvalidNode;
        prevSibling_[childId] = InvalidNode;
        nextSibling_[childId] = InvalidNode;
        bumpTopologyVersion();
    }

    /**
     * @brief Moves all children of `fromId` into the child list of `toId`.
     *
     * @param toId Destination node identifier.
     * @param fromId Source node identifier.
     * @param policy Policy controlling the operation.
     */
    inline void spliceChildrenSlots(NodeId toId, NodeId fromId, ChildSplicePolicy policy = ChildSplicePolicy::AppendToTargetTail) {
        if (toId < 0 || fromId < 0 || toId == fromId)
            return;

        const NodeId firstFrom = firstChild_[fromId];
        const NodeId lastFrom = lastChild_[fromId];
        const int movedCount = numChildrenByNode_[fromId];
        const bool replaceSourceSlot = policy == ChildSplicePolicy::ReplaceSourceSlotWhenDirectChild && nodeParent_[fromId] == toId;

        for (NodeId childId = firstFrom; childId != InvalidNode; childId = nextSibling_[childId]) {
            nodeParent_[childId] = toId;
        }

        if (replaceSourceSlot) {
            const NodeId prev = prevSibling_[fromId];
            const NodeId next = nextSibling_[fromId];

            if (firstFrom == InvalidNode) {
                if (prev == InvalidNode) {
                    firstChild_[toId] = next;
                } else {
                    nextSibling_[prev] = next;
                }

                if (next == InvalidNode) {
                    lastChild_[toId] = prev;
                } else {
                    prevSibling_[next] = prev;
                }
            } else {
                if (prev == InvalidNode) {
                    firstChild_[toId] = firstFrom;
                    prevSibling_[firstFrom] = InvalidNode;
                } else {
                    nextSibling_[prev] = firstFrom;
                    prevSibling_[firstFrom] = prev;
                }

                if (next == InvalidNode) {
                    lastChild_[toId] = lastFrom;
                    nextSibling_[lastFrom] = InvalidNode;
                } else {
                    prevSibling_[next] = lastFrom;
                    nextSibling_[lastFrom] = next;
                }
            }

            nodeParent_[fromId] = InvalidNode;
            prevSibling_[fromId] = InvalidNode;
            nextSibling_[fromId] = InvalidNode;
            numChildrenByNode_[toId] += movedCount - 1;
        } else {
            if (firstFrom == InvalidNode)
                return;

            if (firstChild_[toId] == InvalidNode) {
                firstChild_[toId] = firstFrom;
                lastChild_[toId] = lastFrom;
            } else {
                nextSibling_[lastChild_[toId]] = firstFrom;
                prevSibling_[firstFrom] = lastChild_[toId];
                lastChild_[toId] = lastFrom;
            }

            numChildrenByNode_[toId] += movedCount;
        }

        firstChild_[fromId] = InvalidNode;
        lastChild_[fromId] = InvalidNode;
        numChildrenByNode_[fromId] = 0;
        bumpTopologyVersion();
    }

    /**
     * @brief Releases an isolated detached node and returns its slot to the free pool.
     *
     * @param nodeId Identifier of the node used by the operation.
     */
    inline void releaseNode(NodeId nodeId) {
        if (!isNode(nodeId) || !isAlive(nodeId) || isRoot(nodeId)) {
            return;
        }
        const NodeId nodeSlot = nodeId;
        if (nodeParent_[nodeSlot] != nodeSlot) {
            return;
        }
        if (numChildrenByNode_[nodeSlot] != 0 || numProperPartsByNode_[nodeSlot] != 0) {
            return;
        }
        releaseSlotNode(nodeSlot);
    }

    /**
     * @brief Reuses one free slot as a live detached node without growing storage.
     *
     * @return Reuses one free slot as a live detached node without growing storage.
     */
    inline NodeId allocateNode() {
        if (freeNodeIds_.empty()) {
            return InvalidNode;
        }
        const NodeId nodeSlot = allocateSlot();
        nodeParent_[nodeSlot] = nodeSlot;
        numNodes_++;
        invalidatePrePostOrderCache();
        bumpNodeStructureVersion();
        return nodeSlot;
    }

    /**
     * @brief Creates a live detached node, reusing a free slot or appending a fresh one.
     *
     * Unlike `allocateNode()`, this method can grow the dense node domain when
     * no reusable slot is available. The created node starts self-parented and
     * therefore detached from the connected rooted component.
     *
     * @return The created live detached node, reusing a free slot or appending a fresh one.
     */
    inline NodeId createDetachedNode() {
        const NodeId nodeSlot = allocateSlot();
        nodeParent_[static_cast<size_t>(nodeSlot)] = nodeSlot;
        numNodes_++;
        invalidatePrePostOrderCache();
        bumpNodeStructureVersion();
        return nodeSlot;
    }

    /**
     * @brief Detaches `childId` from `parentNodeId` and optionally releases it.
     *
     * @param parentNodeId Identifier of the parent node.
     * @param childId Identifier of the child node.
     * @param releaseNodeFlag Flag controlling release node flag.
     */
    inline void removeChild(NodeId parentNodeId, NodeId childId, bool releaseNodeFlag) {
        if (!isAlive(parentNodeId) || !isAlive(childId) || !hasChild(parentNodeId, childId)) {
            return;
        }
        const NodeId parentSlotId = parentNodeId;
        const NodeId childSlotId = childId;
        if (releaseNodeFlag && numChildrenByNode_[static_cast<std::size_t>(childSlotId)] == 0 &&
            numProperPartsByNode_[static_cast<std::size_t>(childSlotId)] == 0) {
            freeNodeIds_.reserve(freeNodeIds_.size() + 1);
        }
        unlinkChildSlot(parentSlotId, childSlotId);
        nodeParent_[childSlotId] = childSlotId;
        invalidatePrePostOrderCache();
        if (releaseNodeFlag) {
            releaseNode(childId);
        }
    }

    /**
     * @brief Attaches `nodeId` as the last child of `parentNodeId`.
     *
     * @param parentNodeId Identifier of the parent node.
     * @param nodeId Identifier of the node used by the operation.
     */
    inline void attachNode(NodeId parentNodeId, NodeId nodeId) {
        if (!isAlive(parentNodeId) || !isAlive(nodeId) || isRoot(nodeId) || parentNodeId == nodeId) {
            return;
        }
        const NodeId parentSlotId = parentNodeId;
        const NodeId nodeSlotId = nodeId;
        const NodeId oldParentSlotId = nodeParent_[nodeSlotId];
        if (oldParentSlotId != InvalidNode && oldParentSlotId != nodeSlotId) {
            unlinkChildSlot(oldParentSlotId, nodeSlotId);
        }
        nodeParent_[nodeSlotId] = InvalidNode;
        linkChildSlot(parentSlotId, nodeSlotId);
    }

    /**
     * @brief Detaches `nodeId` from its current parent, leaving it self-parented.
     *
     * @param nodeId Identifier of the node used by the operation.
     */
    inline void detachNode(NodeId nodeId) {
        if (!isAlive(nodeId) || isRoot(nodeId)) {
            return;
        }
        const NodeId nodeSlotId = nodeId;
        const NodeId parentSlotId = nodeParent_[nodeSlotId];
        if (parentSlotId == InvalidNode || parentSlotId == nodeSlotId) {
            return;
        }
        unlinkChildSlot(parentSlotId, nodeSlotId);
        nodeParent_[nodeSlotId] = nodeSlotId;
    }

    /**
     * @brief Moves `nodeId` under `newParentId`.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @param newParentId Parent-node value represented by `newParentId`.
     */
    inline void moveNode(NodeId nodeId, NodeId newParentId) {
        if (!isAlive(nodeId) || !isAlive(newParentId) || isRoot(nodeId) || nodeId == newParentId) {
            return;
        }
        const NodeId nodeSlotId = nodeId;
        const NodeId newParentSlotId = newParentId;
        const NodeId oldParentSlotId = nodeParent_[nodeSlotId];
        if (oldParentSlotId == newParentSlotId) {
            return;
        }
        if (oldParentSlotId != InvalidNode && oldParentSlotId != nodeSlotId) {
            unlinkChildSlot(oldParentSlotId, nodeSlotId);
        }
        nodeParent_[nodeSlotId] = InvalidNode;
        linkChildSlot(newParentSlotId, nodeSlotId);
    }

    /**
     * @brief Transfers every direct child of `sourceId` to `parentNodeId`.
     *
     * @param parentNodeId Identifier of the parent node.
     * @param sourceId Input represented by `sourceId`.
     */
    inline void moveChildren(NodeId parentNodeId, NodeId sourceId) {
        if (!isAlive(parentNodeId) || !isAlive(sourceId) || parentNodeId == sourceId) {
            return;
        }
        spliceChildrenSlots(parentNodeId, sourceId);
    }

    /**
     * @brief Transfers one direct proper part from `sourceNodeId` to `targetNodeId`.
     *
     * @param targetNodeId Node identifier represented by `targetNodeId`.
     * @param sourceNodeId Node identifier represented by `sourceNodeId`.
     * @param properPartId Proper-part identifier used by the operation.
     */
    inline void moveProperPart(NodeId targetNodeId, NodeId sourceNodeId, NodeId properPartId) {
        if (!isAlive(targetNodeId) || !isAlive(sourceNodeId) || !isProperPart(properPartId) || targetNodeId == sourceNodeId) {
            return;
        }
        const NodeId sourceSlotId = sourceNodeId;
        if (properPartOwner_[properPartId] != sourceSlotId) {
            return;
        }
        unlinkProperPartFromOwner(sourceSlotId, properPartId);
        appendDetachedProperPart(targetNodeId, properPartId);
        bumpProperPartVersion();
    }

    /**
     * @brief Transfers all direct proper parts from `sourceNodeId` to `targetNodeId`.
     *
     * @param targetNodeId Node identifier represented by `targetNodeId`.
     * @param sourceNodeId Node identifier represented by `sourceNodeId`.
     */
    inline void moveProperParts(NodeId targetNodeId, NodeId sourceNodeId) {
        if (!isAlive(targetNodeId) || !isAlive(sourceNodeId) || targetNodeId == sourceNodeId) {
            return;
        }
        if (properHead_[static_cast<size_t>(sourceNodeId)] == InvalidNode) {
            return;
        }
        spliceProperPartsSlots(targetNodeId, sourceNodeId);
        bumpProperPartVersion();
    }

    /**
     * @brief Promotes `nodeId` to be the connected root.
     *
     * @param nodeId Identifier of the node used by the operation.
     */
    inline void setRoot(NodeId nodeId) {
        if (!isAlive(nodeId)) {
            return;
        }
        const NodeId nodeSlot = nodeId;
        if (rootNodeId_ == nodeSlot) {
            return;
        }
        if (rootNodeId_ != InvalidNode && !isFreeSlot(rootNodeId_)) {
            nodeParent_[rootNodeId_] = rootNodeId_;
            prevSibling_[rootNodeId_] = InvalidNode;
            nextSibling_[rootNodeId_] = InvalidNode;
        }
        const NodeId oldParentSlot = nodeParent_[nodeSlot];
        if (oldParentSlot != InvalidNode && oldParentSlot != nodeSlot) {
            unlinkChildSlot(oldParentSlot, nodeSlot);
        }
        rootNodeId_ = nodeSlot;
        nodeParent_[nodeSlot] = nodeSlot;
        prevSibling_[nodeSlot] = InvalidNode;
        nextSibling_[nodeSlot] = InvalidNode;
        bumpTopologyVersion();
    }

    /**
     * @brief Imports native node-parent and proper-part-owner buffers.
     *
     * This shared materialization path supports any builder that emits
     * independent node-parent and proper-part-owner domains, including live
     * nodes with no direct proper parts.
     *
     * @param nodeParent Parent-node value represented by `nodeParent`.
     * @param properPartOwner Proper-part data represented by `properPartOwner`.
     * @param root Root node of the traversal.
     * @param gridDomain2D Optional regular-grid proper-part domain.
     * @param semantics Hierarchy semantics validated by the operation.
     */
    void initializeNativeTopologyStorage(std::vector<NodeId> nodeParent, std::vector<NodeId> properPartOwner, NodeId root,
                                         std::optional<GridDomain2D> gridDomain2D, HierarchySemantics semantics) {
        if (gridDomain2D && properPartOwner.size() != gridDomain2D->size("Native topology 2D proper-part domain")) {
            throw std::invalid_argument("Native topology proper-part domain must match the attached 2D grid.");
        }
        if (nodeParent.size() > static_cast<std::size_t>(std::numeric_limits<NodeId>::max())) {
            throw std::invalid_argument("Native topology internal-node domain exceeds NodeId range.");
        }
        switch (semantics.descriptiveKind) {
        case MorphologicalTreeKind::GENERIC:
        case MorphologicalTreeKind::MAX_TREE:
        case MorphologicalTreeKind::MIN_TREE:
        case MorphologicalTreeKind::TREE_OF_SHAPES:
        case MorphologicalTreeKind::SELF_DUAL_RESIDUAL_TREE:
            break;
        default:
            throw std::invalid_argument("Native topology kind is not supported.");
        }
        const NodeId numNodeSlots = static_cast<NodeId>(nodeParent.size());
        const NodeId numProperParts = static_cast<NodeId>(properPartOwner.size());
        if (numNodeSlots <= 0) {
            throw std::invalid_argument("Native topology import requires at least one internal node.");
        }
        if (numProperParts <= 0) {
            throw std::invalid_argument("Native topology import requires at least one proper part.");
        }
        if (root < 0 || root >= numNodeSlots) {
            throw std::invalid_argument("Native topology import requires a valid root node id.");
        }
        semantics_ = std::move(semantics);
        gridDomain2D_ = gridDomain2D;
        const auto* uniformAdjacency = std::get_if<UniformGridAdjacency2D>(&semantics_.adjacency);
        if (uniformAdjacency && !gridDomain2D_) {
            throw std::invalid_argument("Native topology grid adjacency requires an attached 2D grid domain.");
        }
        if (uniformAdjacency &&
            (uniformAdjacency->relation.getNumRows() != gridDomain2D_->rows || uniformAdjacency->relation.getNumCols() != gridDomain2D_->cols)) {
            throw std::invalid_argument("Native topology adjacency must match the attached 2D grid.");
        }
        const auto* directionalAdjacency = std::get_if<DirectionalGridAdjacency2D>(&semantics_.adjacency);
        if (directionalAdjacency && !gridDomain2D_) {
            throw std::invalid_argument("Native topology directional grid adjacency requires an attached 2D grid domain.");
        }
        if (directionalAdjacency &&
            (directionalAdjacency->decreasing.getNumRows() != gridDomain2D_->rows || directionalAdjacency->decreasing.getNumCols() != gridDomain2D_->cols ||
             directionalAdjacency->increasing.getNumRows() != gridDomain2D_->rows || directionalAdjacency->increasing.getNumCols() != gridDomain2D_->cols)) {
            throw std::invalid_argument("Native topology directional adjacency context must match the attached 2D grid.");
        }

        initializeEmptyStorage(static_cast<size_t>(numProperParts));
        nodeParent_ = std::move(nodeParent);
        firstChild_.assign(static_cast<size_t>(numNodeSlots), InvalidNode);
        nextSibling_.assign(static_cast<size_t>(numNodeSlots), InvalidNode);
        prevSibling_.assign(static_cast<size_t>(numNodeSlots), InvalidNode);
        lastChild_.assign(static_cast<size_t>(numNodeSlots), InvalidNode);
        numChildrenByNode_.assign(static_cast<size_t>(numNodeSlots), 0);
        alive_.assign(static_cast<size_t>(numNodeSlots), 1);
        freeNodeIds_.clear();
        properHead_.assign(static_cast<size_t>(numNodeSlots), InvalidNode);
        properTail_.assign(static_cast<size_t>(numNodeSlots), InvalidNode);
        numProperPartsByNode_.assign(static_cast<size_t>(numNodeSlots), 0);
        properPartOwner_ = std::move(properPartOwner);
        initializeProperPartStorage(static_cast<size_t>(numProperParts));

        rootNodeId_ = root;
        numNodes_ = static_cast<int>(numNodeSlots);

        int selfParentedRoots = 0;
        for (NodeId nodeId = 0; nodeId < numNodeSlots; ++nodeId) {
            const NodeId parentId = nodeParent_[static_cast<size_t>(nodeId)];
            if (parentId == nodeId) {
                if (nodeId != root) {
                    throw std::invalid_argument("Native topology import found a detached self-parented non-root node.");
                }
                ++selfParentedRoots;
                continue;
            }
            if (parentId < 0 || parentId >= numNodeSlots) {
                throw std::invalid_argument("Native topology import found a parent outside the internal-node domain.");
            }

            prevSibling_[static_cast<size_t>(nodeId)] = lastChild_[static_cast<size_t>(parentId)];
            if (lastChild_[static_cast<size_t>(parentId)] == InvalidNode) {
                firstChild_[static_cast<size_t>(parentId)] = nodeId;
            } else {
                nextSibling_[static_cast<size_t>(lastChild_[static_cast<size_t>(parentId)])] = nodeId;
            }
            lastChild_[static_cast<size_t>(parentId)] = nodeId;
            ++numChildrenByNode_[static_cast<size_t>(parentId)];
        }
        if (selfParentedRoots != 1) {
            throw std::invalid_argument("Native topology import must encode exactly one self-parented root.");
        }

        rebuildProperPartLinksFromOwnership();
        invalidatePrePostOrderCache();
        invalidateAllIterators();
        preservedExternalNodeIdOffset_.reset();
    }

    /**
     * @brief Copies and fully validates a native topology input.
     *
     * @param nodeParent Parent-node value represented by `nodeParent`.
     * @param properPartOwner Proper-part data represented by `properPartOwner`.
     * @param root Root node of the traversal.
     * @param gridDomain2D Optional regular-grid proper-part domain.
     * @param semantics Hierarchy semantics validated by the operation.
     */
    void initializeNativeTopology(std::span<const NodeId> nodeParent, std::span<const NodeId> properPartOwner, NodeId root,
                                  std::optional<GridDomain2D> gridDomain2D, HierarchySemantics semantics) {
        if (nodeParent.size() > static_cast<std::size_t>(std::numeric_limits<NodeId>::max())) {
            throw std::invalid_argument("Native topology internal-node domain exceeds NodeId range.");
        }
        const NodeId numNodeSlots = static_cast<NodeId>(nodeParent.size());
        for (NodeId owner : properPartOwner) {
            if (owner < 0 || owner >= numNodeSlots) {
                throw std::invalid_argument("Native topology import found a proper-part owner outside the internal-node domain.");
            }
        }
        initializeNativeTopologyStorage(std::vector<NodeId>(nodeParent.begin(), nodeParent.end()),
                                        std::vector<NodeId>(properPartOwner.begin(), properPartOwner.end()), root, gridDomain2D, std::move(semantics));
        validateConnectedRootedTree();
    }

    /**
     * @brief Consumes producer-owned buffers paired with generic structural
     * evidence.
     *
     * @param nodeParent Parent-node value represented by `nodeParent`.
     * @param properPartOwner Proper-part data represented by `properPartOwner`.
     * @param root Root node of the traversal.
     * @param gridDomain2D Optional regular-grid proper-part domain.
     * @param semantics Hierarchy semantics validated by the operation.
     * @param topologyProof Proof that the native topology invariants hold.
     */
    void initializeValidatedNativeTopology(std::vector<NodeId>&& nodeParent, std::vector<NodeId>&& properPartOwner, NodeId root,
                                           std::optional<GridDomain2D> gridDomain2D, HierarchySemantics semantics,
                                           detail::NativeTopologyProof&& topologyProof) {
        topologyProof.requireMatches(nodeParent.size(), properPartOwner.size(), root);
        initializeNativeTopologyStorage(std::move(nodeParent), std::move(properPartOwner), root, gridDomain2D, std::move(semantics));
#ifndef NDEBUG
        validateConnectedRootedTree();
#endif
    }

    /**
     * @brief Factory-only materialization of a proven owning native topology.
     *
     * @param nodeParent Parent-node value represented by `nodeParent`.
     * @param properPartOwner Proper-part data represented by `properPartOwner`.
     * @param root Root node of the traversal.
     * @param gridDomain2D Optional regular-grid proper-part domain.
     * @param semantics Hierarchy semantics validated by the operation.
     * @param topologyProof Proof that the native topology invariants hold.
     */
    MorphologicalTree(detail::MorphologicalTreeConstructionTag, std::vector<NodeId>&& nodeParent, std::vector<NodeId>&& properPartOwner, NodeId root,
                      std::optional<GridDomain2D> gridDomain2D, HierarchySemantics semantics, detail::NativeTopologyProof&& topologyProof) {
        initializeValidatedNativeTopology(std::move(nodeParent), std::move(properPartOwner), root, gridDomain2D, std::move(semantics),
                                          std::move(topologyProof));
    }

    /**
     * @brief Private default constructor used by `clone()`.
     */
    MorphologicalTree() = default;

    /**
     * @brief Moves a committed topology and leaves `other` as an empty tree.
     *
     * The LCA cache is deliberately discarded because it stores a pointer to
     * the object that created it and cannot follow the topology to a new owner.
     *
     * @param other Object to compare with or transfer from.
     */
    void moveCommittedStateFrom(MorphologicalTree&& other) {
        rootNodeId_ = std::exchange(other.rootNodeId_, InvalidNode);
        semantics_ = std::move(other.semantics_);
        other.semantics_ = HierarchySemantics{};
        gridDomain2D_ = std::move(other.gridDomain2D_);
        other.gridDomain2D_.reset();
        numNodes_ = std::exchange(other.numNodes_, 0);
        preservedExternalNodeIdOffset_ = std::move(other.preservedExternalNodeIdOffset_);
        other.preservedExternalNodeIdOffset_.reset();
        editSessionOpen_ = false;
        editValidationStatistics_ = other.editValidationStatistics_;

        properPartOwner_ = std::move(other.properPartOwner_);
        nodeParent_ = std::move(other.nodeParent_);
        firstChild_ = std::move(other.firstChild_);
        nextSibling_ = std::move(other.nextSibling_);
        prevSibling_ = std::move(other.prevSibling_);
        lastChild_ = std::move(other.lastChild_);
        numChildrenByNode_ = std::move(other.numChildrenByNode_);
        alive_ = std::move(other.alive_);
        freeNodeIds_ = std::move(other.freeNodeIds_);
        properHead_ = std::move(other.properHead_);
        properTail_ = std::move(other.properTail_);
        numProperPartsByNode_ = std::move(other.numProperPartsByNode_);
        nextProperPart_ = std::move(other.nextProperPart_);
        prevProperPart_ = std::move(other.prevProperPart_);

        prePostOrderCache_ = std::move(other.prePostOrderCache_);
        other.prePostOrderCache_ = {};
        lcaCache_.reset();
        other.lcaCache_.reset();

        nodeStructureVersion_ = other.nodeStructureVersion_;
        topologyVersion_ = other.topologyVersion_;
        properPartVersion_ = other.properPartVersion_;
        mutationVersion_ = other.mutationVersion_;

        // Iterators into the moved-from owner must not accidentally validate
        // merely because its previous token happened to be zero.
        ++other.nodeStructureVersion_;
        ++other.topologyVersion_;
        ++other.properPartVersion_;
        ++other.mutationVersion_;
    }

  public:
    /**
     * @brief Copying is disabled to keep topology ownership explicit.
     */
    MorphologicalTree(const MorphologicalTree&) = delete;

    /**
     * @brief Copy assignment is disabled to keep topology ownership explicit.
     *
     * @return Reference to the resulting object.
     */
    MorphologicalTree& operator=(const MorphologicalTree&) = delete;

    /**
     * @brief Moves a complete committed topology.
     *
     * Moving an owner with an active editor would leave that editor pointing at
     * the moved-from object, so active edit sessions are rejected.
     *
     * @param other Object to compare with or transfer from.
     */
    MorphologicalTree(MorphologicalTree&& other) {
        other.requireNotEditing("MorphologicalTree move construction");
        moveCommittedStateFrom(std::move(other));
    }

    /**
     * @brief Move-assigns a complete committed topology.
     *
     * Both owners must be outside edit sessions so no editor can retain a
     * pointer to storage that changes owner.
     *
     * @param other Object to compare with or transfer from.
     * @return Reference to the resulting object.
     */
    MorphologicalTree& operator=(MorphologicalTree&& other) {
        requireNotEditing("MorphologicalTree move assignment destination");
        other.requireNotEditing("MorphologicalTree move assignment source");
        if (this != &other) {
            const std::size_t previousNodeVersion = nodeStructureVersion_;
            const std::size_t previousTopologyVersion = topologyVersion_;
            const std::size_t previousProperPartVersion = properPartVersion_;
            const std::size_t previousMutationVersion = mutationVersion_;
            moveCommittedStateFrom(std::move(other));
            nodeStructureVersion_ = std::max(nodeStructureVersion_, previousNodeVersion) + 1;
            topologyVersion_ = std::max(topologyVersion_, previousTopologyVersion) + 1;
            properPartVersion_ = std::max(properPartVersion_, previousProperPartVersion) + 1;
            mutationVersion_ = std::max(mutationVersion_, previousMutationVersion) + 1;
        }
        return *this;
    }

    /**
     * @brief Destroys the topology storage and cached traversal state.
     */
    virtual ~MorphologicalTree() = default;

    /**
     * @brief Tag-protected import from native MAF topology buffers.
     *
     * This path is used by builders that already materialize internal-node
     * parent links and row-major proper-part owners. Nodes are not required to
     * own a direct proper part, so the representation supports arbitrary
     * morphological trees of partial partitions.
     * The buffers are copied into canonical tree storage and validated as one
     * connected rooted hierarchy in which every node has non-empty subtree
     * support.
     *
     * @param nodeParent Parent-node value represented by `nodeParent`.
     * @param properPartOwner Proper-part data represented by `properPartOwner`.
     * @param root Root node of the traversal.
     * @param rows Number of rows in the domain.
     * @param cols Number of columns in the domain.
     * @param semantics Hierarchy semantics validated by the operation.
     */
    MorphologicalTree(detail::MorphologicalTreeConstructionTag, std::span<const NodeId> nodeParent, std::span<const NodeId> properPartOwner, NodeId root,
                      int rows, int cols, HierarchySemantics semantics) {
        initializeNativeTopology(nodeParent, properPartOwner, root, GridDomain2D{rows, cols}, std::move(semantics));
    }

    /**
     * @brief Imports a native hierarchy over an abstract finite proper-part set.
     *
     * No row/column interpretation is attached. Regular-grid reconstruction
     * and geometry-dependent algorithms consequently reject this tree
     * explicitly, while purely topological and support-based algorithms remain
     * available.
     *
     * @param nodeParent Parent-node value represented by `nodeParent`.
     * @param properPartOwner Proper-part data represented by `properPartOwner`.
     * @param root Root node of the traversal.
     * @param semantics Hierarchy semantics validated by the operation.
     */
    MorphologicalTree(detail::MorphologicalTreeConstructionTag, std::span<const NodeId> nodeParent, std::span<const NodeId> properPartOwner, NodeId root,
                      HierarchySemantics semantics) {
        initializeNativeTopology(nodeParent, properPartOwner, root, std::nullopt, std::move(semantics));
    }

    /**
     * @brief Creates an independent copy of the structural tree state.
     *
     * The public copy constructor stays deleted so ownership remains explicit at
     * API boundaries. This method is used by wrappers that need to preserve a
     * caller-owned topology while creating a new tree-backed object. A staged
     * edit cannot be cloned because its topology may be a temporary forest and
     * has not crossed the commit validation boundary.
     *
     * @return The created independent copy of the structural tree state.
     */
    MorphologicalTree clone() const {
        requireNotEditing("MorphologicalTree::clone");
        MorphologicalTree cloned;
        cloned.rootNodeId_ = rootNodeId_;
        cloned.semantics_ = semantics_;
        cloned.gridDomain2D_ = gridDomain2D_;
        cloned.numNodes_ = numNodes_;
        cloned.preservedExternalNodeIdOffset_ = preservedExternalNodeIdOffset_;
        cloned.editSessionOpen_ = false;
        cloned.editValidationStatistics_ = editValidationStatistics_;
        cloned.properPartOwner_ = properPartOwner_;
        cloned.nodeParent_ = nodeParent_;
        cloned.firstChild_ = firstChild_;
        cloned.nextSibling_ = nextSibling_;
        cloned.prevSibling_ = prevSibling_;
        cloned.lastChild_ = lastChild_;
        cloned.numChildrenByNode_ = numChildrenByNode_;
        cloned.alive_ = alive_;
        cloned.freeNodeIds_ = freeNodeIds_;
        cloned.properHead_ = properHead_;
        cloned.properTail_ = properTail_;
        cloned.numProperPartsByNode_ = numProperPartsByNode_;
        cloned.nextProperPart_ = nextProperPart_;
        cloned.prevProperPart_ = prevProperPart_;
        cloned.prePostOrderCache_ = prePostOrderCache_;
        cloned.lcaCache_.reset();
        cloned.nodeStructureVersion_ = nodeStructureVersion_;
        cloned.topologyVersion_ = topologyVersion_;
        cloned.properPartVersion_ = properPartVersion_;
        cloned.mutationVersion_ = mutationVersion_;
        return cloned;
    }

    // Forward declarations for nested iterator/range types whose definitions stay at the end of the class.
    class AliveNodeIterator;
    class AliveNodeRange;
    class ChildrenIterator;
    class ChildrenRange;
    class ProperPartsIterator;
    class ProperPartsRange;
    class ConnectedComponentIterator;
    class ConnectedComponentRange;
    class PostOrderNodeIterator;
    class PostOrderNodeRange;
    class BreadthFirstNodeIterator;
    class BreadthFirstNodeRange;
    class PathToRootIterator;
    class PathToRootRange;
    class PathBetweenNodesIterator;
    class PathBetweenNodesRange;
    class SubtreeNodeIterator;
    class SubtreeNodeRange;
    class DescendantNodeRange;

  public:
    // ========================= Public methods ========================= //

    /**
     * @brief Returns the monotonic mutation counter used by read-only views.
     *
     * @return The monotonic mutation counter used by read-only views.
     */
    std::size_t getMutationVersion() const noexcept { return mutationVersion_; }

    /**
     * @brief Rejects stale read-only views that captured an older mutation version.
     *
     * @param expectedVersion Mutation version required by the operation.
     * @param context Operation context or diagnostic label.
     */
    void requireMutationVersion(std::size_t expectedVersion, const char* context) const {
        if (mutationVersion_ != expectedVersion) {
            throw std::logic_error(std::string(context) + " cannot be used after the referenced tree topology has changed.");
        }
    }

    /**
     * @brief Returns the size of the dense internal-node id domain.
     *
     * Some slots may currently be free and therefore not correspond to live nodes.
     *
     * @return The size of the dense internal-node id domain.
     */
    inline int getNumInternalNodeSlots() const { return static_cast<int>(nodeParent_.size()); }

    /**
     * @brief Returns the size of the proper-part domain.
     *
     * @return The size of the proper-part domain.
     */
    inline int getNumTotalProperParts() const { return static_cast<int>(properPartOwner_.size()); }

    /**
     * @brief Returns the size of the preserved imported Higra node-id domain.
     *
     * @return The size of the preserved imported Higra node-id domain.
     *
     * @throws std::runtime_error if the tree was not imported from Higra or if
     * the preserved Higra node-id space was invalidated by an edit.
     *
     */
    inline int getNumHigraNodes() const {
        if (!preservedExternalNodeIdOffset_) {
            throw std::runtime_error("This tree does not preserve an imported Higra node-id space.");
        }
        return *preservedExternalNodeIdOffset_ + getNumInternalNodeSlots();
    }

    /**
     * @brief Returns the size of the requested node-id domain.
     *
     * `NodeIdSpace::HIGRA` means the preserved imported Higra domain, not the
     * compact domain that would be generated by exporting the current tree.
     *
     * @param outputSpace Node-id domain used to index the output.
     * @return The size of the requested node-id domain.
     */
    inline int getNodeIdSpaceSize(NodeIdSpace outputSpace) const {
        switch (outputSpace) {
        case NodeIdSpace::MORPHOLOGICAL_TREE:
            return getNumInternalNodeSlots();
        case NodeIdSpace::HIGRA:
            return getNumHigraNodes();
        }
        throw std::runtime_error("Unknown NodeIdSpace.");
    }

    /**
     * @brief Returns the preserved imported Higra node id for one live tree node.
     *
     * Returns `InvalidNode` when the original Higra node-id space is not
     * preserved.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @return The preserved imported Higra node id for one live tree node.
     */
    inline NodeId getHigraNodeId(NodeId nodeId) const noexcept {
        if (!preservedExternalNodeIdOffset_ || !isNode(nodeId) || !isAlive(nodeId)) {
            return InvalidNode;
        }
        return *preservedExternalNodeIdOffset_ + nodeId;
    }

    /**
     * @brief Returns the current hierarchy root.
     *
     * @return The current hierarchy root.
     */
    inline NodeId getRoot() const { return rootNodeId_; }

    /**
     * @brief Tests whether `nodeId` belongs to the internal-node id domain.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @return True if nodeId belongs to the internal-node id domain; otherwise false.
     */
    inline bool isNode(NodeId nodeId) const noexcept { return nodeId >= 0 && nodeId < static_cast<int>(nodeParent_.size()); }

    /**
     * @brief Tests whether `id` belongs to the proper-part domain.
     *
     * @param id Identifier used by the operation.
     * @return True if id belongs to the proper-part domain; otherwise false.
     */
    inline bool isProperPart(NodeId id) const noexcept { return id >= 0 && id < static_cast<int>(properPartOwner_.size()); }

    /**
     * @brief Tests whether a node slot currently represents a live node.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @return True if a node slot currently represents a live node; otherwise false.
     */
    inline bool isAlive(NodeId nodeId) const {
        if (!isNode(nodeId)) {
            return false;
        }
        const NodeId localId = nodeId;
        return localId >= 0 && localId < static_cast<NodeId>(nodeParent_.size()) && !isFreeSlot(localId) &&
               (nodeParent_[localId] != InvalidNode || localId == rootNodeId_);
    }

    /**
     * @brief Tests whether `nodeId` is the current root.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @return True if nodeId is the current root; otherwise false.
     */
    inline bool isRoot(NodeId nodeId) const { return nodeId == getRoot(); }

    /**
     * @brief Returns the number of currently reusable node slots.
     *
     * @return The number of currently reusable node slots.
     */
    inline int getNumFreeNodeSlots() const { return static_cast<int>(freeNodeIds_.size()); }

    /**
     * @brief Counts the live nodes that currently have no children.
     *
     * @return Counts the live nodes that currently have no children.
     */
    inline int getNumLeafNodes() const {
        int count = 0;
        for (NodeId nodeId : getAliveNodeIds()) {
            if (isLeaf(nodeId)) {
                ++count;
            }
        }
        return count;
    }

    /**
     * @brief Returns the number of direct children of `nodeId`.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @return The number of direct children of nodeId.
     */
    inline int getNumChildren(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getNumChildren");
        return numChildrenByNode_[nodeId];
    }

    /**
     * @brief Returns the number of internal descendants of `nodeId`.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @return The number of internal descendants of nodeId.
     */
    inline int getNodeNumDescendants(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getNodeNumDescendants");
        ensurePrePostOrderCache();
        const NodeId localId = nodeId;
        return (prePostOrderCache_.timePostOrder[localId] - prePostOrderCache_.timePreOrder[localId] - 1) / 2;
    }

    /**
     * @brief Returns the number of siblings of `nodeId`.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @return The number of siblings of nodeId.
     */
    inline int getNodeNumSiblings(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getNodeNumSiblings");
        if (isRoot(nodeId)) {
            return 0;
        }
        const NodeId parentNodeId = getNodeParent(nodeId);
        return std::max(0, getNumChildren(parentNodeId) - 1);
    }

    /**
     * @brief Returns the preorder time of `nodeId` in the cached DFS traversal.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @return The preorder time of nodeId in the cached DFS traversal.
     */
    inline int getNodeTimePreOrder(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getNodeTimePreOrder");
        ensurePrePostOrderCache();
        return prePostOrderCache_.timePreOrder[nodeId];
    }

    /**
     * @brief Returns the postorder time of `nodeId` in the cached DFS traversal.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @return The postorder time of nodeId in the cached DFS traversal.
     */
    inline int getNodeTimePostOrder(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getNodeTimePostOrder");
        ensurePrePostOrderCache();
        return prePostOrderCache_.timePostOrder[nodeId];
    }

    /**
     * @brief Returns the first direct child of `nodeId`, or `InvalidNode`.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @return The first direct child of nodeId, or InvalidNode.
     */
    inline NodeId getFirstChild(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getFirstChild");
        return firstChild_[nodeId];
    }

    /**
     * @brief Returns the next sibling of `nodeId`, or `InvalidNode`.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @return The next sibling of nodeId, or InvalidNode.
     */
    inline NodeId getNextSibling(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getNextSibling");
        return nextSibling_[nodeId];
    }

    /**
     * @brief Tests whether `nodeId` has no direct children.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @return True if nodeId has no direct children; otherwise false.
     */
    inline bool isLeaf(NodeId nodeId) const { return getFirstChild(nodeId) == InvalidNode; }

    /**
     * @brief Returns the number of direct proper parts owned by `nodeId`.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @return The number of direct proper parts owned by nodeId.
     */
    inline int getNumProperParts(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getNumProperParts");
        return numProperPartsByNode_[nodeId];
    }

    /**
     * @brief Tests whether a node exists only to represent hierarchy structure.
     *
     * Structural nodes own no direct proper parts. In a committed valid tree,
     * they still have non-empty full support through one or more descendants.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @return True if a node exists only to represent hierarchy structure; otherwise false.
     */
    inline bool isStructuralNode(NodeId nodeId) const { return getNumProperParts(nodeId) == 0; }

    /**
     * @brief Tests whether `childId` is a direct child of `parentNodeId`.
     *
     * @param parentNodeId Identifier of the parent node.
     * @param childId Identifier of the child node.
     * @return True if childId is a direct child of parentNodeId; otherwise false.
     */
    inline bool hasChild(NodeId parentNodeId, NodeId childId) const {
        requireAliveNode(parentNodeId, "MorphologicalTree::hasChild");
        requireAliveNode(childId, "MorphologicalTree::hasChild");
        return getNodeParent(childId) == parentNodeId;
    }

    /**
     * @brief Returns the direct parent of `nodeId`.
     *
     * The root and detached nodes report themselves as parent.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @return The direct parent of nodeId.
     */
    inline NodeId getNodeParent(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getNodeParent");
        if (isRoot(nodeId)) {
            return nodeId;
        }
        const NodeId parentSlot = nodeParent_[nodeId];
        if (parentSlot == InvalidNode) {
            return InvalidNode;
        }
        if (parentSlot == nodeId) {
            return nodeId;
        }
        return parentSlot;
    }

    /**
     * @brief Returns the live node that directly owns `properPartId`.
     *
     * @param properPartId Proper-part identifier used by the operation.
     * @return The live node that directly owns properPartId.
     */
    inline NodeId getProperPartOwner(NodeId properPartId) const {
        return isProperPart(properPartId) ? properPartOwner_[static_cast<size_t>(properPartId)] : InvalidNode;
    }

    /**
     * @brief Returns all live leaf nodes in the current hierarchy.
     *
     * @return All live leaf nodes in the current hierarchy.
     */
    std::vector<NodeId> getLeaves() const {
        std::vector<NodeId> leaves;
        if (rootNodeId_ == InvalidNode) {
            return leaves;
        }
        FastQueue<NodeId> s;
        s.push(this->rootNodeId_);

        while (!s.empty()) {
            const NodeId id = s.pop();
            if (numChildrenByNode_[id] == 0) {
                leaves.push_back(id);
            } else {
                for (NodeId c = firstChild_[id]; c != InvalidNode; c = nextSibling_[c]) {
                    s.push(c);
                }
            }
        }
        return leaves;
    }

    /**
     * @brief Returns the optional descriptive hierarchy-family label.
     *
     * Algorithms must use explicit capabilities such as `getAltitudeOrder()`
     * and `getAdjacencyMode()` instead of dispatching on this label.
     *
     * @return The optional descriptive hierarchy-family label.
     */
    inline MorphologicalTreeKind getDescriptiveKind() const noexcept { return semantics_.descriptiveKind; }

    /**
     * @brief Returns all generic semantic capabilities of this hierarchy.
     *
     * @return All generic semantic capabilities of this hierarchy.
     */
    inline const HierarchySemantics& getHierarchySemantics() const noexcept { return semantics_; }

    /**
     * @brief Returns the global parent-to-child altitude ordering constraint.
     *
     * @return The global parent-to-child altitude ordering constraint.
     */
    inline AltitudeOrder getAltitudeOrder() const noexcept { return semantics_.altitudeOrder; }

    /**
     * @brief Returns the shape of the attached adjacency context.
     *
     * @return The shape of the attached adjacency context.
     */
    inline AdjacencyMode getAdjacencyMode() const noexcept { return semantics_.adjacencyMode(); }

    /**
     * @brief Returns the number of currently live nodes.
     *
     * @return The number of currently live nodes.
     */
    inline int getNumNodes() const noexcept { return numNodes_; }

    /**
     * @brief Tests whether a staged edit session is currently open.
     *
     * @return True if a staged edit session is currently open; otherwise false.
     */
    inline bool isEditing() const noexcept { return editSessionOpen_; }

    /**
     * @brief Returns validation-path counters for committed edit sessions.
     *
     * @return Validation-path counters for committed edit sessions.
     */
    [[nodiscard]] const TreeEditValidationStatistics& getEditValidationStatistics() const noexcept { return editValidationStatistics_; }

    /**
     * @brief Rejects operations that require a committed connected topology.
     *
     * @param context Operation context or diagnostic label.
     */
    inline void requireNotEditing(const char* context) const {
        if (isEditing()) {
            throw std::logic_error(std::string(context) + " requires a committed MorphologicalTree; an edit session is still open.");
        }
    }

    /**
     * @brief Returns whether the tree currently contains alive detached nodes.
     *
     * Detached alive nodes are self-parented nodes other than the connected
     * root. They are allowed during explicit edit sessions but make the
     * hierarchy a forest rather than a single rooted tree.
     *
     * @return Whether the tree currently contains alive detached nodes.
     */
    inline bool hasDetachedAliveNodes() const noexcept {
        if (numNodes_ == 0) {
            return false;
        }
        if (rootNodeId_ == InvalidNode || !isAlive(rootNodeId_)) {
            return true;
        }

        for (NodeId nodeId = 0; nodeId < static_cast<NodeId>(nodeParent_.size()); ++nodeId) {
            if (!isAlive(nodeId) || nodeId == rootNodeId_) {
                continue;
            }
            if (nodeParent_[nodeId] == nodeId) {
                return true;
            }
        }
        return false;
    }

    /**
     * @brief Tests whether one uniform regular-grid 2D adjacency is attached.
     *
     * @return True if one uniform regular-grid 2D adjacency is attached; otherwise false.
     */
    inline bool hasUniformGridAdjacency2D() const noexcept { return std::holds_alternative<UniformGridAdjacency2D>(semantics_.adjacency); }

    /**
     * @brief Tests whether directional regular-grid 2D adjacency is attached.
     *
     * @return True if directional regular-grid 2D adjacency is attached; otherwise false.
     */
    inline bool hasDirectionalGridAdjacency2D() const noexcept { return std::holds_alternative<DirectionalGridAdjacency2D>(semantics_.adjacency); }

    /**
     * @brief Returns the directional adjacency context, or `nullptr`.
     *
     * @return The directional adjacency context, or nullptr.
     */
    inline const DirectionalGridAdjacency2D* getDirectionalGridAdjacency2D() const noexcept {
        return std::get_if<DirectionalGridAdjacency2D>(&semantics_.adjacency);
    }

    /**
     * @brief Returns the decreasing-arc adjacency relation, or `nullptr`.
     *
     * @return The decreasing-arc adjacency relation, or nullptr.
     */
    inline const RegularGridAdjacency2D* getDecreasingGridAdjacency2D() const noexcept {
        const auto* directional = getDirectionalGridAdjacency2D();
        return directional ? &directional->decreasing : nullptr;
    }

    /**
     * @brief Returns the increasing-arc adjacency relation, or `nullptr`.
     *
     * @return The increasing-arc adjacency relation, or nullptr.
     */
    inline const RegularGridAdjacency2D* getIncreasingGridAdjacency2D() const noexcept {
        const auto* directional = getDirectionalGridAdjacency2D();
        return directional ? &directional->increasing : nullptr;
    }

    /**
     * @brief Tests whether proper-part ids have an attached row/column layout.
     *
     * @return True if proper-part ids have an attached row/column layout; otherwise false.
     */
    [[nodiscard]] inline bool hasGridDomain2D() const noexcept { return gridDomain2D_.has_value(); }

    /**
     * @brief Returns the optional regular 2D proper-part domain.
     *
     * @return The optional regular 2D proper-part domain.
     */
    [[nodiscard]] inline const std::optional<GridDomain2D>& getGridDomain2D() const noexcept { return gridDomain2D_; }

    /**
     * @brief Returns the regular 2D domain or rejects a geometry-dependent call.
     *
     * @param context Operation context or diagnostic label.
     * @return The regular 2D domain or rejects a geometry-dependent call.
     */
    inline const GridDomain2D& requireGridDomain2D(const char* context) const {
        if (!gridDomain2D_) {
            throw std::invalid_argument(std::string(context) + " requires a regular 2D proper-part domain.");
        }
        return *gridDomain2D_;
    }

    /**
     * @brief Returns the number of rows in the regular 2D proper-part domain.
     *
     * @return The number of rows in the regular 2D proper-part domain.
     */
    inline int getNumRowsOfGridDomain2D() const { return requireGridDomain2D("MorphologicalTree::getNumRowsOfGridDomain2D").rows; }

    /**
     * @brief Returns the number of columns in the regular 2D proper-part domain.
     *
     * @return The number of columns in the regular 2D proper-part domain.
     */
    inline int getNumColsOfGridDomain2D() const { return requireGridDomain2D("MorphologicalTree::getNumColsOfGridDomain2D").cols; }

    /**
     * @brief Returns the immutable uniform regular-grid 2D adjacency.
     *
     * @return The immutable uniform regular-grid 2D adjacency.
     */
    inline const RegularGridAdjacency2D* getUniformGridAdjacency2D() const noexcept {
        const auto* uniform = std::get_if<UniformGridAdjacency2D>(&semantics_.adjacency);
        return uniform ? &uniform->relation : nullptr;
    }

    /**
     * @brief Validates one rooted tree with complete ownership and non-empty supports.
     *
     * This validation is intended for edit-session commits and debug/invariant
     * checks. It is intentionally stronger than the local runtime guards used
     * by low-level mutators.
     */
    void validateConnectedRootedTree() const {
        if (numNodes_ <= 0) {
            throw std::runtime_error("Connected-tree validation requires at least one live node.");
        }
        if (!isAlive(rootNodeId_)) {
            throw std::runtime_error("Connected-tree validation requires a live root.");
        }
        if (nodeParent_[rootNodeId_] != rootNodeId_) {
            throw std::runtime_error("Connected-tree validation requires the root to point to itself.");
        }

        const int numSlots = getNumInternalNodeSlots();
        const int numProperParts = getNumTotalProperParts();
        std::vector<int> expectedChildrenByNode(static_cast<size_t>(numSlots), 0);
        std::vector<int> expectedProperPartsByNode(static_cast<size_t>(numSlots), 0);
        int aliveCount = 0;

        for (NodeId nodeId = 0; nodeId < numSlots; ++nodeId) {
            if (!isAlive(nodeId)) {
                continue;
            }
            ++aliveCount;
            const NodeId parentNodeId = nodeParent_[static_cast<size_t>(nodeId)];
            if (nodeId == rootNodeId_) {
                continue;
            }
            if (parentNodeId == InvalidNode) {
                throw std::runtime_error("Connected-tree validation found an alive node with no parent.");
            }
            if (parentNodeId == nodeId) {
                throw std::runtime_error("Connected-tree validation found a detached alive node.");
            }
            if (!isNode(parentNodeId) || !isAlive(parentNodeId)) {
                throw std::runtime_error("Connected-tree validation found an alive node whose parent is outside the alive node domain.");
            }
            expectedChildrenByNode[static_cast<size_t>(parentNodeId)] += 1;
        }

        if (aliveCount != numNodes_) {
            throw std::runtime_error("Connected-tree validation found getNumNodes() out of sync with the alive node slots.");
        }

        std::vector<uint8_t> seenAsChild(static_cast<size_t>(numSlots), 0);
        for (NodeId nodeId = 0; nodeId < numSlots; ++nodeId) {
            if (!isAlive(nodeId)) {
                continue;
            }

            int actualChildCount = 0;
            NodeId previousChildId = InvalidNode;
            for (NodeId childNodeId = firstChild_[static_cast<size_t>(nodeId)]; childNodeId != InvalidNode;
                 childNodeId = nextSibling_[static_cast<size_t>(childNodeId)]) {
                if (!isNode(childNodeId) || !isAlive(childNodeId)) {
                    throw std::runtime_error("Connected-tree validation found a child list that references a non-alive node.");
                }
                if (nodeParent_[static_cast<size_t>(childNodeId)] != nodeId) {
                    throw std::runtime_error("Connected-tree validation found a child whose parent pointer disagrees with the child list.");
                }
                if (prevSibling_[static_cast<size_t>(childNodeId)] != previousChildId) {
                    throw std::runtime_error("Connected-tree validation found broken previous-sibling links.");
                }
                if (seenAsChild[static_cast<size_t>(childNodeId)] != 0) {
                    throw std::runtime_error("Connected-tree validation found a node referenced by multiple child lists.");
                }
                seenAsChild[static_cast<size_t>(childNodeId)] = 1;
                previousChildId = childNodeId;
                ++actualChildCount;
                if (actualChildCount > aliveCount) {
                    throw std::runtime_error("Connected-tree validation found a cycle in a child list.");
                }
            }

            if (firstChild_[static_cast<size_t>(nodeId)] == InvalidNode) {
                if (lastChild_[static_cast<size_t>(nodeId)] != InvalidNode) {
                    throw std::runtime_error("Connected-tree validation found an empty child list with a non-empty tail pointer.");
                }
            } else {
                if (lastChild_[static_cast<size_t>(nodeId)] != previousChildId) {
                    throw std::runtime_error("Connected-tree validation found an incorrect last-child pointer.");
                }
                if (nextSibling_[static_cast<size_t>(previousChildId)] != InvalidNode) {
                    throw std::runtime_error("Connected-tree validation found a child list whose tail still points to a next sibling.");
                }
            }

            if (actualChildCount != numChildrenByNode_[static_cast<size_t>(nodeId)]) {
                throw std::runtime_error("Connected-tree validation found an incorrect child count cache.");
            }
            if (actualChildCount != expectedChildrenByNode[static_cast<size_t>(nodeId)]) {
                throw std::runtime_error("Connected-tree validation found child lists out of sync with parent pointers.");
            }
        }

        if (seenAsChild[static_cast<size_t>(rootNodeId_)] != 0) {
            throw std::runtime_error("Connected-tree validation found the root inside a child list.");
        }
        for (NodeId nodeId = 0; nodeId < numSlots; ++nodeId) {
            if (!isAlive(nodeId) || nodeId == rootNodeId_) {
                continue;
            }
            if (seenAsChild[static_cast<size_t>(nodeId)] == 0) {
                throw std::runtime_error("Connected-tree validation found a non-root alive node missing from the child lists.");
            }
        }

        std::vector<NodeId> traversal;
        traversal.reserve(static_cast<size_t>(aliveCount));
        traversal.push_back(rootNodeId_);
        for (size_t index = 0; index < traversal.size(); ++index) {
            const NodeId nodeId = traversal[index];
            for (NodeId childId = firstChild_[static_cast<size_t>(nodeId)]; childId != InvalidNode; childId = nextSibling_[static_cast<size_t>(childId)]) {
                traversal.push_back(childId);
            }
        }
        if (traversal.size() != static_cast<size_t>(aliveCount)) {
            throw std::runtime_error("Connected-tree validation found nodes outside the rooted component or a cycle.");
        }

        for (int properPartId = 0; properPartId < numProperParts; ++properPartId) {
            const NodeId ownerNodeId = properPartOwner_[static_cast<size_t>(properPartId)];
            if (!isAlive(ownerNodeId)) {
                throw std::runtime_error("Connected-tree validation found a proper part without an alive owner.");
            }
            expectedProperPartsByNode[static_cast<size_t>(ownerNodeId)] += 1;
        }

        std::vector<uint8_t> seenProperPart(static_cast<size_t>(numProperParts), 0);
        for (NodeId nodeId = 0; nodeId < numSlots; ++nodeId) {
            if (!isAlive(nodeId)) {
                continue;
            }

            int actualProperPartCount = 0;
            NodeId previousProperPartId = InvalidNode;
            for (NodeId properPartId = properHead_[static_cast<size_t>(nodeId)]; properPartId != InvalidNode;
                 properPartId = nextProperPart_[static_cast<size_t>(properPartId)]) {
                if (!isProperPart(properPartId)) {
                    throw std::runtime_error("Connected-tree validation found an invalid proper-part id in the ownership list.");
                }
                if (properPartOwner_[static_cast<size_t>(properPartId)] != nodeId) {
                    throw std::runtime_error("Connected-tree validation found a proper-part list that disagrees with the ownership map.");
                }
                if (prevProperPart_[static_cast<size_t>(properPartId)] != previousProperPartId) {
                    throw std::runtime_error("Connected-tree validation found broken previous-proper-part links.");
                }
                if (seenProperPart[static_cast<size_t>(properPartId)] != 0) {
                    throw std::runtime_error("Connected-tree validation found a proper part referenced multiple times.");
                }
                seenProperPart[static_cast<size_t>(properPartId)] = 1;
                previousProperPartId = properPartId;
                ++actualProperPartCount;
                if (actualProperPartCount > numProperParts) {
                    throw std::runtime_error("Connected-tree validation found a cycle in a proper-part list.");
                }
            }

            if (properHead_[static_cast<size_t>(nodeId)] == InvalidNode) {
                if (properTail_[static_cast<size_t>(nodeId)] != InvalidNode) {
                    throw std::runtime_error("Connected-tree validation found an empty proper-part list with a non-empty tail pointer.");
                }
            } else {
                if (properTail_[static_cast<size_t>(nodeId)] != previousProperPartId) {
                    throw std::runtime_error("Connected-tree validation found an incorrect proper-part tail pointer.");
                }
                if (nextProperPart_[static_cast<size_t>(previousProperPartId)] != InvalidNode) {
                    throw std::runtime_error("Connected-tree validation found a proper-part list whose tail still points forward.");
                }
            }

            if (actualProperPartCount != numProperPartsByNode_[static_cast<size_t>(nodeId)]) {
                throw std::runtime_error("Connected-tree validation found an incorrect direct proper-part count cache.");
            }
            if (actualProperPartCount != expectedProperPartsByNode[static_cast<size_t>(nodeId)]) {
                throw std::runtime_error("Connected-tree validation found proper-part lists out of sync with the ownership map.");
            }
        }

        for (int properPartId = 0; properPartId < numProperParts; ++properPartId) {
            if (seenProperPart[static_cast<size_t>(properPartId)] == 0) {
                throw std::runtime_error("Connected-tree validation found a proper part missing from the ownership lists.");
            }
        }

        std::vector<int> subtreeSupport = expectedProperPartsByNode;
        for (auto it = traversal.rbegin(); it != traversal.rend(); ++it) {
            const NodeId nodeId = *it;
            if (subtreeSupport[static_cast<size_t>(nodeId)] == 0) {
                throw std::runtime_error("Connected-tree validation found a live node whose subtree support is empty.");
            }
            if (nodeId != rootNodeId_) {
                const NodeId parentId = nodeParent_[static_cast<size_t>(nodeId)];
                subtreeSupport[static_cast<size_t>(parentId)] += subtreeSupport[static_cast<size_t>(nodeId)];
            }
        }
    }

    /**
     * @brief Runs strong validation and returns the result instead of throwing.
     *
     * @return Result produced by running strong validation and returns the result instead of throwing.
     */
    [[nodiscard]] TreeValidationResult validateConnectedRootedTreeResult() const noexcept {
        try {
            validateConnectedRootedTree();
            return {true, ""};
        } catch (const std::exception& ex) {
            return {false, ex.what()};
        } catch (...) {
            return {false, "Connected-tree validation failed with an unknown error."};
        }
    }

    /**
     * @brief Tests whether `u` is an ancestor of `v`.
     *
     * @param u First endpoint or element used by the operation.
     * @param v Second endpoint or element used by the operation.
     * @return True if u is an ancestor of v; otherwise false.
     */
    inline bool isAncestor(NodeId u, NodeId v) const {
        requireAliveNode(u, "MorphologicalTree::isAncestor");
        requireAliveNode(v, "MorphologicalTree::isAncestor");
        ensurePrePostOrderCache();
        const NodeId slotU = u;
        const NodeId slotV = v;
        return prePostOrderCache_.timePreOrder[slotU] <= prePostOrderCache_.timePreOrder[slotV] &&
               prePostOrderCache_.timePostOrder[slotU] >= prePostOrderCache_.timePostOrder[slotV];
    }

    /**
     * @brief Tests whether `u` is a descendant of `v`.
     *
     * @param u First endpoint or element used by the operation.
     * @param v Second endpoint or element used by the operation.
     * @return True if u is a descendant of v; otherwise false.
     */
    inline bool isDescendant(NodeId u, NodeId v) const {
        requireAliveNode(u, "MorphologicalTree::isDescendant");
        requireAliveNode(v, "MorphologicalTree::isDescendant");
        ensurePrePostOrderCache();
        const NodeId slotU = u;
        const NodeId slotV = v;
        return prePostOrderCache_.timePreOrder[slotV] <= prePostOrderCache_.timePreOrder[slotU] &&
               prePostOrderCache_.timePostOrder[slotV] >= prePostOrderCache_.timePostOrder[slotU];
    }
    /**
     * @brief Tests whether `u` and `v` are comparable in the ancestry order.
     *
     * @param u First endpoint or element used by the operation.
     * @param v Second endpoint or element used by the operation.
     * @return True if u and v are comparable in the ancestry order; otherwise false.
     */
    inline bool isComparable(NodeId u, NodeId v) const { return isAncestor(u, v) || isAncestor(v, u); }

    /**
     * @brief Tests whether `u` is a strict ancestor of `v`.
     *
     * @param u First endpoint or element used by the operation.
     * @param v Second endpoint or element used by the operation.
     * @return True if u is a strict ancestor of v; otherwise false.
     */
    inline bool isStrictAncestor(NodeId u, NodeId v) const { return u != v && isAncestor(u, v); }

    /**
     * @brief Tests whether `u` is a strict descendant of `v`.
     *
     * @param u First endpoint or element used by the operation.
     * @param v Second endpoint or element used by the operation.
     * @return True if u is a strict descendant of v; otherwise false.
     */
    inline bool isStrictDescendant(NodeId u, NodeId v) const { return u != v && isDescendant(u, v); }

    /**
     * @brief Tests whether `u` and `v` are strictly comparable in the ancestry order.
     *
     * @param u First endpoint or element used by the operation.
     * @param v Second endpoint or element used by the operation.
     * @return True if u and v are strictly comparable in the ancestry order; otherwise false.
     */
    inline bool isStrictComparable(NodeId u, NodeId v) const { return isStrictAncestor(u, v) || isStrictAncestor(v, u); }

    /**
     * @brief Returns the lowest common ancestor of `u` and `v`.
     *
     * @param u First endpoint or element used by the operation.
     * @param v Second endpoint or element used by the operation.
     * @return The lowest common ancestor of u and v.
     */
    inline NodeId getLowestCommonAncestor(NodeId u, NodeId v) const {
        if (!isAlive(u) || !isAlive(v)) {
            return InvalidNode;
        }
        if (u == v) {
            return u;
        }
        if (isAncestor(u, v)) {
            return u;
        }
        if (isAncestor(v, u)) {
            return v;
        }
        if (getNodeParent(u) == u || getNodeParent(v) == v) {
            return InvalidNode;
        }
        return ensureLcaCache().findLowestCommonAncestor(u, v);
    }

    /**
     * @brief Returns a fail-fast range over all live node ids.
     *
     * @return A fail-fast range over all live node ids.
     */
    inline AliveNodeRange getAliveNodeIds() const { return AliveNodeRange(this, 0, getNumInternalNodeSlots(), nodeStructureVersion_); }

    /**
     * @brief Returns a fail-fast range over the direct children of `nodeId`.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @return A fail-fast range over the direct children of nodeId.
     */
    inline ChildrenRange getChildren(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getChildren");
        return ChildrenRange(this, firstChild_[nodeId], topologyVersion_);
    }

    /**
     * @brief Returns a fail-fast range over the direct proper parts of `nodeId`.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @return A fail-fast range over the direct proper parts of nodeId.
     */
    inline ProperPartsRange getProperParts(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getProperParts");
        return ProperPartsRange(this, properHead_[nodeId], properPartVersion_);
    }

    /**
     * @brief Returns a fail-fast range over all proper parts in the connected
     * component represented by `nodeId`.
     *
     * The range walks the subtree rooted at `nodeId` and yields every direct
     * proper part owned by those nodes, without materialising a vector.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @return A fail-fast range over all proper parts in the connected component represented by nodeId.
     */
    inline ConnectedComponentRange getConnectedComponent(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getConnectedComponent");
        return ConnectedComponentRange(this, nodeId, topologyVersion_, properPartVersion_);
    }

    /**
     * @brief Returns a post-order traversal range rooted at the connected root.
     *
     * @return A post-order traversal range rooted at the connected root.
     */
    inline PostOrderNodeRange getPostOrderNodes() const { return PostOrderNodeRange(this, getRoot(), topologyVersion_); }

    /**
     * @brief Returns a post-order traversal range rooted at `rootNodeId`.
     *
     * @param rootNodeId Identifier of the traversal root.
     * @return A post-order traversal range rooted at rootNodeId.
     */
    inline PostOrderNodeRange getPostOrderNodes(NodeId rootNodeId) const {
        requireAliveNode(rootNodeId, "MorphologicalTree::getPostOrderNodes");
        return PostOrderNodeRange(this, rootNodeId, topologyVersion_);
    }

    /**
     * @brief Returns a breadth-first traversal range rooted at the connected root.
     *
     * @return A breadth-first traversal range rooted at the connected root.
     */
    inline BreadthFirstNodeRange getIteratorBreadthFirstTraversal() const { return BreadthFirstNodeRange(this, getRoot(), topologyVersion_); }

    /**
     * @brief Returns a breadth-first traversal range rooted at `rootNodeId`.
     *
     * @param rootNodeId Identifier of the traversal root.
     * @return A breadth-first traversal range rooted at rootNodeId.
     */
    inline BreadthFirstNodeRange getIteratorBreadthFirstTraversal(NodeId rootNodeId) const {
        requireAliveNode(rootNodeId, "MorphologicalTree::getIteratorBreadthFirstTraversal");
        return BreadthFirstNodeRange(this, rootNodeId, topologyVersion_);
    }

    /**
     * @brief Returns the path from `nodeId` to the connected root.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @return The path from nodeId to the connected root.
     */
    inline PathToRootRange getPathToRootNodes(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getPathToRootNodes");
        return PathToRootRange(this, nodeId, topologyVersion_);
    }

    /**
     * @brief Returns the path that connects `sourceNodeId` and `targetNodeId`.
     *
     * @param sourceNodeId Node identifier represented by `sourceNodeId`.
     * @param targetNodeId Node identifier represented by `targetNodeId`.
     * @return The path that connects sourceNodeId and targetNodeId.
     */
    inline PathBetweenNodesRange getPathBetweenNodes(NodeId sourceNodeId, NodeId targetNodeId) const {
        requireAliveNode(sourceNodeId, "MorphologicalTree::getPathBetweenNodes");
        requireAliveNode(targetNodeId, "MorphologicalTree::getPathBetweenNodes");
        return PathBetweenNodesRange(this, sourceNodeId, targetNodeId, topologyVersion_);
    }

    /**
     * @brief Returns a pre-order traversal range over the subtree of `nodeId`.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @return A pre-order traversal range over the subtree of nodeId.
     */
    inline SubtreeNodeRange getNodeSubtree(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getNodeSubtree");
        return SubtreeNodeRange(this, nodeId, topologyVersion_);
    }

    /**
     * @brief Returns a range over all proper descendants of `nodeId`.
     *
     * @param nodeId Identifier of the node used by the operation.
     * @return A range over all proper descendants of nodeId.
     */
    inline DescendantNodeRange getDescendants(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getDescendants");
        return DescendantNodeRange(this, nodeId, topologyVersion_);
    }

    /**
     * @brief Opens the only public entrypoint for staged structural mutations.
     *
     * Use this for multi-step topology rewiring that may temporarily detach
     * nodes, move children/proper parts, or create intermediate nodes.
     * `TreeEditor::commit()` runs the strong connected-rooted-tree validation.
     *
     * Defined in `TreeEditor.hpp` after `TreeEditor` is complete; at this point
     * only the forward declaration is available.
     *
     * @return The opened only public entrypoint for staged structural mutations.
     */
    [[nodiscard]] TreeEditor edit();

    /**
     * @brief Prunes the subtree of `nodeId`, moving all its support to the parent.
     *
     * This is a safe committed edit: it is a complete local operation and does
     * not intentionally leave the tree in a staged disconnected state, but it
     * still advances the tree mutation version and invalidates derived state.
     *
     * @param nodeId Identifier of the node used by the operation.
     */
    inline void pruneNode(NodeId nodeId) {
        requireNotEditing("MorphologicalTree::pruneNode");
        requireAliveNonRootNode(nodeId, "MorphologicalTree::pruneNode");
        const NodeId parentNodeId = getNodeParent(nodeId);
        if (parentNodeId == InvalidNode || parentNodeId == nodeId) {
            throw std::invalid_argument("MorphologicalTree::pruneNode requires an attached non-root node.");
        }
        const NodeId parentSlotId = parentNodeId;
        std::size_t releaseCount = 0;
        for ([[maybe_unused]] NodeId subtreeNode : getNodeSubtree(nodeId)) {
            ++releaseCount;
        }
        freeNodeIds_.reserve(freeNodeIds_.size() + releaseCount);

        const auto descendToDeepestLastChild = [this](NodeId startId) {
            NodeId currentId = startId;
            while (firstChild_[currentId] != InvalidNode) {
                currentId = lastChild_[currentId];
            }
            return currentId;
        };

        NodeId currentId = descendToDeepestLastChild(nodeId);
        while (true) {
            const NodeId currentParentSlotId = nodeParent_[currentId];
            moveProperParts(parentSlotId, currentId);
            if (currentParentSlotId != InvalidNode && currentParentSlotId != currentId) {
                unlinkChildSlot(currentParentSlotId, currentId);
            }
            nodeParent_[currentId] = currentId;
            releaseNode(currentId);

            if (currentId == nodeId) {
                break;
            }
            currentId = firstChild_[currentParentSlotId] == InvalidNode ? currentParentSlotId : descendToDeepestLastChild(currentParentSlotId);
        }
    }

    /**
     * @brief Merges `nodeId` into its parent and releases the emptied slot.
     *
     * This is a safe committed edit: children and direct proper parts are
     * transferred to the parent before the node slot is released. The operation
     * advances the tree mutation version and invalidates derived state.
     *
     * @param nodeId Identifier of the node used by the operation.
     */
    inline void mergeNodeIntoParent(NodeId nodeId) {
        requireNotEditing("MorphologicalTree::mergeNodeIntoParent");
        requireAliveNonRootNode(nodeId, "MorphologicalTree::mergeNodeIntoParent");
        const NodeId parentNodeId = getNodeParent(nodeId);
        if (parentNodeId == InvalidNode || parentNodeId == nodeId) {
            throw std::invalid_argument("MorphologicalTree::mergeNodeIntoParent requires an attached non-root node.");
        }
        const NodeId parentSlotId = parentNodeId;
        const NodeId nodeSlotId = nodeId;
        freeNodeIds_.reserve(freeNodeIds_.size() + 1);

        spliceProperPartsSlots(parentSlotId, nodeSlotId);
        spliceChildrenSlots(parentSlotId, nodeSlotId, ChildSplicePolicy::ReplaceSourceSlotWhenDirectChild);
        nodeParent_[nodeSlotId] = nodeSlotId;
        bumpProperPartVersion();
        releaseNode(nodeId);
    }

  private:
    // ========================= Internal classes ========================= //
    /**
     * @brief Euler-tour + RMQ cache for lowest-common-ancestor queries.
     *
     * Stage 1: Euler tour
     * Performs a DFS on the tree and records:
     * 1. The order of visited nodes -> `euler_[]`
     * 2. The depth of each visited node -> `depth_[]`
     * 3. The first Euler-tour position of each node -> `firstOccurrence_[]`
     *
     * Stage 2: RMQ over `depth_[]`
     * To answer `LCA(u, v)`:
     * 1. Let `i = firstOccurrence_[u]` and `j = firstOccurrence_[v]`
     * 2. Run an RMQ over `depth_[min(i, j)..max(i, j)]`
     * 3. The Euler-tour position with minimum depth corresponds to the LCA
     *
     * Example:
     * ```text
     *   0
     *  / \
     * 1   2
     * |
     * 3
     *
     * Indices:          0  1  2  3  4  5  6
     * euler_ =          [0, 1, 3, 1, 0, 2, 0]
     * depth_ =          [0, 1, 2, 1, 0, 1, 0]
     * firstOccurrence_ =[0, 1, 5, 2]
     *
     * LCA(3, 2) = 0
     * i = firstOccurrence_[3] = 2
     * j = firstOccurrence_[2] = 5
     * RMQ over depth_[2..5] = [2, 1, 0, 1]
     * The minimum depth is 0 at Euler-tour position 4
     * euler_[4] = 0, therefore the LCA is node 0
     * ```
     *
     * The sparse table stores Euler-tour positions, not node ids directly. The
     * cache is built lazily on demand and invalidated whenever topology changes.
     */
    class LCAEulerRMQ {
      private:
        /** @brief Stores the euler. */
        std::vector<NodeId> euler_;
        /** @brief Stores the depth. */
        std::vector<int> depth_;
        /** @brief Stores the first occurrence. */
        std::vector<int> firstOccurrence_;
        /** @brief Stores the log2. */
        std::vector<int> log2_;
        /** @brief Stores the sparse table. */
        std::vector<int> sparseTable_;
        /** @brief Stores the sparse table stride. */
        int sparseTableStride_ = 0;
        /** @brief References the tree used by the component. */
        const MorphologicalTree* tree_ = nullptr;

        /**
         * @brief Appends one DFS step to the Euler tour and recurses into children.
         *
         * @param nodeId Identifier of the node used by the operation.
         * @param currentDepth Depth of the current traversal node.
         */
        void depthFirstTraversal(NodeId nodeId, int currentDepth) {
            const NodeId slotId = nodeId;
            if (firstOccurrence_[static_cast<size_t>(slotId)] == -1) {
                firstOccurrence_[static_cast<size_t>(slotId)] = static_cast<int>(euler_.size());
            }

            euler_.push_back(nodeId);
            depth_.push_back(currentDepth);

            for (NodeId childNodeId : tree_->getChildren(nodeId)) {
                depthFirstTraversal(childNodeId, currentDepth + 1);
                euler_.push_back(nodeId);
                depth_.push_back(currentDepth);
            }
        }

        /**
         * @brief Builds the sparse table used by the RMQ query on `depth_[]`.
         */
        void buildSparseTable() {
            const int n = static_cast<int>(depth_.size());
            if (n == 0) {
                log2_.clear();
                sparseTable_.clear();
                sparseTableStride_ = 0;
                return;
            }

            log2_.assign(static_cast<size_t>(n + 1), 0);
            for (int i = 2; i <= n; ++i) {
                log2_[static_cast<size_t>(i)] = log2_[static_cast<size_t>(i / 2)] + 1;
            }

            sparseTableStride_ = 1;
            while ((1 << sparseTableStride_) <= n) {
                ++sparseTableStride_;
            }

            sparseTable_.assign(static_cast<size_t>(n) * static_cast<size_t>(sparseTableStride_), 0);

            for (int i = 0; i < n; ++i) {
                sparseTable_[sparseTableIndex(i, 0)] = i;
            }

            for (int j = 1; j < sparseTableStride_; ++j) {
                const int blockSize = 1 << j;
                const int halfBlock = blockSize >> 1;
                for (int i = 0; i + blockSize <= n; ++i) {
                    const int leftIndex = sparseTable_[sparseTableIndex(i, j - 1)];
                    const int rightIndex = sparseTable_[sparseTableIndex(i + halfBlock, j - 1)];
                    sparseTable_[sparseTableIndex(i, j)] =
                        depth_[static_cast<size_t>(leftIndex)] <= depth_[static_cast<size_t>(rightIndex)] ? leftIndex : rightIndex;
                }
            }
        }

        /**
         * @brief Maps a sparse-table row/column pair to the flat storage index.
         *
         * @param row Zero-based row coordinate.
         * @param col Zero-based column coordinate.
         * @return The mapped sparse-table row/column pair to the flat storage index.
         */
        std::size_t sparseTableIndex(int row, int col) const {
            return static_cast<std::size_t>(row) * static_cast<std::size_t>(sparseTableStride_) + static_cast<std::size_t>(col);
        }

        /**
         * @brief Returns the Euler-tour position of minimum depth on `[left, right]`.
         *
         * @param left Left-hand sparse-table index.
         * @param right Right-hand sparse-table index.
         * @return The Euler-tour position of minimum depth on [left, right].
         */
        int rmq(int left, int right) const {
            const int length = right - left + 1;
            const int logLength = log2_[static_cast<size_t>(length)];
            const int leftIndex = sparseTable_[sparseTableIndex(left, logLength)];
            const int rightIndex = sparseTable_[sparseTableIndex(right - (1 << logLength) + 1, logLength)];
            return depth_[static_cast<size_t>(leftIndex)] <= depth_[static_cast<size_t>(rightIndex)] ? leftIndex : rightIndex;
        }

      public:
        /**
         * @brief Builds the Euler tour and the RMQ sparse table for one rooted tree.
         *
         * @param tree Tree topology used by the operation.
         */
        explicit LCAEulerRMQ(const MorphologicalTree* tree) : tree_(tree) {
            if (tree_ == nullptr || tree_->getRoot() == InvalidNode) {
                return;
            }

            euler_.reserve(static_cast<size_t>(std::max(0, 2 * tree_->getNumNodes() - 1)));
            depth_.reserve(static_cast<size_t>(std::max(0, 2 * tree_->getNumNodes() - 1)));
            firstOccurrence_.assign(static_cast<size_t>(tree_->getNumInternalNodeSlots()), -1);

            depthFirstTraversal(tree_->getRoot(), 0);
            buildSparseTable();
        }

        /**
         * @brief Returns the lowest common ancestor of `u` and `v`.
         *
         * @param u First endpoint or element used by the operation.
         * @param v Second endpoint or element used by the operation.
         * @return The lowest common ancestor of u and v.
         */
        NodeId findLowestCommonAncestor(NodeId u, NodeId v) const {
            const int firstU = firstOccurrence_[static_cast<size_t>(u)];
            const int firstV = firstOccurrence_[static_cast<size_t>(v)];
            if (firstU < 0 || firstV < 0) {
                return InvalidNode;
            }

            int left = firstU;
            int right = firstV;
            if (left > right) {
                std::swap(left, right);
            }
            return euler_[static_cast<size_t>(rmq(left, right))];
        }
    };

  public:
    // ========================= Internal classes and iterators ========================= //

    /**
     * @brief Iterator over live node ids in the dense internal-node domain.
     */
    class AliveNodeIterator {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Stores the current. */
        NodeId current_ = InvalidNode;
        /** @brief Stores the end. */
        NodeId end_ = InvalidNode;
        /** @brief Stores the expected version. */
        std::size_t expectedVersion_ = 0;

        /**
         * @brief Advances to the next live slot or to the end sentinel.
         */
        void settle_() {
            if (!T_) {
                current_ = InvalidNode;
                return;
            }
            T_->checkNodeIteratorVersion(expectedVersion_);
            while (current_ != InvalidNode && current_ < end_ && !T_->isAlive(current_)) {
                ++current_;
            }
            if (current_ >= end_) {
                current_ = InvalidNode;
            }
        }

      public:
        /// Standard iterator category exposed for STL compatibility.
        using iterator_category = std::input_iterator_tag;
        /// Node id yielded by the iterator.
        using value_type = NodeId;
        /// Signed distance type exposed for STL compatibility.
        using difference_type = std::ptrdiff_t;
        /// Pointer type exposed for STL compatibility.
        using pointer = const NodeId*;
        /// Reference type exposed for STL compatibility.
        using reference = const NodeId&;

        /**
         * @brief Creates the end/sentinel iterator.
         */
        AliveNodeIterator() = default;

        /**
         * @brief Creates an iterator over `[current, end)` with fail-fast version checking.
         *
         * @param tree Tree topology used by the operation.
         * @param current Current iterator or traversal position.
         * @param end Exclusive end position.
         * @param expectedVersion Mutation version required by the operation.
         */
        AliveNodeIterator(const MorphologicalTree* tree, NodeId current, NodeId end, std::size_t expectedVersion)
            : T_(tree), current_(current), end_(end), expectedVersion_(expectedVersion) {
            settle_();
        }

        /**
         * @brief Advances to the next live node slot.
         *
         * @return Reference to the resulting object.
         */
        AliveNodeIterator& operator++() {
            T_->checkNodeIteratorVersion(expectedVersion_);
            if (current_ != InvalidNode) {
                ++current_;
                settle_();
            }
            return *this;
        }

        /**
         * @brief Returns the current live node id.
         *
         * @return The current live node id.
         */
        NodeId operator*() const {
            T_->checkNodeIteratorVersion(expectedVersion_);
            return current_;
        }

        /**
         * @brief Compares iterator positions.
         *
         * @param other Object to compare with or transfer from.
         * @return True when the documented condition holds; otherwise false.
         */
        bool operator==(const AliveNodeIterator& other) const { return current_ == other.current_; }

        /**
         * @brief Compares iterator positions for inequality.
         *
         * @param other Object to compare with or transfer from.
         * @return True when the documented condition holds; otherwise false.
         */
        bool operator!=(const AliveNodeIterator& other) const { return !(*this == other); }
    };

    /**
     * @brief Range wrapper for iterating over live node ids.
     */
    class AliveNodeRange {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Stores the begin. */
        NodeId begin_ = InvalidNode;
        /** @brief Stores the end. */
        NodeId end_ = InvalidNode;
        /** @brief Stores the expected version. */
        std::size_t expectedVersion_ = 0;

      public:
        /**
         * @brief Creates an empty live-node range.
         */
        AliveNodeRange() = default;

        /**
         * @brief Creates a fail-fast live-node range over dense slots.
         *
         * @param tree Tree topology used by the operation.
         * @param begin Inclusive beginning of the range.
         * @param end Exclusive end position.
         * @param expectedVersion Mutation version required by the operation.
         */
        AliveNodeRange(const MorphologicalTree* tree, NodeId begin, NodeId end, std::size_t expectedVersion)
            : T_(tree), begin_(begin), end_(end), expectedVersion_(expectedVersion) {}

        /**
         * @brief Returns an iterator positioned at the first live slot.
         *
         * @return An iterator positioned at the first live slot.
         */
        AliveNodeIterator begin() const { return AliveNodeIterator(T_, begin_, end_, expectedVersion_); }

        /**
         * @brief Returns the sentinel iterator for the live-node range.
         *
         * @return The sentinel iterator for the live-node range.
         */
        AliveNodeIterator end() const { return AliveNodeIterator(T_, InvalidNode, end_, expectedVersion_); }
    };

    /**
     * @brief Iterator over the direct children of one node.
     */
    class ChildrenIterator {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Stores the current local. */
        NodeId currentLocal_ = InvalidNode;
        /** @brief Stores the expected version. */
        std::size_t expectedVersion_ = 0;

      public:
        /// Standard iterator category exposed for STL compatibility.
        using iterator_category = std::forward_iterator_tag;
        /// Node id yielded by the iterator.
        using value_type = NodeId;
        /// Signed distance type exposed for STL compatibility.
        using difference_type = std::ptrdiff_t;
        /// Pointer type exposed for STL compatibility.
        using pointer = const NodeId*;
        /// Reference type exposed for STL compatibility.
        using reference = const NodeId&;

        /**
         * @brief Creates the end/sentinel child iterator.
         */
        ChildrenIterator() = default;

        /**
         * @brief Creates a child iterator starting from a linked-list slot.
         *
         * @param tree Tree topology used by the operation.
         * @param currentLocal Current local iterator position.
         * @param expectedVersion Mutation version required by the operation.
         */
        ChildrenIterator(const MorphologicalTree* tree, NodeId currentLocal, std::size_t expectedVersion)
            : T_(tree), currentLocal_(currentLocal), expectedVersion_(expectedVersion) {}

        /**
         * @brief Advances to the next sibling in the child list.
         *
         * @return Reference to the resulting object.
         */
        ChildrenIterator& operator++() {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            if (T_ && currentLocal_ != InvalidNode) {
                currentLocal_ = T_->nextSibling_[currentLocal_];
            }
            return *this;
        }

        /**
         * @brief Returns the current child node id.
         *
         * @return The current child node id.
         */
        NodeId operator*() const {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            return currentLocal_;
        }

        /**
         * @brief Compares child iterator positions.
         *
         * @param other Object to compare with or transfer from.
         * @return True when the documented condition holds; otherwise false.
         */
        bool operator==(const ChildrenIterator& other) const { return currentLocal_ == other.currentLocal_; }

        /**
         * @brief Compares child iterator positions for inequality.
         *
         * @param other Object to compare with or transfer from.
         * @return True when the documented condition holds; otherwise false.
         */
        bool operator!=(const ChildrenIterator& other) const { return !(*this == other); }
    };

    /**
     * @brief Range wrapper for direct-child iteration.
     */
    class ChildrenRange {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Stores the first local. */
        NodeId firstLocal_ = InvalidNode;
        /** @brief Stores the expected version. */
        std::size_t expectedVersion_ = 0;

      public:
        /**
         * @brief Creates an empty child range.
         */
        ChildrenRange() = default;

        /**
         * @brief Creates a range over a linked list of direct children.
         *
         * @param tree Tree topology used by the operation.
         * @param firstLocal First local iterator position.
         * @param expectedVersion Mutation version required by the operation.
         */
        ChildrenRange(const MorphologicalTree* tree, NodeId firstLocal, std::size_t expectedVersion)
            : T_(tree), firstLocal_(firstLocal), expectedVersion_(expectedVersion) {}

        /**
         * @brief Returns an iterator at the first child.
         *
         * @return An iterator at the first child.
         */
        ChildrenIterator begin() const { return ChildrenIterator(T_, firstLocal_, expectedVersion_); }

        /**
         * @brief Returns the child-range sentinel iterator.
         *
         * @return The child-range sentinel iterator.
         */
        ChildrenIterator end() const { return ChildrenIterator(T_, InvalidNode, expectedVersion_); }
    };

    /**
     * @brief Iterator over the direct proper parts owned by one node.
     */
    class ProperPartsIterator {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Stores the current proper part. */
        NodeId currentProperPart_ = InvalidNode;
        /** @brief Stores the expected version. */
        std::size_t expectedVersion_ = 0;

      public:
        /// Standard iterator category exposed for STL compatibility.
        using iterator_category = std::forward_iterator_tag;
        /// Proper-part id yielded by the iterator.
        using value_type = NodeId;
        /// Signed distance type exposed for STL compatibility.
        using difference_type = std::ptrdiff_t;
        /// Pointer type exposed for STL compatibility.
        using pointer = const NodeId*;
        /// Reference type exposed for STL compatibility.
        using reference = const NodeId&;

        /**
         * @brief Creates the end/sentinel proper-part iterator.
         */
        ProperPartsIterator() = default;

        /**
         * @brief Creates a proper-part iterator starting from a linked-list entry.
         *
         * @param tree Tree topology used by the operation.
         * @param currentProperPart Proper-part data represented by `currentProperPart`.
         * @param expectedVersion Mutation version required by the operation.
         */
        ProperPartsIterator(const MorphologicalTree* tree, NodeId currentProperPart, std::size_t expectedVersion)
            : T_(tree), currentProperPart_(currentProperPart), expectedVersion_(expectedVersion) {}

        /**
         * @brief Advances to the next proper part in the owner list.
         *
         * @return Reference to the resulting object.
         */
        ProperPartsIterator& operator++() {
            T_->checkProperPartIteratorVersion(expectedVersion_);
            if (T_ && currentProperPart_ != InvalidNode) {
                currentProperPart_ = T_->nextProperPart_[currentProperPart_];
            }
            return *this;
        }

        /**
         * @brief Returns the current proper-part id.
         *
         * @return The current proper-part id.
         */
        NodeId operator*() const {
            T_->checkProperPartIteratorVersion(expectedVersion_);
            return currentProperPart_;
        }

        /**
         * @brief Compares proper-part iterator positions.
         *
         * @param other Object to compare with or transfer from.
         * @return True when the documented condition holds; otherwise false.
         */
        bool operator==(const ProperPartsIterator& other) const { return currentProperPart_ == other.currentProperPart_; }

        /**
         * @brief Compares proper-part iterator positions for inequality.
         *
         * @param other Object to compare with or transfer from.
         * @return True when the documented condition holds; otherwise false.
         */
        bool operator!=(const ProperPartsIterator& other) const { return !(*this == other); }
    };

    /**
     * @brief Range wrapper for direct proper-part iteration.
     */
    class ProperPartsRange {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Stores the first proper part. */
        NodeId firstProperPart_ = InvalidNode;
        /** @brief Stores the expected version. */
        std::size_t expectedVersion_ = 0;

      public:
        /**
         * @brief Creates an empty direct proper-part range.
         */
        ProperPartsRange() = default;

        /**
         * @brief Creates a range over one node's direct proper-part list.
         *
         * @param tree Tree topology used by the operation.
         * @param firstProperPart Proper-part data represented by `firstProperPart`.
         * @param expectedVersion Mutation version required by the operation.
         */
        ProperPartsRange(const MorphologicalTree* tree, NodeId firstProperPart, std::size_t expectedVersion)
            : T_(tree), firstProperPart_(firstProperPart), expectedVersion_(expectedVersion) {}

        /**
         * @brief Returns an iterator at the first direct proper part.
         *
         * @return An iterator at the first direct proper part.
         */
        ProperPartsIterator begin() const { return ProperPartsIterator(T_, firstProperPart_, expectedVersion_); }

        /**
         * @brief Returns the direct proper-part range sentinel.
         *
         * @return The direct proper-part range sentinel.
         */
        ProperPartsIterator end() const { return ProperPartsIterator(T_, InvalidNode, expectedVersion_); }
    };

    /**
     * @brief Iterator over all proper parts in one connected component.
     */
    class ConnectedComponentIterator {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Stores the node stack. */
        std::vector<NodeId> nodeStack_;
        /** @brief Stores the current proper part. */
        NodeId currentProperPart_ = InvalidNode;
        /** @brief Stores the expected topology version. */
        std::size_t expectedTopologyVersion_ = 0;
        /** @brief Stores the expected proper part version. */
        std::size_t expectedProperPartVersion_ = 0;

        /**
         * @brief Checks that topology and proper-part storage have not changed.
         */
        void checkVersions() const {
            T_->checkTopologyIteratorVersion(expectedTopologyVersion_);
            T_->checkProperPartIteratorVersion(expectedProperPartVersion_);
        }

        /**
         * @brief Pushes direct children in reverse order so traversal remains left-to-right.
         *
         * @param nodeId Identifier of the node used by the operation.
         */
        void pushChildren(NodeId nodeId) {
            std::vector<NodeId> children;
            for (NodeId childId : T_->getChildren(nodeId)) {
                children.push_back(childId);
            }
            for (auto it = children.rbegin(); it != children.rend(); ++it) {
                nodeStack_.push_back(*it);
            }
        }

        /**
         * @brief Advances to the next available proper part in DFS subtree order.
         */
        void settle() {
            checkVersions();
            while (currentProperPart_ == InvalidNode && !nodeStack_.empty()) {
                const NodeId nodeId = nodeStack_.back();
                nodeStack_.pop_back();
                pushChildren(nodeId);
                currentProperPart_ = T_->properHead_[static_cast<size_t>(nodeId)];
            }
        }

      public:
        /// Standard iterator category exposed for STL compatibility.
        using iterator_category = std::forward_iterator_tag;
        /// Proper-part id yielded by the iterator.
        using value_type = NodeId;
        /// Signed distance type exposed for STL compatibility.
        using difference_type = std::ptrdiff_t;
        /// Pointer type exposed for STL compatibility.
        using pointer = const NodeId*;
        /// Reference type exposed for STL compatibility.
        using reference = const NodeId&;

        /**
         * @brief Creates the end/sentinel connected-component iterator.
         */
        ConnectedComponentIterator() = default;

        /**
         * @brief Creates an iterator over every proper part in a rooted subtree.
         *
         * @param tree Tree topology used by the operation.
         * @param rootNodeId Identifier of the traversal root.
         * @param expectedTopologyVersion Topology version captured by the iterator.
         * @param expectedProperPartVersion Proper-part data represented by `expectedProperPartVersion`.
         */
        ConnectedComponentIterator(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedTopologyVersion, std::size_t expectedProperPartVersion)
            : T_(tree), expectedTopologyVersion_(expectedTopologyVersion), expectedProperPartVersion_(expectedProperPartVersion) {
            if (T_ && rootNodeId != InvalidNode) {
                nodeStack_.push_back(rootNodeId);
                settle();
            }
        }

        /**
         * @brief Advances to the next proper part in the connected component.
         *
         * @return Reference to the resulting object.
         */
        ConnectedComponentIterator& operator++() {
            checkVersions();
            if (currentProperPart_ != InvalidNode) {
                currentProperPart_ = T_->nextProperPart_[static_cast<size_t>(currentProperPart_)];
            }
            settle();
            return *this;
        }

        /**
         * @brief Returns the current proper-part id.
         *
         * @return The current proper-part id.
         */
        NodeId operator*() const {
            checkVersions();
            return currentProperPart_;
        }

        /**
         * @brief Compares connected-component iterator positions.
         *
         * @param other Object to compare with or transfer from.
         * @return True when the documented condition holds; otherwise false.
         */
        bool operator==(const ConnectedComponentIterator& other) const { return currentProperPart_ == other.currentProperPart_; }

        /**
         * @brief Compares connected-component iterator positions for inequality.
         *
         * @param other Object to compare with or transfer from.
         * @return True when the documented condition holds; otherwise false.
         */
        bool operator!=(const ConnectedComponentIterator& other) const { return !(*this == other); }
    };

    /**
     * @brief Range wrapper for connected-component proper-part iteration.
     */
    class ConnectedComponentRange {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Stores the root node identifier. */
        NodeId rootNodeId_ = InvalidNode;
        /** @brief Stores the expected topology version. */
        std::size_t expectedTopologyVersion_ = 0;
        /** @brief Stores the expected proper part version. */
        std::size_t expectedProperPartVersion_ = 0;

      public:
        /**
         * @brief Creates an empty connected-component range.
         */
        ConnectedComponentRange() = default;

        /**
         * @brief Creates a range over all proper parts in the subtree of `rootNodeId`.
         *
         * @param tree Tree topology used by the operation.
         * @param rootNodeId Identifier of the traversal root.
         * @param expectedTopologyVersion Topology version captured by the iterator.
         * @param expectedProperPartVersion Proper-part data represented by `expectedProperPartVersion`.
         */
        ConnectedComponentRange(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedTopologyVersion, std::size_t expectedProperPartVersion)
            : T_(tree), rootNodeId_(rootNodeId), expectedTopologyVersion_(expectedTopologyVersion), expectedProperPartVersion_(expectedProperPartVersion) {}

        /**
         * @brief Returns an iterator at the first proper part in the component.
         *
         * @return An iterator at the first proper part in the component.
         */
        ConnectedComponentIterator begin() const { return ConnectedComponentIterator(T_, rootNodeId_, expectedTopologyVersion_, expectedProperPartVersion_); }

        /**
         * @brief Returns the connected-component range sentinel.
         *
         * @return The connected-component range sentinel.
         */
        ConnectedComponentIterator end() const { return ConnectedComponentIterator(); }
    };

    /**
     * @brief Post-order iterator over one subtree.
     */
    class PostOrderNodeIterator {
      private:
        /** @brief Stores one traversal stack entry and its expansion state. */
        struct Item {
            /** @brief Stores the identifier. */
            NodeId id;
            /** @brief Indicates whether the node's children were already expanded. */
            bool expanded;
        };
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Stores the stack. */
        std::vector<Item> stack_;
        /** @brief Stores the current. */
        NodeId current_ = InvalidNode;
        /** @brief Stores the expected version. */
        std::size_t expectedVersion_ = 0;

        /**
         * @brief Expands the DFS stack until the next post-order node is available.
         */
        void settle_() {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            while (!stack_.empty()) {
                Item& top = stack_.back();
                if (!top.expanded) {
                    top.expanded = true;
                    std::vector<NodeId> children;
                    for (NodeId childId : T_->getChildren(top.id)) {
                        children.push_back(childId);
                    }
                    for (auto it = children.rbegin(); it != children.rend(); ++it) {
                        stack_.push_back({*it, false});
                    }
                } else {
                    current_ = top.id;
                    return;
                }
            }
            current_ = InvalidNode;
        }

      public:
        /// Standard iterator category exposed for STL compatibility.
        using iterator_category = std::input_iterator_tag;
        /// Node id yielded by the iterator.
        using value_type = NodeId;
        /// Signed distance type exposed for STL compatibility.
        using difference_type = std::ptrdiff_t;
        /// Pointer type exposed for STL compatibility.
        using pointer = const NodeId*;
        /// Reference type exposed for STL compatibility.
        using reference = const NodeId&;

        /**
         * @brief Creates the end/sentinel post-order iterator.
         */
        PostOrderNodeIterator() = default;

        /**
         * @brief Creates a post-order iterator rooted at `rootNodeId`.
         *
         * @param tree Tree topology used by the operation.
         * @param rootNodeId Identifier of the traversal root.
         * @param expectedVersion Mutation version required by the operation.
         */
        PostOrderNodeIterator(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedVersion) : T_(tree), expectedVersion_(expectedVersion) {
            if (T_ && rootNodeId != InvalidNode) {
                stack_.push_back({rootNodeId, false});
                settle_();
            }
        }

        /**
         * @brief Advances to the next node in post-order.
         *
         * @return Reference to the resulting object.
         */
        PostOrderNodeIterator& operator++() {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            if (!stack_.empty()) {
                stack_.pop_back();
                settle_();
            }
            return *this;
        }

        /**
         * @brief Returns the current post-order node id.
         *
         * @return The current post-order node id.
         */
        NodeId operator*() const {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            return current_;
        }

        /**
         * @brief Compares post-order iterator positions.
         *
         * @param other Object to compare with or transfer from.
         * @return True when the documented condition holds; otherwise false.
         */
        bool operator==(const PostOrderNodeIterator& other) const { return current_ == other.current_; }

        /**
         * @brief Compares post-order iterator positions for inequality.
         *
         * @param other Object to compare with or transfer from.
         * @return True when the documented condition holds; otherwise false.
         */
        bool operator!=(const PostOrderNodeIterator& other) const { return !(*this == other); }
    };

    /**
     * @brief Range wrapper for post-order subtree traversal.
     */
    class PostOrderNodeRange {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Stores the root node identifier. */
        NodeId rootNodeId_ = InvalidNode;
        /** @brief Stores the expected version. */
        std::size_t expectedVersion_ = 0;

      public:
        /**
         * @brief Creates an empty post-order range.
         */
        PostOrderNodeRange() = default;

        /**
         * @brief Creates a post-order range rooted at `rootNodeId`.
         *
         * @param tree Tree topology used by the operation.
         * @param rootNodeId Identifier of the traversal root.
         * @param expectedVersion Mutation version required by the operation.
         */
        PostOrderNodeRange(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedVersion)
            : T_(tree), rootNodeId_(rootNodeId), expectedVersion_(expectedVersion) {}

        /**
         * @brief Returns an iterator at the first post-order node.
         *
         * @return An iterator at the first post-order node.
         */
        PostOrderNodeIterator begin() const { return PostOrderNodeIterator(T_, rootNodeId_, expectedVersion_); }

        /**
         * @brief Returns the post-order range sentinel.
         *
         * @return The post-order range sentinel.
         */
        PostOrderNodeIterator end() const { return PostOrderNodeIterator(); }
    };

    /**
     * @brief Breadth-first iterator over one subtree.
     */
    class BreadthFirstNodeIterator {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Stores the queue. */
        FastQueue<NodeId> queue_;
        /** @brief Stores the expected version. */
        std::size_t expectedVersion_ = 0;

      public:
        /// Standard iterator category exposed for STL compatibility.
        using iterator_category = std::input_iterator_tag;
        /// Node id yielded by the iterator.
        using value_type = NodeId;
        /// Signed distance type exposed for STL compatibility.
        using difference_type = std::ptrdiff_t;
        /// Pointer type exposed for STL compatibility.
        using pointer = const NodeId*;
        /// Reference type exposed for STL compatibility.
        using reference = const NodeId&;

        /**
         * @brief Creates the end/sentinel breadth-first iterator.
         */
        BreadthFirstNodeIterator() = default;

        /**
         * @brief Creates a breadth-first iterator rooted at `rootNodeId`.
         *
         * @param tree Tree topology used by the operation.
         * @param rootNodeId Identifier of the traversal root.
         * @param expectedVersion Mutation version required by the operation.
         */
        BreadthFirstNodeIterator(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedVersion) : T_(tree), expectedVersion_(expectedVersion) {
            if (T_ && rootNodeId != InvalidNode) {
                queue_.push(rootNodeId);
            }
        }

        /**
         * @brief Advances to the next node in breadth-first order.
         *
         * @return Reference to the resulting object.
         */
        BreadthFirstNodeIterator& operator++() {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            if (!queue_.empty()) {
                NodeId current = queue_.pop();
                for (NodeId childId : T_->getChildren(current)) {
                    queue_.push(childId);
                }
            }
            return *this;
        }

        /**
         * @brief Returns the current breadth-first node id.
         *
         * @return The current breadth-first node id.
         */
        NodeId operator*() const {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            return queue_.front();
        }

        /**
         * @brief Compares breadth-first iterator exhaustion state.
         *
         * @param other Object to compare with or transfer from.
         * @return True when the documented condition holds; otherwise false.
         */
        bool operator==(const BreadthFirstNodeIterator& other) const { return queue_.empty() == other.queue_.empty(); }

        /**
         * @brief Compares breadth-first iterator exhaustion state for inequality.
         *
         * @param other Object to compare with or transfer from.
         * @return True when the documented condition holds; otherwise false.
         */
        bool operator!=(const BreadthFirstNodeIterator& other) const { return !(*this == other); }
    };

    /**
     * @brief Range wrapper for breadth-first subtree traversal.
     */
    class BreadthFirstNodeRange {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Stores the root node identifier. */
        NodeId rootNodeId_ = InvalidNode;
        /** @brief Stores the expected version. */
        std::size_t expectedVersion_ = 0;

      public:
        /**
         * @brief Creates an empty breadth-first range.
         */
        BreadthFirstNodeRange() = default;

        /**
         * @brief Creates a breadth-first range rooted at `rootNodeId`.
         *
         * @param tree Tree topology used by the operation.
         * @param rootNodeId Identifier of the traversal root.
         * @param expectedVersion Mutation version required by the operation.
         */
        BreadthFirstNodeRange(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedVersion)
            : T_(tree), rootNodeId_(rootNodeId), expectedVersion_(expectedVersion) {}

        /**
         * @brief Returns an iterator at the first breadth-first node.
         *
         * @return An iterator at the first breadth-first node.
         */
        BreadthFirstNodeIterator begin() const { return BreadthFirstNodeIterator(T_, rootNodeId_, expectedVersion_); }

        /**
         * @brief Returns the breadth-first range sentinel.
         *
         * @return The breadth-first range sentinel.
         */
        BreadthFirstNodeIterator end() const { return BreadthFirstNodeIterator(); }
    };

    /**
     * @brief Iterator that walks from a node towards the root.
     */
    class PathToRootIterator {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Stores the current. */
        NodeId current_ = InvalidNode;
        /** @brief Stores the expected version. */
        std::size_t expectedVersion_ = 0;

      public:
        /// Standard iterator category exposed for STL compatibility.
        using iterator_category = std::input_iterator_tag;
        /// Node id yielded by the iterator.
        using value_type = NodeId;
        /// Signed distance type exposed for STL compatibility.
        using difference_type = std::ptrdiff_t;
        /// Pointer type exposed for STL compatibility.
        using pointer = const NodeId*;
        /// Reference type exposed for STL compatibility.
        using reference = const NodeId&;

        /**
         * @brief Creates the end/sentinel rootward-path iterator.
         */
        PathToRootIterator() = default;

        /**
         * @brief Creates an iterator starting at one node and walking to the root.
         *
         * @param tree Tree topology used by the operation.
         * @param current Current iterator or traversal position.
         * @param expectedVersion Mutation version required by the operation.
         */
        PathToRootIterator(const MorphologicalTree* tree, NodeId current, std::size_t expectedVersion)
            : T_(tree), current_(current), expectedVersion_(expectedVersion) {}

        /**
         * @brief Advances one step toward the root.
         *
         * @return Reference to the resulting object.
         */
        PathToRootIterator& operator++() {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            if (T_ && current_ != InvalidNode) {
                if (T_->isRoot(current_)) {
                    current_ = InvalidNode;
                } else {
                    current_ = T_->getNodeParent(current_);
                }
            }
            return *this;
        }

        /**
         * @brief Returns the current node on the rootward path.
         *
         * @return The current node on the rootward path.
         */
        NodeId operator*() const {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            return current_;
        }

        /**
         * @brief Compares path-to-root iterator positions.
         *
         * @param other Object to compare with or transfer from.
         * @return True when the documented condition holds; otherwise false.
         */
        bool operator==(const PathToRootIterator& other) const { return current_ == other.current_; }

        /**
         * @brief Compares path-to-root iterator positions for inequality.
         *
         * @param other Object to compare with or transfer from.
         * @return True when the documented condition holds; otherwise false.
         */
        bool operator!=(const PathToRootIterator& other) const { return !(*this == other); }
    };

    /**
     * @brief Range wrapper for rootward path traversal.
     */
    class PathToRootRange {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Stores the start. */
        NodeId start_ = InvalidNode;
        /** @brief Stores the expected version. */
        std::size_t expectedVersion_ = 0;

      public:
        /**
         * @brief Creates an empty rootward-path range.
         */
        PathToRootRange() = default;

        /**
         * @brief Creates a range from `start` to the connected root.
         *
         * @param tree Tree topology used by the operation.
         * @param start Inclusive start position.
         * @param expectedVersion Mutation version required by the operation.
         */
        PathToRootRange(const MorphologicalTree* tree, NodeId start, std::size_t expectedVersion)
            : T_(tree), start_(start), expectedVersion_(expectedVersion) {}

        /**
         * @brief Returns an iterator at the path start node.
         *
         * @return An iterator at the path start node.
         */
        PathToRootIterator begin() const { return PathToRootIterator(T_, start_, expectedVersion_); }

        /**
         * @brief Returns the rootward-path range sentinel.
         *
         * @return The rootward-path range sentinel.
         */
        PathToRootIterator end() const { return PathToRootIterator(); }
    };

    /**
     * @brief Iterator over a materialised path between two nodes.
     */
    class PathBetweenNodesIterator {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Stores the path. */
        const std::vector<NodeId>* path_ = nullptr;
        /** @brief Stores the index. */
        std::size_t index_ = 0;
        /** @brief Stores the expected version. */
        std::size_t expectedVersion_ = 0;

      public:
        /// Standard iterator category exposed for STL compatibility.
        using iterator_category = std::input_iterator_tag;
        /// Node id yielded by the iterator.
        using value_type = NodeId;
        /// Signed distance type exposed for STL compatibility.
        using difference_type = std::ptrdiff_t;
        /// Pointer type exposed for STL compatibility.
        using pointer = const NodeId*;
        /// Reference type exposed for STL compatibility.
        using reference = const NodeId&;

        /**
         * @brief Creates the end/sentinel path-between-nodes iterator.
         */
        PathBetweenNodesIterator() = default;

        /**
         * @brief Creates an iterator over a materialised node path.
         *
         * @param tree Tree topology used by the operation.
         * @param path Node path updated by the edit.
         * @param index Zero-based index used by the operation.
         * @param expectedVersion Mutation version required by the operation.
         */
        PathBetweenNodesIterator(const MorphologicalTree* tree, const std::vector<NodeId>* path, std::size_t index, std::size_t expectedVersion)
            : T_(tree), path_(path), index_(index), expectedVersion_(expectedVersion) {}

        /**
         * @brief Advances to the next node in the materialised path.
         *
         * @return Reference to the resulting object.
         */
        PathBetweenNodesIterator& operator++() {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            if (path_ && index_ < path_->size()) {
                ++index_;
            }
            return *this;
        }

        /**
         * @brief Returns the current node id in the materialised path.
         *
         * @return The current node id in the materialised path.
         */
        NodeId operator*() const {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            return (*path_)[index_];
        }

        /**
         * @brief Compares path iterator positions.
         *
         * @param other Object to compare with or transfer from.
         * @return True when the documented condition holds; otherwise false.
         */
        bool operator==(const PathBetweenNodesIterator& other) const { return path_ == other.path_ && index_ == other.index_; }

        /**
         * @brief Compares path iterator positions for inequality.
         *
         * @param other Object to compare with or transfer from.
         * @return True when the documented condition holds; otherwise false.
         */
        bool operator!=(const PathBetweenNodesIterator& other) const { return !(*this == other); }
    };

    /**
     * @brief Range wrapper for the path connecting two nodes in the same component.
     */
    class PathBetweenNodesRange {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Stores the path. */
        std::vector<NodeId> path_;
        /** @brief Stores the expected version. */
        std::size_t expectedVersion_ = 0;

        /**
         * @brief Materialises the path joining two nodes in the same connected component.
         *
         * @param tree Tree topology used by the operation.
         * @param sourceNodeId Node identifier represented by `sourceNodeId`.
         * @param targetNodeId Node identifier represented by `targetNodeId`.
         * @return Values produced by the operation.
         */
        static std::vector<NodeId> buildPath(const MorphologicalTree* tree, NodeId sourceNodeId, NodeId targetNodeId) {
            if (tree == nullptr || !tree->isAlive(sourceNodeId) || !tree->isAlive(targetNodeId)) {
                return {};
            }

            auto componentAnchor = [tree](NodeId nodeId) {
                NodeId currentNodeId = nodeId;
                while (true) {
                    const NodeId parentNodeId = tree->getNodeParent(currentNodeId);
                    if (parentNodeId == InvalidNode || parentNodeId == currentNodeId) {
                        return currentNodeId;
                    }
                    currentNodeId = parentNodeId;
                }
            };

            if (componentAnchor(sourceNodeId) != componentAnchor(targetNodeId)) {
                return {};
            }

            const NodeId lcaNodeId = tree->getLowestCommonAncestor(sourceNodeId, targetNodeId);
            if (lcaNodeId == InvalidNode) {
                return {};
            }

            std::vector<NodeId> path;
            for (NodeId currentNodeId = sourceNodeId;; currentNodeId = tree->getNodeParent(currentNodeId)) {
                path.push_back(currentNodeId);
                if (currentNodeId == lcaNodeId) {
                    break;
                }

                const NodeId parentNodeId = tree->getNodeParent(currentNodeId);
                if (parentNodeId == InvalidNode || parentNodeId == currentNodeId) {
                    return {};
                }
            }

            std::vector<NodeId> descendingTail;
            for (NodeId currentNodeId = targetNodeId; currentNodeId != lcaNodeId; currentNodeId = tree->getNodeParent(currentNodeId)) {
                descendingTail.push_back(currentNodeId);

                const NodeId parentNodeId = tree->getNodeParent(currentNodeId);
                if (parentNodeId == InvalidNode || parentNodeId == currentNodeId) {
                    return {};
                }
            }

            path.insert(path.end(), descendingTail.rbegin(), descendingTail.rend());
            return path;
        }

      public:
        /**
         * @brief Creates an empty path-between-nodes range.
         */
        PathBetweenNodesRange() = default;

        /**
         * @brief Materialises and owns the path between two nodes.
         *
         * @param tree Tree topology used by the operation.
         * @param sourceNodeId Node identifier represented by `sourceNodeId`.
         * @param targetNodeId Node identifier represented by `targetNodeId`.
         * @param expectedVersion Mutation version required by the operation.
         */
        PathBetweenNodesRange(const MorphologicalTree* tree, NodeId sourceNodeId, NodeId targetNodeId, std::size_t expectedVersion)
            : T_(tree), path_(buildPath(tree, sourceNodeId, targetNodeId)), expectedVersion_(expectedVersion) {}

        /**
         * @brief Returns an iterator at the first node in the materialised path.
         *
         * @return An iterator at the first node in the materialised path.
         */
        PathBetweenNodesIterator begin() const { return PathBetweenNodesIterator(T_, &path_, 0, expectedVersion_); }

        /**
         * @brief Returns the materialised path range sentinel.
         *
         * @return The materialised path range sentinel.
         */
        PathBetweenNodesIterator end() const { return PathBetweenNodesIterator(T_, &path_, path_.size(), expectedVersion_); }
    };

    /**
     * @brief Depth-first iterator over a subtree in pre-order.
     */
    class SubtreeNodeIterator {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Stores the stack. */
        std::vector<NodeId> stack_;
        /** @brief Stores the expected version. */
        std::size_t expectedVersion_ = 0;

      public:
        /// Standard iterator category exposed for STL compatibility.
        using iterator_category = std::input_iterator_tag;
        /// Node id yielded by the iterator.
        using value_type = NodeId;
        /// Signed distance type exposed for STL compatibility.
        using difference_type = std::ptrdiff_t;
        /// Pointer type exposed for STL compatibility.
        using pointer = const NodeId*;
        /// Reference type exposed for STL compatibility.
        using reference = const NodeId&;

        /**
         * @brief Creates the end/sentinel subtree iterator.
         */
        SubtreeNodeIterator() = default;

        /**
         * @brief Creates a pre-order subtree iterator rooted at `rootNodeId`.
         *
         * @param tree Tree topology used by the operation.
         * @param rootNodeId Identifier of the traversal root.
         * @param expectedVersion Mutation version required by the operation.
         */
        SubtreeNodeIterator(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedVersion) : T_(tree), expectedVersion_(expectedVersion) {
            if (T_ && rootNodeId != InvalidNode) {
                stack_.push_back(rootNodeId);
            }
        }

        /**
         * @brief Advances to the next node in pre-order subtree traversal.
         *
         * @return Reference to the resulting object.
         */
        SubtreeNodeIterator& operator++() {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            if (!stack_.empty()) {
                NodeId current = stack_.back();
                stack_.pop_back();
                std::vector<NodeId> children;
                for (NodeId childId : T_->getChildren(current)) {
                    children.push_back(childId);
                }
                for (auto it = children.rbegin(); it != children.rend(); ++it) {
                    stack_.push_back(*it);
                }
            }
            return *this;
        }

        /**
         * @brief Returns the current pre-order subtree node id.
         *
         * @return The current pre-order subtree node id.
         */
        NodeId operator*() const {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            return stack_.back();
        }

        /**
         * @brief Compares subtree iterator exhaustion state.
         *
         * @param other Object to compare with or transfer from.
         * @return True when the documented condition holds; otherwise false.
         */
        bool operator==(const SubtreeNodeIterator& other) const { return stack_.empty() == other.stack_.empty(); }

        /**
         * @brief Compares subtree iterator exhaustion state for inequality.
         *
         * @param other Object to compare with or transfer from.
         * @return True when the documented condition holds; otherwise false.
         */
        bool operator!=(const SubtreeNodeIterator& other) const { return !(*this == other); }
    };

    /**
     * @brief Range wrapper for pre-order subtree traversal.
     */
    class SubtreeNodeRange {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Stores the root node identifier. */
        NodeId rootNodeId_ = InvalidNode;
        /** @brief Stores the expected version. */
        std::size_t expectedVersion_ = 0;

      public:
        /**
         * @brief Creates an empty subtree range.
         */
        SubtreeNodeRange() = default;

        /**
         * @brief Creates a pre-order range over the subtree rooted at `rootNodeId`.
         *
         * @param tree Tree topology used by the operation.
         * @param rootNodeId Identifier of the traversal root.
         * @param expectedVersion Mutation version required by the operation.
         */
        SubtreeNodeRange(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedVersion)
            : T_(tree), rootNodeId_(rootNodeId), expectedVersion_(expectedVersion) {}

        /**
         * @brief Returns an iterator at the subtree root.
         *
         * @return An iterator at the subtree root.
         */
        SubtreeNodeIterator begin() const { return SubtreeNodeIterator(T_, rootNodeId_, expectedVersion_); }

        /**
         * @brief Returns the subtree range sentinel.
         *
         * @return The subtree range sentinel.
         */
        SubtreeNodeIterator end() const { return SubtreeNodeIterator(); }
    };

    /**
     * @brief Range wrapper over the proper descendants of one node.
     */
    class DescendantNodeRange {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Stores the root node identifier. */
        NodeId rootNodeId_ = InvalidNode;
        /** @brief Stores the expected version. */
        std::size_t expectedVersion_ = 0;

      public:
        /**
         * @brief Creates an empty descendant range.
         */
        DescendantNodeRange() = default;

        /**
         * @brief Creates a range over proper descendants of `rootNodeId`.
         *
         * @param tree Tree topology used by the operation.
         * @param rootNodeId Identifier of the traversal root.
         * @param expectedVersion Mutation version required by the operation.
         */
        DescendantNodeRange(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedVersion)
            : T_(tree), rootNodeId_(rootNodeId), expectedVersion_(expectedVersion) {}

        /**
         * @brief Returns an iterator at the first proper descendant.
         *
         * @return An iterator at the first proper descendant.
         */
        SubtreeNodeIterator begin() const {
            auto it = SubtreeNodeIterator(T_, rootNodeId_, expectedVersion_);
            ++it;
            return it;
        }

        /**
         * @brief Returns the descendant range sentinel.
         *
         * @return The descendant range sentinel.
         */
        SubtreeNodeIterator end() const { return SubtreeNodeIterator(); }
    };
};

} // namespace mmcfilters
