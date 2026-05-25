#pragma once

#include "../utils/AdjacencyRelation.hpp"
#include "../utils/Assert.hpp"
#include "../utils/Altitude.hpp"
#include "../utils/Image.hpp"
#include "../utils/Common.hpp"
#include "../dataStructure/FastQueue.hpp"
#include "BuilderMorphologicalTreeByUnionFind.hpp"
#include "detail/MorphologicalTreeConstructionTag.hpp"
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
enum class NodeIdSpace {
    MORPHOLOGICAL_TREE,
    HIGRA
};

/**
 * @brief Component-tree polarity used by image-based max/min construction.
 *
 * This enum is intentionally limited to component-tree builders. Use
 * `MorphologicalTreeKind` for generic factory APIs that may also import trees
 * of shapes.
 */
enum class ComponentTreeKind {
    MAX_TREE,
    MIN_TREE
};

/**
 * @brief Semantic kind of a morphological hierarchy accepted by the factory.
 *
 * Construction APIs use this strongly typed value instead of public integer
 * tree-type parameters.
 */
enum class MorphologicalTreeKind {
    MAX_TREE,
    MIN_TREE,
    TREE_OF_SHAPES,
    SELF_DUAL_RESIDUAL_TREE
};

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
    explicit operator bool() const noexcept {
        return ok;
    }
};

// Forward declaration for the edit-session wrapper.
class TreeEditor;

/**
 * @brief Mutable morphological tree built directly on proper parts and dense node ids.
 *
 * `MorphologicalTree` is the central mutable hierarchy of this project. It can
 * represent a max-tree, min-tree, or tree of shapes and exposes a dense
 * `NodeId` domain together with direct proper-part ownership. The class keeps
 * explicit parent/child links, linked lists of direct proper parts, optional
 * 2D image-domain metadata, and a small set of structural caches used by the
 * public traversal and ancestry queries.
 *
 * Data model:
 *
 * - proper parts are indexed by `NodeId` in the range `[0, getNumTotalProperParts())`;
 * - internal nodes are indexed by `NodeId` in the range `[0, getNumInternalNodeSlots())`;
 * - each live node owns zero or more direct proper parts and may have direct children;
 * - the full support of a node is the union of the direct proper parts in its subtree;
 * - the root and detached nodes point to themselves as parent.
 *
 * Main responsibilities:
 *
 * - construct the hierarchy from an image or a static Higra parent representation;
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

    class LCAEulerRMQ; // Forward declaration for the LCA cache implementation.

    // ========================= Private attributes ========================= //
    NodeId rootNodeId_ = InvalidNode;
    MorphologicalTreeKind treeType_ = MorphologicalTreeKind::MAX_TREE;
    int numRows_ = 0;
    int numCols_ = 0;
    std::optional<AdjacencyRelation> adj_; // optional adjacency context; not part of the structural tree contract
    std::optional<std::pair<AdjacencyRelation, AdjacencyRelation>> tosAdjacencyPolicy_; // min-tree adjacency, max-tree adjacency
    int numNodes_ = 0;
    bool preservesHigraNodeIdSpace_ = false;
    bool editSessionOpen_ = false;

    // Proper-part ownership, indexed by proper-part global id [0, getNumTotalProperParts()).
    std::vector<NodeId> properPartOwner_;

    // Parent links, indexed by local node-slot id [0, getNumInternalNodeSlots()).
    std::vector<NodeId> nodeParent_;

    // Internal hierarchy linked structure, indexed by local node-slot id [0, getNumInternalNodeSlots()).
    std::vector<NodeId> firstChild_;
    std::vector<NodeId> nextSibling_;
    std::vector<NodeId> prevSibling_;
    std::vector<NodeId> lastChild_;
    std::vector<int> numChildrenByNode_;

    // Free slot management and alive-node iteration, indexed by local node-slot id [0, getNumInternalNodeSlots()).
    std::vector<uint8_t> alive_;
    std::vector<NodeId> freeNodeIds_;

    // Proper-parts linked lists, indexed by local node-slot id [0, getNumTotalProperParts()).
    std::vector<NodeId> properHead_;
    std::vector<NodeId> properTail_;
    std::vector<int> numProperPartsByNode_;
    std::vector<NodeId> nextProperPart_;
    std::vector<NodeId> prevProperPart_;

    // Structural caches for traversal-based queries, indexed by local node-slot id [0, getNumInternalNodeSlots()).
    struct PrePostOrderCache {
        std::vector<int> timePreOrder;
        std::vector<int> timePostOrder;
        bool valid = false;

        /**
         * @brief Marks the cached traversal timestamps as stale.
         */
        void invalidate() noexcept { valid = false; }
    };
    mutable PrePostOrderCache prePostOrderCache_;
    mutable std::unique_ptr<LCAEulerRMQ> lcaCache_;


    // Version counters for iterator invalidation.
    std::size_t nodeStructureVersion_ = 0;
    std::size_t topologyVersion_ = 0;
    std::size_t properPartVersion_ = 0;
    std::size_t mutationVersion_ = 0;

    enum class ChildSplicePolicy {
        AppendToTargetTail,
        ReplaceSourceSlotWhenDirectChild
    };

    // ========================= Private methods ========================= //
    /**
     * @brief Allocates one live node slot and optionally attaches it to a parent.
     */
    inline NodeId makeNode(NodeId parentSlotId) {
        const NodeId id = this->allocateSlot();

        if (parentSlotId >= 0) {
            linkChildSlot(parentSlotId, id);
        } else {
            nodeParent_[id] = id;
        }
        this->numNodes_++;
        return id;
    }

    /**
     * @brief Returns a reusable slot or appends a fresh one to the dense storage.
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
     * @brief Reserves node-indexed storage ahead of a builder-driven construction.
     */
    inline void reserveNodes(int expected) {
        nodeParent_.reserve(expected);
        firstChild_.reserve(expected);
        nextSibling_.reserve(expected);
        prevSibling_.reserve(expected);
        lastChild_.reserve(expected);
        numChildrenByNode_.reserve(expected);
        alive_.reserve(expected);
        properHead_.reserve(expected);
        properTail_.reserve(expected);
        numProperPartsByNode_.reserve(expected);
    }

    /**
     * @brief Resets the proper-part linked-list storage to the requested size.
     */
    inline void initializeProperPartStorage(size_t numProperParts) {
        nextProperPart_.assign(numProperParts, InvalidNode);
        prevProperPart_.assign(numProperParts, InvalidNode);
    }

    /**
     * @brief Drops the preserved imported Higra id-domain marker after topology changes.
     */
    inline void invalidateHigraNodeIdSpace() noexcept { preservesHigraNodeIdSpace_ = false; }

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
    inline void endEditSession() noexcept {
        editSessionOpen_ = false;
    }

    /**
     * @brief Clears the current topology and keeps only an empty proper-part domain.
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
     */
    inline void releaseSlotStorage(NodeId slotId) noexcept {
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
        freeNodeIds_.push_back(slotId);
    }

    /**
     * @brief Tests whether a dense slot currently belongs to the free list.
     */
    inline bool isFreeSlot(NodeId slotId) const noexcept {
        return slotId >= 0 && slotId < static_cast<NodeId>(alive_.size()) && alive_[slotId] == 0;
    }

    /**
     * @brief Rejects null or empty image domains before builder code indexes them.
     */
    template<class PixelType>
    static void requireNonEmptyImageDomain(const ImagePtr<PixelType>& img, const char* context) {
        if (!img) {
            throw std::invalid_argument(std::string(context) + " requires a non-null image.");
        }
        if (img->getNumRows() <= 0 || img->getNumCols() <= 0 || img->getSize() <= 0) {
            throw std::invalid_argument(std::string(context) + " requires a non-empty 2D image.");
        }
    }

    /**
     * @brief Returns the auxiliary min/max component-tree adjacencies implied by a ToS interpolation policy.
     */
    static std::pair<AdjacencyRelation, AdjacencyRelation> treeOfShapesAdjacencyPolicy(ToSInterpolation interpolation, int rows, int cols) {
        switch (interpolation) {
            case ToSInterpolation::SelfDual:
                return {AdjacencyRelation(rows, cols, 1.0), AdjacencyRelation(rows, cols, 1.0)};
            case ToSInterpolation::Min4cMax8c:
                return {AdjacencyRelation(rows, cols, 1.0), AdjacencyRelation(rows, cols, 1.5)};
            case ToSInterpolation::Min8cMax4c:
                return {AdjacencyRelation(rows, cols, 1.5), AdjacencyRelation(rows, cols, 1.0)};
        }
        throw std::invalid_argument("Unsupported tree-of-shapes interpolation.");
    }

    /**
     * @brief Rejects invalid or released internal-node ids before indexed reads.
     */
    inline void requireAliveNode(NodeId nodeId, const char* context) const {
        if (!isAlive(nodeId)) {
            throw std::invalid_argument(std::string(context) + " requires a live internal NodeId.");
        }
    }

    /**
     * @brief Rejects invalid, released, or root node ids for non-root-only edits.
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

        for (NodeId pixelId = 0; pixelId < static_cast<NodeId>(properPartOwner_.size()); ++pixelId) {
            const NodeId ownerSlotId = properPartOwner_[pixelId];
            if (ownerSlotId == InvalidNode || ownerSlotId >= static_cast<NodeId>(nodeParent_.size()) || isFreeSlot(ownerSlotId)) {
                continue;
            }

            const NodeId tailPixelId = properTail_[ownerSlotId];
            if (tailPixelId == InvalidNode) {
                properHead_[ownerSlotId] = pixelId;
                properTail_[ownerSlotId] = pixelId;
            } else {
                nextProperPart_[tailPixelId] = pixelId;
                prevProperPart_[pixelId] = tailPixelId;
                properTail_[ownerSlotId] = pixelId;
            }
            numProperPartsByNode_[ownerSlotId]++;
        }
    }

    /**
     * @brief Removes one direct proper part from the linked list of its current owner.
     */
    inline void unlinkProperPartFromOwner(NodeId ownerSlotId, NodeId pixelId) noexcept {
        const NodeId prev = prevProperPart_[static_cast<size_t>(pixelId)];
        const NodeId next = nextProperPart_[static_cast<size_t>(pixelId)];

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

        nextProperPart_[static_cast<size_t>(pixelId)] = InvalidNode;
        prevProperPart_[static_cast<size_t>(pixelId)] = InvalidNode;
        --numProperPartsByNode_[static_cast<size_t>(ownerSlotId)];
    }

    /**
     * @brief Appends one detached proper part to the target owner's direct list.
     */
    inline void appendDetachedProperPart(NodeId targetSlotId, NodeId pixelId) noexcept {
        properPartOwner_[static_cast<size_t>(pixelId)] = targetSlotId;
        const NodeId tail = properTail_[static_cast<size_t>(targetSlotId)];
        if (tail == InvalidNode) {
            properHead_[static_cast<size_t>(targetSlotId)] = pixelId;
            properTail_[static_cast<size_t>(targetSlotId)] = pixelId;
        } else {
            nextProperPart_[static_cast<size_t>(tail)] = pixelId;
            prevProperPart_[static_cast<size_t>(pixelId)] = tail;
            properTail_[static_cast<size_t>(targetSlotId)] = pixelId;
        }
        ++numProperPartsByNode_[static_cast<size_t>(targetSlotId)];
    }

    /**
     * @brief Splices all direct proper parts from `sourceSlotId` to the tail of `targetSlotId`.
     */
    inline void spliceProperPartsSlots(NodeId targetSlotId, NodeId sourceSlotId) noexcept {
        const NodeId sourceHead = properHead_[static_cast<size_t>(sourceSlotId)];
        if (sourceHead == InvalidNode) {
            return;
        }

        for (NodeId pixelId = sourceHead; pixelId != InvalidNode; pixelId = nextProperPart_[static_cast<size_t>(pixelId)]) {
            properPartOwner_[static_cast<size_t>(pixelId)] = targetSlotId;
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
     */
    inline void releaseSlotNode(NodeId slotNodeId) noexcept {
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
     */
    inline void checkNodeIteratorVersion([[maybe_unused]] std::size_t expectedVersion) const {
        assert(expectedVersion == nodeStructureVersion_ && "Alive-node iterator invalidated by node-structure mutation.");
    }

    /**
     * @brief Debug-only fail-fast check for topology iterators.
     */
    inline void checkTopologyIteratorVersion([[maybe_unused]] std::size_t expectedVersion) const {
        assert(expectedVersion == topologyVersion_ && "Topology iterator invalidated by tree-structure mutation.");
    }

    /**
     * @brief Debug-only fail-fast check for direct proper-part iterators.
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
            const_cast<MorphologicalTree*>(this),
            getRoot(),
            [&](NodeId nodeId) -> void {
                prePostOrderCache_.timePreOrder[nodeId] = timer++;
            },
            [&](NodeId, NodeId) -> void {},
            [&](NodeId nodeId) -> void {
                prePostOrderCache_.timePostOrder[nodeId] = timer++;
            }
        );

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
     */
    inline const LCAEulerRMQ& ensureLcaCache() const {
        if (!lcaCache_) {
            lcaCache_ = std::make_unique<LCAEulerRMQ>(this);
        }
        return *lcaCache_;
    }

    /**
     * @brief Generic depth-first accumulation helper used by incremental computations.
     */
    template<class PreProcessing, class MergeProcessing, class PostProcessing>
    static void computeIncrementalAttributes(MorphologicalTree* tree, NodeId rootNodeId, PreProcessing&& preProcessing, MergeProcessing&& mergeProcessing, PostProcessing&& postProcessing) {
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
     */
    void linkChildSlot(NodeId parentSlotId, NodeId childId) {
        if (parentSlotId < 0 || childId < 0) return;

        // If the child already has a parent, detach it before reading sibling links.
        if (nodeParent_[childId] != InvalidNode) {
            unlinkChildSlot(nodeParent_[childId], childId, false);
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
     * @brief Removes one child from its parent list and optionally releases the slot.
     */
    inline void unlinkChildSlot(NodeId parentSlotId, NodeId childId, bool release) {
        if (parentSlotId < 0 || childId < 0) return;
        if (nodeParent_[childId] != parentSlotId) return;

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
        if(release){
            releaseSlotNode(childId);
        } else {
            bumpTopologyVersion();
        }
    }


    /**
     * @brief Moves all children of `fromId` into the child list of `toId`.
     */
    inline void spliceChildrenSlots(NodeId toId, NodeId fromId, ChildSplicePolicy policy = ChildSplicePolicy::AppendToTargetTail) {
        if (toId < 0 || fromId < 0 || toId == fromId) return;

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
            if (firstFrom == InvalidNode) return;

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
     */
    inline void releaseNode(NodeId nodeId) noexcept {
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
     */
    inline void removeChild(NodeId parentNodeId, NodeId childId, bool releaseNodeFlag) {
        if (!isAlive(parentNodeId) || !isAlive(childId) || !hasChild(parentNodeId, childId)) {
            return;
        }
        const NodeId parentSlotId = parentNodeId;
        const NodeId childSlotId = childId;
        unlinkChildSlot(parentSlotId, childSlotId, false);
        nodeParent_[childSlotId] = childSlotId;
        invalidatePrePostOrderCache();
        if (releaseNodeFlag) {
            releaseNode(childId);
        }
    }

    /**
     * @brief Attaches `nodeId` as the last child of `parentNodeId`.
     */
    inline void attachNode(NodeId parentNodeId, NodeId nodeId) {
        if (!isAlive(parentNodeId) || !isAlive(nodeId) || isRoot(nodeId) || parentNodeId == nodeId) {
            return;
        }
        const NodeId parentSlotId = parentNodeId;
        const NodeId nodeSlotId = nodeId;
        const NodeId oldParentSlotId = nodeParent_[nodeSlotId];
        if (oldParentSlotId != InvalidNode && oldParentSlotId != nodeSlotId) {
            unlinkChildSlot(oldParentSlotId, nodeSlotId, false);
        }
        nodeParent_[nodeSlotId] = InvalidNode;
        linkChildSlot(parentSlotId, nodeSlotId);
    }

    /**
     * @brief Detaches `nodeId` from its current parent, leaving it self-parented.
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
        unlinkChildSlot(parentSlotId, nodeSlotId, false);
        nodeParent_[nodeSlotId] = nodeSlotId;
    }

    /**
     * @brief Moves `nodeId` under `newParentId`.
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
            unlinkChildSlot(oldParentSlotId, nodeSlotId, false);
        }
        nodeParent_[nodeSlotId] = InvalidNode;
        linkChildSlot(newParentSlotId, nodeSlotId);
    }

    /**
     * @brief Transfers every direct child of `sourceId` to `parentNodeId`.
     */
    inline void moveChildren(NodeId parentNodeId, NodeId sourceId) {
        if (!isAlive(parentNodeId) || !isAlive(sourceId) || parentNodeId == sourceId) {
            return;
        }
        spliceChildrenSlots(parentNodeId, sourceId);
    }

    /**
     * @brief Transfers one direct proper part from `sourceNodeId` to `targetNodeId`.
     */
    inline void moveProperPart(NodeId targetNodeId, NodeId sourceNodeId, NodeId pixelId) {
        if (!isAlive(targetNodeId) || !isAlive(sourceNodeId) || !isProperPart(pixelId) || targetNodeId == sourceNodeId) {
            return;
        }
        const NodeId sourceSlotId = sourceNodeId;
        if (properPartOwner_[pixelId] != sourceSlotId) {
            return;
        }
        unlinkProperPartFromOwner(sourceSlotId, pixelId);
        appendDetachedProperPart(targetNodeId, pixelId);
        bumpProperPartVersion();
    }

    /**
     * @brief Transfers all direct proper parts from `sourceNodeId` to `targetNodeId`.
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
            unlinkChildSlot(oldParentSlot, nodeSlot, false);
        }
        rootNodeId_ = nodeSlot;
        nodeParent_[nodeSlot] = nodeSlot;
        prevSibling_[nodeSlot] = InvalidNode;
        nextSibling_[nodeSlot] = InvalidNode;
        bumpTopologyVersion();
    }

    /**
     * @brief Private default constructor used by `clone()`.
     */
    MorphologicalTree() = default;


public:
    /**
     * @brief Copying is disabled to keep topology ownership explicit.
     */
    MorphologicalTree(const MorphologicalTree&) = delete;

    /**
     * @brief Copy assignment is disabled to keep topology ownership explicit.
     */
    MorphologicalTree& operator=(const MorphologicalTree&) = delete;

    /**
     * @brief Moving transfers the complete topology state.
     */
    MorphologicalTree(MorphologicalTree&&) noexcept = default;

    /**
     * @brief Move assignment transfers the complete topology state.
     */
    MorphologicalTree& operator=(MorphologicalTree&&) noexcept = default;

    /**
     * @brief Destroys the topology storage and cached traversal state.
     */
    virtual ~MorphologicalTree() = default;

    /**
     * @brief Tag-protected constructor for image-based max/min component trees.
     *
     * The factory owns this path. The resulting topology uses the image pixels
     * as proper parts, stores the selected component-tree type internally, and
     * keeps the adjacency relation used by later topology-aware operations.
     */
    template<AltitudeValue PixelType>
    MorphologicalTree(detail::MorphologicalTreeConstructionTag, ImagePtr<PixelType> imgPtr, ComponentTreeKind kind, double radius) {
        requireNonEmptyImageDomain(imgPtr, "MorphologicalTree component-tree construction");

        const bool isMaxtree = kind == ComponentTreeKind::MAX_TREE;
        treeType_ = isMaxtree ? MorphologicalTreeKind::MAX_TREE : MorphologicalTreeKind::MIN_TREE;
        numRows_ = imgPtr->getNumRows();
        numCols_ = imgPtr->getNumCols();
        adj_.emplace(numRows_, numCols_, radius);
        tosAdjacencyPolicy_ = std::nullopt;
        initializeEmptyStorage(static_cast<size_t>(imgPtr->getSize()));

        BuilderComponentTree builderUF(&*adj_, isMaxtree);
        auto [parent, orderedPixels, numBuiltNodes] = builderUF.createTreeByUnionFind<PixelType>(imgPtr);

        const int numPixels = imgPtr->getSize();
        auto img = imgPtr->rawData();
        this->reserveNodes(numBuiltNodes);
        for (int i = 0; i < numPixels; i++) {
            const int p = orderedPixels[i];

            if (p == parent[p]) {
                properPartOwner_[p] = this->rootNodeId_ = this->makeNode(InvalidNode);
            } else if (img[p] != img[parent[p]]) {
                properPartOwner_[p] = this->makeNode(properPartOwner_[parent[p]]);
            } else {
                properPartOwner_[p] = properPartOwner_[parent[p]];
            }
        }

        properHead_.assign(nodeParent_.size(), InvalidNode);
        properTail_.assign(nodeParent_.size(), InvalidNode);
        numProperPartsByNode_.assign(nodeParent_.size(), 0);
        initializeProperPartStorage(numPixels);
        rebuildProperPartLinksFromOwnership();
        invalidateAllIterators();

    }

    /**
     * @brief Tag-protected constructor for image-based tree-of-shapes topology.
     *
     * The factory owns this path. The topology uses the image pixels as proper
     * parts and stores the ToS interpolation policy as adjacency metadata for
     * operations that need to preserve the same interpretation.
     */
    MorphologicalTree(detail::MorphologicalTreeConstructionTag, ImageUInt8Ptr imgPtr, ToSInterpolation interpolation, int infinitySeedRow, int infinitySeedCol) {
        requireNonEmptyImageDomain(imgPtr, "MorphologicalTree tree-of-shapes construction");

        treeType_ = MorphologicalTreeKind::TREE_OF_SHAPES;
        numRows_ = imgPtr->getNumRows();
        numCols_ = imgPtr->getNumCols();
        adj_ = std::nullopt;
        tosAdjacencyPolicy_ = treeOfShapesAdjacencyPolicy(interpolation, numRows_, numCols_);
        initializeEmptyStorage(static_cast<size_t>(imgPtr->getSize()));

        BuilderTreeOfShape builderUF(interpolation, infinitySeedRow, infinitySeedCol);
        auto [parent, orderedPixels, numBuiltNodes] = builderUF.createTreeByUnionFind(imgPtr);

        const int numPixels = imgPtr->getSize();
        auto img = imgPtr->rawData();

        this->reserveNodes(numBuiltNodes);
        for (int i = 0; i < numPixels; i++) {
            const int p = orderedPixels[i];

            if (p == parent[p]) {
                properPartOwner_[p] = this->rootNodeId_ = this->makeNode(InvalidNode);
            } else if (img[p] != img[parent[p]]) {
                properPartOwner_[p] = this->makeNode(properPartOwner_[parent[p]]);
            } else {
                properPartOwner_[p] = properPartOwner_[parent[p]];
            }
        }

        properHead_.assign(nodeParent_.size(), InvalidNode);
        properTail_.assign(nodeParent_.size(), InvalidNode);
        numProperPartsByNode_.assign(nodeParent_.size(), 0);
        initializeProperPartStorage(numPixels);
        rebuildProperPartLinksFromOwnership();
        invalidateAllIterators();
    }

    /**
     * @brief Tag-protected import from compact Higra parent representation.
     *
     * Higra ids are leaves first and internal nodes after the image-domain
     * leaves. The constructor converts internal Higra ids to dense internal
     * `NodeId` slots while preserving enough mapping metadata for weighted
     * altitude import. Max/min imports require an adjacency relation; tree of
     * shapes imports can omit it.
     */
    MorphologicalTree(detail::MorphologicalTreeConstructionTag, std::span<const NodeId> parent, int rows, int cols, MorphologicalTreeKind kind, std::optional<AdjacencyRelation> adjacency) {
        if (kind == MorphologicalTreeKind::SELF_DUAL_RESIDUAL_TREE) {
            throw std::invalid_argument("Higra import does not accept SELF_DUAL_RESIDUAL_TREE.");
        }
        if (!adjacency && kind != MorphologicalTreeKind::TREE_OF_SHAPES) {
            throw std::invalid_argument("Higra import of max/min trees requires an explicit adjacency relation.");
        }

        treeType_ = kind;
        numRows_ = rows;
        numCols_ = cols;
        adj_ = std::move(adjacency);
        tosAdjacencyPolicy_ = std::nullopt;

        const NodeId numProperParts = static_cast<NodeId>(rows * cols);
        if (adj_) {
            if (adj_->getNumRows() * adj_->getNumCols() != numProperParts) {
                throw std::invalid_argument("Higra leaf count must match image domain size.");
            }
        }

        const NodeId numHigraNodes = static_cast<NodeId>(parent.size());
        if (numProperParts <= 0 || numProperParts >= numHigraNodes) {
            throw std::invalid_argument("Higra parent array must contain leaves followed by at least one internal node.");
        }
        const NodeId numNodeSlots = numHigraNodes - numProperParts;

        initializeEmptyStorage(static_cast<size_t>(numProperParts));
        nodeParent_.assign(static_cast<size_t>(numNodeSlots), InvalidNode);
        firstChild_.assign(static_cast<size_t>(numNodeSlots), InvalidNode);
        nextSibling_.assign(static_cast<size_t>(numNodeSlots), InvalidNode);
        prevSibling_.assign(static_cast<size_t>(numNodeSlots), InvalidNode);
        lastChild_.assign(static_cast<size_t>(numNodeSlots), InvalidNode);
        numChildrenByNode_.assign(static_cast<size_t>(numNodeSlots), 0);
        alive_.assign(static_cast<size_t>(numNodeSlots), 1);
        properHead_.assign(static_cast<size_t>(numNodeSlots), InvalidNode);
        properTail_.assign(static_cast<size_t>(numNodeSlots), InvalidNode);
        numProperPartsByNode_.assign(static_cast<size_t>(numNodeSlots), 0);
        initializeProperPartStorage(static_cast<size_t>(numProperParts));

        numNodes_ = static_cast<int>(numNodeSlots);
        rootNodeId_ = InvalidNode;

        for (NodeId properPartId = 0; properPartId < numProperParts; ++properPartId) {
            const NodeId ownerHigraNodeId = parent[static_cast<size_t>(properPartId)];
            if (ownerHigraNodeId < numProperParts || ownerHigraNodeId >= numHigraNodes) {
                throw std::invalid_argument("Each Higra leaf must point to an internal node.");
            }
            properPartOwner_[static_cast<size_t>(properPartId)] = ownerHigraNodeId - numProperParts;
        }

        for (NodeId slotId = 0; slotId < numNodeSlots; ++slotId) {
            const NodeId higraNodeId = numProperParts + slotId;
            const NodeId parentHigraNodeId = parent[static_cast<size_t>(higraNodeId)];
            if (parentHigraNodeId == higraNodeId) {
                if (rootNodeId_ != InvalidNode) {
                    throw std::invalid_argument("A Higra hierarchy must encode exactly one root.");
                }
                nodeParent_[static_cast<size_t>(slotId)] = slotId;
                rootNodeId_ = slotId;
                continue;
            }

            if (parentHigraNodeId < numProperParts || parentHigraNodeId >= numHigraNodes) {
                throw std::invalid_argument("Each Higra internal node must point to another internal node or to itself.");
            }
            const NodeId parentSlotId = parentHigraNodeId - numProperParts;
            linkChildSlot(parentSlotId, slotId);
        }

        if (rootNodeId_ == InvalidNode) {
            throw std::invalid_argument("A Higra hierarchy must encode exactly one root.");
        }
        rebuildProperPartLinksFromOwnership();
        invalidatePrePostOrderCache();
        invalidateAllIterators();
        preservesHigraNodeIdSpace_ = true;
    }

    /**
     * @brief Tag-protected import from native MAF topology buffers.
     *
     * This path is used by builders that already materialize internal-node
     * parent links and row-major proper-part owners, such as the SDRT builder.
     * The buffers are copied into canonical tree storage and validated as one
     * connected rooted hierarchy.
     */
    MorphologicalTree(detail::MorphologicalTreeConstructionTag, std::span<const NodeId> nodeParent, std::span<const NodeId> properPartOwner, NodeId root, int rows, int cols) {
        if (properPartOwner.size() != static_cast<size_t>(rows * cols)) {
            throw std::invalid_argument("Self-dual residual tree proper-part domain must match rows * cols.");
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

        treeType_ = MorphologicalTreeKind::SELF_DUAL_RESIDUAL_TREE;
        numRows_ = rows;
        numCols_ = cols;
        adj_ = std::nullopt;
        tosAdjacencyPolicy_ = std::nullopt;

        initializeEmptyStorage(static_cast<size_t>(numProperParts));
        nodeParent_.assign(nodeParent.begin(), nodeParent.end());
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
        properPartOwner_.assign(properPartOwner.begin(), properPartOwner.end());
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

        for (NodeId properPartId = 0; properPartId < numProperParts; ++properPartId) {
            const NodeId ownerId = properPartOwner_[static_cast<size_t>(properPartId)];
            if (ownerId < 0 || ownerId >= numNodeSlots) {
                throw std::invalid_argument("Native topology import found a proper-part owner outside the internal-node domain.");
            }
        }

        rebuildProperPartLinksFromOwnership();
        invalidatePrePostOrderCache();
        invalidateAllIterators();
        preservesHigraNodeIdSpace_ = false;
    }

    /**
     * @brief Creates an independent copy of the structural tree state.
     *
     * The public copy constructor stays deleted so ownership remains explicit at
     * API boundaries. This method is used by wrappers that need to preserve a
     * caller-owned topology while creating a new tree-backed object.
     */
    MorphologicalTree clone() const {
        MorphologicalTree cloned;
        cloned.rootNodeId_ = rootNodeId_;
        cloned.treeType_ = treeType_;
        cloned.numRows_ = numRows_;
        cloned.numCols_ = numCols_;
        cloned.adj_ = adj_;
        cloned.tosAdjacencyPolicy_ = tosAdjacencyPolicy_;
        cloned.numNodes_ = numNodes_;
        cloned.preservesHigraNodeIdSpace_ = preservesHigraNodeIdSpace_;
        cloned.editSessionOpen_ = false;
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
     */
    std::size_t getMutationVersion() const noexcept { return mutationVersion_;}

    /**
     * @brief Rejects stale read-only views that captured an older mutation version.
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
     */
    inline int getNumInternalNodeSlots() const { return static_cast<int>(nodeParent_.size()); }

    /**
     * @brief Returns the size of the proper-part domain.
     */
    inline int getNumTotalProperParts() const { return static_cast<int>(properPartOwner_.size()); }

    /**
     * @brief Returns the size of the preserved imported Higra node-id domain.
     *
     * @throws std::runtime_error if the tree was not imported from Higra or if
     * the preserved Higra node-id space was invalidated by an edit.
     */
    inline int getNumHigraNodes() const {
        if (!preservesHigraNodeIdSpace_) {
            throw std::runtime_error("This tree does not preserve an imported Higra node-id space.");
        }
        return getNumTotalProperParts() + getNumInternalNodeSlots();
    }

    /**
     * @brief Returns the size of the requested node-id domain.
     *
     * `NodeIdSpace::HIGRA` means the preserved imported Higra domain, not the
     * compact domain that would be generated by exporting the current tree.
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
     */
    inline NodeId getHigraNodeId(NodeId nodeId) const noexcept {
        if (!preservesHigraNodeIdSpace_ || !isNode(nodeId) || !isAlive(nodeId)) {
            return InvalidNode;
        }
        return getNumTotalProperParts() + nodeId;
    }

    /**
     * @brief Returns the current hierarchy root.
     */
    inline NodeId getRoot() const { return rootNodeId_; }

    /**
     * @brief Tests whether `nodeId` belongs to the internal-node id domain.
     */
    inline bool isNode(NodeId nodeId) const noexcept { return nodeId >= 0 && nodeId < static_cast<int>(nodeParent_.size()); }

    /**
     * @brief Tests whether `id` belongs to the proper-part domain.
     */
    inline bool isProperPart(NodeId id) const noexcept { return id >= 0 && id < static_cast<int>(properPartOwner_.size()); }

    /**
     * @brief Tests whether a node slot currently represents a live node.
     */
    inline bool isAlive(NodeId nodeId) const {
        if (!isNode(nodeId)) {
            return false;
        }
        const NodeId localId = nodeId;
        return localId >= 0
            && localId < static_cast<NodeId>(nodeParent_.size())
            && !isFreeSlot(localId)
            && (nodeParent_[localId] != InvalidNode || localId == rootNodeId_);
    }

    /**
     * @brief Tests whether `nodeId` is the current root.
     */
    inline bool isRoot(NodeId nodeId) const { return nodeId == getRoot(); }

    /**
     * @brief Returns the number of currently reusable node slots.
     */
    inline int getNumFreeNodeSlots() const { return static_cast<int>(freeNodeIds_.size()); }

    /**
     * @brief Counts the live nodes that currently have no children.
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
     */
    inline int getNumChildren(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getNumChildren");
        return numChildrenByNode_[nodeId];
    }

    /**
     * @brief Returns the number of internal descendants of `nodeId`.
     */
    inline int getNodeNumDescendants(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getNodeNumDescendants");
        ensurePrePostOrderCache();
        const NodeId localId = nodeId;
        return (prePostOrderCache_.timePostOrder[localId] - prePostOrderCache_.timePreOrder[localId] - 1) / 2;
    }

    /**
     * @brief Returns the number of siblings of `nodeId`.
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
     */
    inline int getNodeTimePreOrder(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getNodeTimePreOrder");
        ensurePrePostOrderCache();
        return prePostOrderCache_.timePreOrder[nodeId];
    }

    /**
     * @brief Returns the postorder time of `nodeId` in the cached DFS traversal.
     */
    inline int getNodeTimePostOrder(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getNodeTimePostOrder");
        ensurePrePostOrderCache();
        return prePostOrderCache_.timePostOrder[nodeId];
    }

    /**
     * @brief Returns the first direct child of `nodeId`, or `InvalidNode`.
     */
    inline NodeId getFirstChild(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getFirstChild");
        return firstChild_[nodeId];
    }

    /**
     * @brief Returns the next sibling of `nodeId`, or `InvalidNode`.
     */
    inline NodeId getNextSibling(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getNextSibling");
        return nextSibling_[nodeId];
    }

    /**
     * @brief Tests whether `nodeId` has no direct children.
     */
    inline bool isLeaf(NodeId nodeId) const { return getFirstChild(nodeId) == InvalidNode; }

    /**
     * @brief Returns the number of direct proper parts owned by `nodeId`.
     */
    inline int getNumProperParts(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getNumProperParts");
        return numProperPartsByNode_[nodeId];
    }

    /**
     * @brief Tests whether `childId` is a direct child of `parentNodeId`.
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
     */
    inline NodeId getProperPartOwner(NodeId properPartId) const {
        return isProperPart(properPartId)
            ? properPartOwner_[static_cast<size_t>(properPartId)]
            : InvalidNode;
    }

    /**
     * @brief Returns all live leaf nodes in the current hierarchy.
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
     * @brief Returns the current tree type.
     */
    inline MorphologicalTreeKind getTreeType() const noexcept{return treeType_;}

    /**
     * @brief Returns the number of currently live nodes.
     */
    inline int getNumNodes()const noexcept{ return numNodes_; }

    /**
     * @brief Tests whether a staged edit session is currently open.
     */
    inline bool isEditing() const noexcept { return editSessionOpen_; }

    /**
     * @brief Rejects operations that require a committed connected topology.
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
     * @brief Tests whether an adjacency relation is attached to the tree.
     */
    inline bool hasAdjacencyRelation() const noexcept { return static_cast<bool>(adj_); }

    /**
     * @brief Tests whether this image-built tree of shapes carries min/max auxiliary adjacency metadata.
     */
    inline bool hasTreeOfShapesAdjacencyPolicy() const noexcept { return static_cast<bool>(tosAdjacencyPolicy_); }

    /**
     * @brief Returns the auxiliary min-tree adjacency radius used by the ToS interpolation policy.
     */
    inline double getTreeOfShapesMinTreeAdjacencyRadius() const {
        if (!tosAdjacencyPolicy_) {
            throw std::runtime_error("Tree-of-shapes adjacency policy is not available.");
        }
        return tosAdjacencyPolicy_->first.getRadius();
    }

    /**
     * @brief Returns the auxiliary max-tree adjacency radius used by the ToS interpolation policy.
     */
    inline double getTreeOfShapesMaxTreeAdjacencyRadius() const {
        if (!tosAdjacencyPolicy_) {
            throw std::runtime_error("Tree-of-shapes adjacency policy is not available.");
        }
        return tosAdjacencyPolicy_->second.getRadius();
    }

    /**
     * @brief Returns the auxiliary min-tree adjacency relation used by the ToS interpolation policy.
     */
    inline const AdjacencyRelation* getTreeOfShapesMinTreeAdjacencyRelation() const noexcept {
        return tosAdjacencyPolicy_ ? &tosAdjacencyPolicy_->first : nullptr;
    }

    /**
     * @brief Returns the auxiliary max-tree adjacency relation used by the ToS interpolation policy.
     */
    inline const AdjacencyRelation* getTreeOfShapesMaxTreeAdjacencyRelation() const noexcept {
        return tosAdjacencyPolicy_ ? &tosAdjacencyPolicy_->second : nullptr;
    }

    /**
     * @brief Returns the number of image rows in the attached 2D domain.
     */
    inline int getNumRowsOfImage() const {
        return numRows_;
    }

    /**
     * @brief Returns the number of image columns in the attached 2D domain.
     */
    inline int getNumColsOfImage() const {
        return numCols_;
    }

    /**
     * @brief Returns the mutable adjacency relation attached to the tree, if any.
     */
    inline AdjacencyRelation* getAdjacencyRelation() noexcept {return adj_ ? &*adj_ : nullptr;}

    /**
     * @brief Returns the read-only adjacency relation attached to the tree, if any.
     */
    inline const AdjacencyRelation* getAdjacencyRelation() const noexcept {return adj_ ? &*adj_ : nullptr;}

    /**
     * @brief Validates that the current structure is one connected rooted tree.
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
        std::vector<int> parentChainStamp(static_cast<size_t>(numSlots), 0);
        int aliveCount = 0;
        int stamp = 0;

        for (NodeId nodeId = 0; nodeId < numSlots; ++nodeId) {
            if (!isAlive(nodeId)) {
                continue;
            }
            ++aliveCount;
            ++stamp;

            NodeId currentNodeId = nodeId;
            while (true) {
                if (!isNode(currentNodeId) || !isAlive(currentNodeId)) {
                    throw std::runtime_error("Connected-tree validation found a parent chain that leaves the alive node domain.");
                }
                if (parentChainStamp[static_cast<size_t>(currentNodeId)] == stamp) {
                    throw std::runtime_error("Connected-tree validation found a cycle in the parent chain.");
                }
                parentChainStamp[static_cast<size_t>(currentNodeId)] = stamp;

                const NodeId parentNodeId = nodeParent_[static_cast<size_t>(currentNodeId)];
                if (parentNodeId == InvalidNode) {
                    throw std::runtime_error("Connected-tree validation found an alive node with no parent.");
                }
                if (parentNodeId == currentNodeId) {
                    if (currentNodeId != rootNodeId_) {
                        throw std::runtime_error("Connected-tree validation found a detached alive node.");
                    }
                    break;
                }
                if (!isNode(parentNodeId) || !isAlive(parentNodeId)) {
                    throw std::runtime_error("Connected-tree validation found an alive node whose parent is outside the alive node domain.");
                }
                currentNodeId = parentNodeId;
            }

            if (nodeId != rootNodeId_) {
                expectedChildrenByNode[static_cast<size_t>(nodeParent_[static_cast<size_t>(nodeId)])] += 1;
            }
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
            for (NodeId childNodeId = firstChild_[static_cast<size_t>(nodeId)];
                 childNodeId != InvalidNode;
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
            for (NodeId properPartId = properHead_[static_cast<size_t>(nodeId)];
                 properPartId != InvalidNode;
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
    }

    /**
     * @brief Runs strong validation and returns the result instead of throwing.
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
     */
    inline bool isAncestor(NodeId u, NodeId v) const {
        requireAliveNode(u, "MorphologicalTree::isAncestor");
        requireAliveNode(v, "MorphologicalTree::isAncestor");
        ensurePrePostOrderCache();
        const NodeId slotU = u;
        const NodeId slotV = v;
        return prePostOrderCache_.timePreOrder[slotU] <= prePostOrderCache_.timePreOrder[slotV]
            && prePostOrderCache_.timePostOrder[slotU] >= prePostOrderCache_.timePostOrder[slotV];
    }

    /**
     * @brief Tests whether `u` is a descendant of `v`.
     */
    inline bool isDescendant(NodeId u, NodeId v) const {
        requireAliveNode(u, "MorphologicalTree::isDescendant");
        requireAliveNode(v, "MorphologicalTree::isDescendant");
        ensurePrePostOrderCache();
        const NodeId slotU = u;
        const NodeId slotV = v;
        return prePostOrderCache_.timePreOrder[slotV] <= prePostOrderCache_.timePreOrder[slotU]
            && prePostOrderCache_.timePostOrder[slotV] >= prePostOrderCache_.timePostOrder[slotU];
    }
    /**
     * @brief Tests whether `u` and `v` are comparable in the ancestry order.
     */
    inline bool isComparable(NodeId u, NodeId v) const {return isAncestor(u, v) || isAncestor(v, u);}

    /**
     * @brief Tests whether `u` is a strict ancestor of `v`.
     */
    inline bool isStrictAncestor(NodeId u, NodeId v) const {return u != v && isAncestor(u, v);}

    /**
     * @brief Tests whether `u` is a strict descendant of `v`.
     */
    inline bool isStrictDescendant(NodeId u, NodeId v) const { return u != v && isDescendant(u, v);}

    /**
     * @brief Tests whether `u` and `v` are strictly comparable in the ancestry order.
     */
    inline bool isStrictComparable(NodeId u, NodeId v) const { return isStrictAncestor(u, v) || isStrictAncestor(v, u);}

    /**
     * @brief Returns the lowest common ancestor of `u` and `v`.
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
     */
    inline AliveNodeRange getAliveNodeIds() const {
        return AliveNodeRange(this, 0, getNumInternalNodeSlots(), nodeStructureVersion_);
    }

    /**
     * @brief Returns a fail-fast range over the direct children of `nodeId`.
     */
    inline ChildrenRange getChildren(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getChildren");
        return ChildrenRange(this, firstChild_[nodeId], topologyVersion_);
    }

    /**
     * @brief Returns a fail-fast range over the direct proper parts of `nodeId`.
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
     */
    inline ConnectedComponentRange getConnectedComponent(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getConnectedComponent");
        return ConnectedComponentRange(this, nodeId, topologyVersion_, properPartVersion_);
    }

    /**
     * @brief Returns a post-order traversal range rooted at the connected root.
     */
    inline PostOrderNodeRange getPostOrderNodes() const {
        return PostOrderNodeRange(this, getRoot(), topologyVersion_);
    }

    /**
     * @brief Returns a post-order traversal range rooted at `rootNodeId`.
     */
    inline PostOrderNodeRange getPostOrderNodes(NodeId rootNodeId) const {
        requireAliveNode(rootNodeId, "MorphologicalTree::getPostOrderNodes");
        return PostOrderNodeRange(this, rootNodeId, topologyVersion_);
    }

    /**
     * @brief Returns a breadth-first traversal range rooted at the connected root.
     */
    inline BreadthFirstNodeRange getIteratorBreadthFirstTraversal() const {
        return BreadthFirstNodeRange(this, getRoot(), topologyVersion_);
    }

    /**
     * @brief Returns a breadth-first traversal range rooted at `rootNodeId`.
     */
    inline BreadthFirstNodeRange getIteratorBreadthFirstTraversal(NodeId rootNodeId) const {
        requireAliveNode(rootNodeId, "MorphologicalTree::getIteratorBreadthFirstTraversal");
        return BreadthFirstNodeRange(this, rootNodeId, topologyVersion_);
    }

    /**
     * @brief Returns the path from `nodeId` to the connected root.
     */
    inline PathToRootRange getPathToRootNodes(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getPathToRootNodes");
        return PathToRootRange(this, nodeId, topologyVersion_);
    }

    /**
     * @brief Returns the path that connects `sourceNodeId` and `targetNodeId`.
     */
    inline PathBetweenNodesRange getPathBetweenNodes(NodeId sourceNodeId, NodeId targetNodeId) const {
        requireAliveNode(sourceNodeId, "MorphologicalTree::getPathBetweenNodes");
        requireAliveNode(targetNodeId, "MorphologicalTree::getPathBetweenNodes");
        return PathBetweenNodesRange(this, sourceNodeId, targetNodeId, topologyVersion_);
    }

    /**
     * @brief Returns a pre-order traversal range over the subtree of `nodeId`.
     */
    inline SubtreeNodeRange getNodeSubtree(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getNodeSubtree");
        return SubtreeNodeRange(this, nodeId, topologyVersion_);
    }

    /**
     * @brief Returns a range over all proper descendants of `nodeId`.
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
     */
    [[nodiscard("Discarding the editor leaves the tree edit session open")]] TreeEditor edit();

    /**
     * @brief Prunes the subtree of `nodeId`, moving all its support to the parent.
     *
     * This is a safe committed edit: it is a complete local operation and does
     * not intentionally leave the tree in a staged disconnected state, but it
     * still advances the tree mutation version and invalidates derived state.
     */
    inline void pruneNode(NodeId nodeId) {
        requireNotEditing("MorphologicalTree::pruneNode");
        requireAliveNonRootNode(nodeId, "MorphologicalTree::pruneNode");
        const NodeId parentNodeId = getNodeParent(nodeId);
        if (parentNodeId == InvalidNode || parentNodeId == nodeId) {
            throw std::invalid_argument("MorphologicalTree::pruneNode requires an attached non-root node.");
        }
        const NodeId parentSlotId = parentNodeId;

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
                unlinkChildSlot(currentParentSlotId, currentId, false);
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
        std::vector<NodeId> euler_;
        std::vector<int> depth_;
        std::vector<int> firstOccurrence_;
        std::vector<int> log2_;
        std::vector<int> sparseTable_;
        int sparseTableStride_ = 0;
        const MorphologicalTree* tree_ = nullptr;

        /**
         * @brief Appends one DFS step to the Euler tour and recurses into children.
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
         */
        std::size_t sparseTableIndex(int row, int col) const {
            return static_cast<std::size_t>(row) * static_cast<std::size_t>(sparseTableStride_) +
                   static_cast<std::size_t>(col);
        }

        /**
         * @brief Returns the Euler-tour position of minimum depth on `[left, right]`.
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
        const MorphologicalTree* T_ = nullptr;
        NodeId current_ = InvalidNode;
        NodeId end_ = InvalidNode;
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
         */
        AliveNodeIterator(const MorphologicalTree* tree, NodeId current, NodeId end, std::size_t expectedVersion)
            : T_(tree), current_(current), end_(end), expectedVersion_(expectedVersion) {
            settle_();
        }

        /**
         * @brief Advances to the next live node slot.
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
         */
        NodeId operator*() const {
            T_->checkNodeIteratorVersion(expectedVersion_);
            return current_;
        }

        /**
         * @brief Compares iterator positions.
         */
        bool operator==(const AliveNodeIterator& other) const { return current_ == other.current_; }

        /**
         * @brief Compares iterator positions for inequality.
         */
        bool operator!=(const AliveNodeIterator& other) const { return !(*this == other); }
    };

    /**
     * @brief Range wrapper for iterating over live node ids.
     */
    class AliveNodeRange {
    private:
        const MorphologicalTree* T_ = nullptr;
        NodeId begin_ = InvalidNode;
        NodeId end_ = InvalidNode;
        std::size_t expectedVersion_ = 0;

    public:
        /**
         * @brief Creates an empty live-node range.
         */
        AliveNodeRange() = default;

        /**
         * @brief Creates a fail-fast live-node range over dense slots.
         */
        AliveNodeRange(const MorphologicalTree* tree, NodeId begin, NodeId end, std::size_t expectedVersion)
            : T_(tree), begin_(begin), end_(end), expectedVersion_(expectedVersion) {}

        /**
         * @brief Returns an iterator positioned at the first live slot.
         */
        AliveNodeIterator begin() const { return AliveNodeIterator(T_, begin_, end_, expectedVersion_); }

        /**
         * @brief Returns the sentinel iterator for the live-node range.
         */
        AliveNodeIterator end() const { return AliveNodeIterator(T_, InvalidNode, end_, expectedVersion_); }
    };

    /**
     * @brief Iterator over the direct children of one node.
     */
    class ChildrenIterator {
    private:
        const MorphologicalTree* T_ = nullptr;
        NodeId currentLocal_ = InvalidNode;
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
         */
        ChildrenIterator(const MorphologicalTree* tree, NodeId currentLocal, std::size_t expectedVersion)
            : T_(tree), currentLocal_(currentLocal), expectedVersion_(expectedVersion) {}

        /**
         * @brief Advances to the next sibling in the child list.
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
         */
        NodeId operator*() const {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            return currentLocal_;
        }

        /**
         * @brief Compares child iterator positions.
         */
        bool operator==(const ChildrenIterator& other) const { return currentLocal_ == other.currentLocal_; }

        /**
         * @brief Compares child iterator positions for inequality.
         */
        bool operator!=(const ChildrenIterator& other) const { return !(*this == other); }
    };

    /**
     * @brief Range wrapper for direct-child iteration.
     */
    class ChildrenRange {
    private:
        const MorphologicalTree* T_ = nullptr;
        NodeId firstLocal_ = InvalidNode;
        std::size_t expectedVersion_ = 0;

    public:
        /**
         * @brief Creates an empty child range.
         */
        ChildrenRange() = default;

        /**
         * @brief Creates a range over a linked list of direct children.
         */
        ChildrenRange(const MorphologicalTree* tree, NodeId firstLocal, std::size_t expectedVersion)
            : T_(tree), firstLocal_(firstLocal), expectedVersion_(expectedVersion) {}

        /**
         * @brief Returns an iterator at the first child.
         */
        ChildrenIterator begin() const { return ChildrenIterator(T_, firstLocal_, expectedVersion_); }

        /**
         * @brief Returns the child-range sentinel iterator.
         */
        ChildrenIterator end() const { return ChildrenIterator(T_, InvalidNode, expectedVersion_); }
    };

    /**
     * @brief Iterator over the direct proper parts owned by one node.
     */
    class ProperPartsIterator {
    private:
        const MorphologicalTree* T_ = nullptr;
        NodeId currentPixel_ = InvalidNode;
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
         */
        ProperPartsIterator(const MorphologicalTree* tree, NodeId currentPixel, std::size_t expectedVersion)
            : T_(tree), currentPixel_(currentPixel), expectedVersion_(expectedVersion) {}

        /**
         * @brief Advances to the next proper part in the owner list.
         */
        ProperPartsIterator& operator++() {
            T_->checkProperPartIteratorVersion(expectedVersion_);
            if (T_ && currentPixel_ != InvalidNode) {
                currentPixel_ = T_->nextProperPart_[currentPixel_];
            }
            return *this;
        }

        /**
         * @brief Returns the current proper-part id.
         */
        NodeId operator*() const {
            T_->checkProperPartIteratorVersion(expectedVersion_);
            return currentPixel_;
        }

        /**
         * @brief Compares proper-part iterator positions.
         */
        bool operator==(const ProperPartsIterator& other) const { return currentPixel_ == other.currentPixel_; }

        /**
         * @brief Compares proper-part iterator positions for inequality.
         */
        bool operator!=(const ProperPartsIterator& other) const { return !(*this == other); }
    };

    /**
     * @brief Range wrapper for direct proper-part iteration.
     */
    class ProperPartsRange {
    private:
        const MorphologicalTree* T_ = nullptr;
        NodeId firstPixel_ = InvalidNode;
        std::size_t expectedVersion_ = 0;

    public:
        /**
         * @brief Creates an empty direct proper-part range.
         */
        ProperPartsRange() = default;

        /**
         * @brief Creates a range over one node's direct proper-part list.
         */
        ProperPartsRange(const MorphologicalTree* tree, NodeId firstPixel, std::size_t expectedVersion)
            : T_(tree), firstPixel_(firstPixel), expectedVersion_(expectedVersion) {}

        /**
         * @brief Returns an iterator at the first direct proper part.
         */
        ProperPartsIterator begin() const { return ProperPartsIterator(T_, firstPixel_, expectedVersion_); }

        /**
         * @brief Returns the direct proper-part range sentinel.
         */
        ProperPartsIterator end() const { return ProperPartsIterator(T_, InvalidNode, expectedVersion_); }
    };

    /**
     * @brief Iterator over all proper parts in one connected component.
     */
    class ConnectedComponentIterator {
    private:
        const MorphologicalTree* T_ = nullptr;
        std::vector<NodeId> nodeStack_;
        NodeId currentPixel_ = InvalidNode;
        std::size_t expectedTopologyVersion_ = 0;
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
            while (currentPixel_ == InvalidNode && !nodeStack_.empty()) {
                const NodeId nodeId = nodeStack_.back();
                nodeStack_.pop_back();
                pushChildren(nodeId);
                currentPixel_ = T_->properHead_[static_cast<size_t>(nodeId)];
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
         */
        ConnectedComponentIterator(
            const MorphologicalTree* tree,
            NodeId rootNodeId,
            std::size_t expectedTopologyVersion,
            std::size_t expectedProperPartVersion)
            : T_(tree),
              expectedTopologyVersion_(expectedTopologyVersion),
              expectedProperPartVersion_(expectedProperPartVersion)
        {
            if (T_ && rootNodeId != InvalidNode) {
                nodeStack_.push_back(rootNodeId);
                settle();
            }
        }

        /**
         * @brief Advances to the next proper part in the connected component.
         */
        ConnectedComponentIterator& operator++() {
            checkVersions();
            if (currentPixel_ != InvalidNode) {
                currentPixel_ = T_->nextProperPart_[static_cast<size_t>(currentPixel_)];
            }
            settle();
            return *this;
        }

        /**
         * @brief Returns the current proper-part id.
         */
        NodeId operator*() const {
            checkVersions();
            return currentPixel_;
        }

        /**
         * @brief Compares connected-component iterator positions.
         */
        bool operator==(const ConnectedComponentIterator& other) const { return currentPixel_ == other.currentPixel_; }

        /**
         * @brief Compares connected-component iterator positions for inequality.
         */
        bool operator!=(const ConnectedComponentIterator& other) const { return !(*this == other); }
    };

    /**
     * @brief Range wrapper for connected-component proper-part iteration.
     */
    class ConnectedComponentRange {
    private:
        const MorphologicalTree* T_ = nullptr;
        NodeId rootNodeId_ = InvalidNode;
        std::size_t expectedTopologyVersion_ = 0;
        std::size_t expectedProperPartVersion_ = 0;

    public:
        /**
         * @brief Creates an empty connected-component range.
         */
        ConnectedComponentRange() = default;

        /**
         * @brief Creates a range over all proper parts in the subtree of `rootNodeId`.
         */
        ConnectedComponentRange(
            const MorphologicalTree* tree,
            NodeId rootNodeId,
            std::size_t expectedTopologyVersion,
            std::size_t expectedProperPartVersion)
            : T_(tree),
              rootNodeId_(rootNodeId),
              expectedTopologyVersion_(expectedTopologyVersion),
              expectedProperPartVersion_(expectedProperPartVersion) {}

        /**
         * @brief Returns an iterator at the first proper part in the component.
         */
        ConnectedComponentIterator begin() const {
            return ConnectedComponentIterator(T_, rootNodeId_, expectedTopologyVersion_, expectedProperPartVersion_);
        }

        /**
         * @brief Returns the connected-component range sentinel.
         */
        ConnectedComponentIterator end() const { return ConnectedComponentIterator(); }
    };

    /**
     * @brief Post-order iterator over one subtree.
     */
    class PostOrderNodeIterator {
    private:
        struct Item { NodeId id; bool expanded; };
        const MorphologicalTree* T_ = nullptr;
        std::vector<Item> stack_;
        NodeId current_ = InvalidNode;
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
         */
        PostOrderNodeIterator(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedVersion) : T_(tree), expectedVersion_(expectedVersion) {
            if (T_ && rootNodeId != InvalidNode) {
                stack_.push_back({rootNodeId, false});
                settle_();
            }
        }

        /**
         * @brief Advances to the next node in post-order.
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
         */
        NodeId operator*() const {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            return current_;
        }

        /**
         * @brief Compares post-order iterator positions.
         */
        bool operator==(const PostOrderNodeIterator& other) const { return current_ == other.current_; }

        /**
         * @brief Compares post-order iterator positions for inequality.
         */
        bool operator!=(const PostOrderNodeIterator& other) const { return !(*this == other); }
    };

    /**
     * @brief Range wrapper for post-order subtree traversal.
     */
    class PostOrderNodeRange {
    private:
        const MorphologicalTree* T_ = nullptr;
        NodeId rootNodeId_ = InvalidNode;
        std::size_t expectedVersion_ = 0;

    public:
        /**
         * @brief Creates an empty post-order range.
         */
        PostOrderNodeRange() = default;

        /**
         * @brief Creates a post-order range rooted at `rootNodeId`.
         */
        PostOrderNodeRange(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedVersion)
            : T_(tree), rootNodeId_(rootNodeId), expectedVersion_(expectedVersion) {}

        /**
         * @brief Returns an iterator at the first post-order node.
         */
        PostOrderNodeIterator begin() const { return PostOrderNodeIterator(T_, rootNodeId_, expectedVersion_); }

        /**
         * @brief Returns the post-order range sentinel.
         */
        PostOrderNodeIterator end() const { return PostOrderNodeIterator(); }
    };

    /**
     * @brief Breadth-first iterator over one subtree.
     */
    class BreadthFirstNodeIterator {
    private:
        const MorphologicalTree* T_ = nullptr;
        FastQueue<NodeId> queue_;
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
         */
        BreadthFirstNodeIterator(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedVersion) : T_(tree), expectedVersion_(expectedVersion) {
            if (T_ && rootNodeId != InvalidNode) {
                queue_.push(rootNodeId);
            }
        }

        /**
         * @brief Advances to the next node in breadth-first order.
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
         */
        NodeId operator*() const {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            return queue_.front();
        }

        /**
         * @brief Compares breadth-first iterator exhaustion state.
         */
        bool operator==(const BreadthFirstNodeIterator& other) const { return queue_.empty() == other.queue_.empty(); }

        /**
         * @brief Compares breadth-first iterator exhaustion state for inequality.
         */
        bool operator!=(const BreadthFirstNodeIterator& other) const { return !(*this == other); }
    };

    /**
     * @brief Range wrapper for breadth-first subtree traversal.
     */
    class BreadthFirstNodeRange {
    private:
        const MorphologicalTree* T_ = nullptr;
        NodeId rootNodeId_ = InvalidNode;
        std::size_t expectedVersion_ = 0;

    public:
        /**
         * @brief Creates an empty breadth-first range.
         */
        BreadthFirstNodeRange() = default;

        /**
         * @brief Creates a breadth-first range rooted at `rootNodeId`.
         */
        BreadthFirstNodeRange(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedVersion)
            : T_(tree), rootNodeId_(rootNodeId), expectedVersion_(expectedVersion) {}

        /**
         * @brief Returns an iterator at the first breadth-first node.
         */
        BreadthFirstNodeIterator begin() const { return BreadthFirstNodeIterator(T_, rootNodeId_, expectedVersion_); }

        /**
         * @brief Returns the breadth-first range sentinel.
         */
        BreadthFirstNodeIterator end() const { return BreadthFirstNodeIterator(); }
    };

    /**
     * @brief Iterator that walks from a node towards the root.
     */
    class PathToRootIterator {
    private:
        const MorphologicalTree* T_ = nullptr;
        NodeId current_ = InvalidNode;
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
         */
        PathToRootIterator(const MorphologicalTree* tree, NodeId current, std::size_t expectedVersion)
            : T_(tree), current_(current), expectedVersion_(expectedVersion) {}

        /**
         * @brief Advances one step toward the root.
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
         */
        NodeId operator*() const {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            return current_;
        }

        /**
         * @brief Compares path-to-root iterator positions.
         */
        bool operator==(const PathToRootIterator& other) const { return current_ == other.current_; }

        /**
         * @brief Compares path-to-root iterator positions for inequality.
         */
        bool operator!=(const PathToRootIterator& other) const { return !(*this == other); }
    };

    /**
     * @brief Range wrapper for rootward path traversal.
     */
    class PathToRootRange {
    private:
        const MorphologicalTree* T_ = nullptr;
        NodeId start_ = InvalidNode;
        std::size_t expectedVersion_ = 0;

    public:
        /**
         * @brief Creates an empty rootward-path range.
         */
        PathToRootRange() = default;

        /**
         * @brief Creates a range from `start` to the connected root.
         */
        PathToRootRange(const MorphologicalTree* tree, NodeId start, std::size_t expectedVersion)
            : T_(tree), start_(start), expectedVersion_(expectedVersion) {}

        /**
         * @brief Returns an iterator at the path start node.
         */
        PathToRootIterator begin() const { return PathToRootIterator(T_, start_, expectedVersion_); }

        /**
         * @brief Returns the rootward-path range sentinel.
         */
        PathToRootIterator end() const { return PathToRootIterator(); }
    };

    /**
     * @brief Iterator over a materialised path between two nodes.
     */
    class PathBetweenNodesIterator {
    private:
        const MorphologicalTree* T_ = nullptr;
        const std::vector<NodeId>* path_ = nullptr;
        std::size_t index_ = 0;
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
         */
        PathBetweenNodesIterator(
            const MorphologicalTree* tree,
            const std::vector<NodeId>* path,
            std::size_t index,
            std::size_t expectedVersion)
            : T_(tree), path_(path), index_(index), expectedVersion_(expectedVersion) {}

        /**
         * @brief Advances to the next node in the materialised path.
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
         */
        NodeId operator*() const {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            return (*path_)[index_];
        }

        /**
         * @brief Compares path iterator positions.
         */
        bool operator==(const PathBetweenNodesIterator& other) const {
            return path_ == other.path_ && index_ == other.index_;
        }

        /**
         * @brief Compares path iterator positions for inequality.
         */
        bool operator!=(const PathBetweenNodesIterator& other) const { return !(*this == other); }
    };

    /**
     * @brief Range wrapper for the path connecting two nodes in the same component.
     */
    class PathBetweenNodesRange {
    private:
        const MorphologicalTree* T_ = nullptr;
        std::vector<NodeId> path_;
        std::size_t expectedVersion_ = 0;

        /**
         * @brief Materialises the path joining two nodes in the same connected component.
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
         */
        PathBetweenNodesRange(
            const MorphologicalTree* tree,
            NodeId sourceNodeId,
            NodeId targetNodeId,
            std::size_t expectedVersion)
            : T_(tree),
              path_(buildPath(tree, sourceNodeId, targetNodeId)),
              expectedVersion_(expectedVersion) {}

        /**
         * @brief Returns an iterator at the first node in the materialised path.
         */
        PathBetweenNodesIterator begin() const {
            return PathBetweenNodesIterator(T_, &path_, 0, expectedVersion_);
        }

        /**
         * @brief Returns the materialised path range sentinel.
         */
        PathBetweenNodesIterator end() const {
            return PathBetweenNodesIterator(T_, &path_, path_.size(), expectedVersion_);
        }
    };

    /**
     * @brief Depth-first iterator over a subtree in pre-order.
     */
    class SubtreeNodeIterator {
    private:
        const MorphologicalTree* T_ = nullptr;
        std::vector<NodeId> stack_;
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
         */
        SubtreeNodeIterator(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedVersion) : T_(tree), expectedVersion_(expectedVersion) {
            if (T_ && rootNodeId != InvalidNode) {
                stack_.push_back(rootNodeId);
            }
        }

        /**
         * @brief Advances to the next node in pre-order subtree traversal.
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
         */
        NodeId operator*() const {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            return stack_.back();
        }

        /**
         * @brief Compares subtree iterator exhaustion state.
         */
        bool operator==(const SubtreeNodeIterator& other) const { return stack_.empty() == other.stack_.empty(); }

        /**
         * @brief Compares subtree iterator exhaustion state for inequality.
         */
        bool operator!=(const SubtreeNodeIterator& other) const { return !(*this == other); }
    };

    /**
     * @brief Range wrapper for pre-order subtree traversal.
     */
    class SubtreeNodeRange {
    private:
        const MorphologicalTree* T_ = nullptr;
        NodeId rootNodeId_ = InvalidNode;
        std::size_t expectedVersion_ = 0;

    public:
        /**
         * @brief Creates an empty subtree range.
         */
        SubtreeNodeRange() = default;

        /**
         * @brief Creates a pre-order range over the subtree rooted at `rootNodeId`.
         */
        SubtreeNodeRange(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedVersion)
            : T_(tree), rootNodeId_(rootNodeId), expectedVersion_(expectedVersion) {}

        /**
         * @brief Returns an iterator at the subtree root.
         */
        SubtreeNodeIterator begin() const { return SubtreeNodeIterator(T_, rootNodeId_, expectedVersion_); }

        /**
         * @brief Returns the subtree range sentinel.
         */
        SubtreeNodeIterator end() const { return SubtreeNodeIterator(); }
    };

    /**
     * @brief Range wrapper over the proper descendants of one node.
     */
    class DescendantNodeRange {
    private:
        const MorphologicalTree* T_ = nullptr;
        NodeId rootNodeId_ = InvalidNode;
        std::size_t expectedVersion_ = 0;

    public:
        /**
         * @brief Creates an empty descendant range.
         */
        DescendantNodeRange() = default;

        /**
         * @brief Creates a range over proper descendants of `rootNodeId`.
         */
        DescendantNodeRange(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedVersion)
            : T_(tree), rootNodeId_(rootNodeId), expectedVersion_(expectedVersion) {}

        /**
         * @brief Returns an iterator at the first proper descendant.
         */
        SubtreeNodeIterator begin() const {
            auto it = SubtreeNodeIterator(T_, rootNodeId_, expectedVersion_);
            ++it;
            return it;
        }

        /**
         * @brief Returns the descendant range sentinel.
         */
        SubtreeNodeIterator end() const { return SubtreeNodeIterator(); }
    };



};

} // namespace mmcfilters
