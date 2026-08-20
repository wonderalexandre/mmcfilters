#pragma once

#include "../utils/RegularGridAdjacency2D.hpp"
#include "../utils/Assert.hpp"
#include "../utils/Altitude.hpp"
#include "../utils/Common.hpp"
#include "../utils/Contract.hpp"
#include "../dataStructure/FastQueue.hpp"
#include "MorphologicalTreeSemantics.hpp"
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

namespace detail {
class CommittedTreeAccess;
}

/**
 * @brief Selects the node-id domain used by attribute buffers exposed to callers.
 */
enum class NodeIdSpace { MorphologicalTree, Higra };

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
 * @brief Mutable connected-subset tree on a finite pixel domain.
 *
 * `MorphologicalTree` is the central mutable hierarchy of this project. It can
 * represent any rooted connected-subset tree model supported by the project.
 * A morphological tree of partial partitions is the stricter case in which
 * every node has a non-empty proper part. The class exposes a dense `NodeId`
 * domain together with the smallest-node map, explicit parent/child links,
 * linked pixel lists for the nodes' proper parts, optional regular 2D pixel
 * metadata, and a small set of
 * structural caches used by the
 * public traversal and ancestry queries.
 *
 * Data model:
 *
 * - pixels are indexed in the range `[0, numPixels())`;
 * - internal nodes are indexed by `NodeId` in the range `[0, numInternalNodeSlots())`;
 * - `properPart(node)` contains the pixels whose smallest node is `node` and may be empty;
 * - `nodeSupport(node)` is the union of the proper parts in the node's subtree;
 * - every node in a checked committed hierarchy has non-empty node support;
 * - the root and detached nodes point to themselves as parent.
 *
 * Main responsibilities:
 *
 * - expose mutable topology operations used by filters and adjusters;
 * - provide structural traversals over live nodes, children, proper parts, node supports, and subtrees;
 * - serve as the topology/smallest-node core used by valued wrappers such as
 *   `ValuedMorphologicalTree<std::uint8_t>`.
 *
 * The class intentionally keeps only the canonical structural state and a small
 * number of derived caches. Higher-level attribute computation is delegated to
 * the incremental attribute computers in `mmcfilters/attributes`.
 */
class MorphologicalTree {
  private:
    friend class TreeEditor;
    friend class MorphologicalTreeFactory;
    friend class detail::CommittedTreeAccess;

    class LCAEulerRMQ; // Forward declaration for the LCA cache implementation.

    // ========================= Private attributes ========================= //
    /** @brief Dense node identifier of the root node identifier. */
    NodeId rootNodeId_ = InvalidNode;
    /** @brief Semantics. */
    MorphologicalTreeSemantics semantics_;
    /** @brief Grid domain2 d. */
    std::optional<GridDomain2D> gridDomain2D_;
    /** @brief Number of nodes. */
    int numNodes_ = 0;
    /** @brief Dense node identifier of the preserved external node identifier offset. */
    std::optional<NodeId> preservedExternalNodeIdOffset_;
    /** @brief Indicates whether a guarded edit session is open. */
    bool editSessionOpen_ = false;
    /** @brief Edit validation statistics. */
    TreeEditValidationStatistics editValidationStatistics_;

    // Smallest-node map, indexed by pixel id [0, numPixels()).
    /** @brief Dense node identifier of the smallest node containing each pixel. */
    std::vector<NodeId> smallestNodeMap_;

    // Parent links, indexed by local node-slot id [0, numInternalNodeSlots()).
    /** @brief Dense node identifier of the node parent. */
    std::vector<NodeId> nodeParent_;

    // Internal hierarchy linked structure, indexed by local node-slot id [0, numInternalNodeSlots()).
    /** @brief Dense node identifier of the first child. */
    std::vector<NodeId> firstChild_;
    /** @brief Dense node identifier of the next sibling. */
    std::vector<NodeId> nextSibling_;
    /** @brief Dense node identifier of the previous sibling. */
    std::vector<NodeId> prevSibling_;
    /** @brief Dense node identifier of the last child. */
    std::vector<NodeId> lastChild_;
    /** @brief Number of children by node. */
    std::vector<int> numChildrenByNode_;

    // Free slot management and alive-node iteration, indexed by local node-slot id [0, numInternalNodeSlots()).
    /** @brief Alive buffer. */
    std::vector<uint8_t> alive_;
    /** @brief Dense node identifier of the free node identifiers. */
    std::vector<NodeId> freeNodeIds_;

    // Per-node proper-part pixel lists, indexed by local node-slot id.
    /** @brief Pixel identifier of the proper head. */
    std::vector<PixelId> properHead_;
    /** @brief Pixel identifier of the proper tail. */
    std::vector<PixelId> properTail_;
    /** @brief Number of proper parts by node. */
    std::vector<int> properPartCardinalityByNode_;
    /** @brief Pixel identifier of the next proper part. */
    std::vector<PixelId> nextProperPart_;
    /** @brief Pixel identifier of the previous proper part. */
    std::vector<PixelId> prevProperPart_;

    // Structural caches for traversal-based queries, indexed by local node-slot id [0, numInternalNodeSlots()).
    /** @brief Caches interleaved DFS entry/exit event indices for the current topology version. */
    struct DfsIntervalCache {
        /** @brief Stores each node's zero-based DFS entry-event index. */
        std::vector<int> entryIndex;
        /** @brief Stores each node's zero-based DFS exit-event index. */
        std::vector<int> exitIndex;
        /** @brief Indicates whether the cached DFS intervals are valid. */
        bool valid = false;
        /** @brief Marks the cached DFS intervals as stale. */
        void invalidate() noexcept { valid = false; }
    };
    /** @brief DFS-interval cache. */
    mutable DfsIntervalCache dfsIntervalCache_;
    /** @brief Lowest-common-ancestor cache. */
    mutable std::unique_ptr<LCAEulerRMQ> lcaCache_;

    /** @brief Caches support cardinalities and row-major spatial keys for the current committed tree. */
    struct NodeSupportMetadataCache {
        /** @brief Number of pixels in every live node support. */
        std::vector<std::int32_t> cardinalityByNode;
        /** @brief Smallest row-major pixel in every live node support. */
        std::vector<PixelId> smallestPixelByNode;
        /** @brief Node-structure version represented by this cache. */
        std::size_t nodeStructureVersion = std::numeric_limits<std::size_t>::max();
        /** @brief Topology version represented by this cache. */
        std::size_t topologyVersion = std::numeric_limits<std::size_t>::max();
        /** @brief Proper-part version represented by this cache. */
        std::size_t properPartVersion = std::numeric_limits<std::size_t>::max();

        /** @brief Marks the cache as stale while retaining reusable storage. */
        void invalidate() noexcept {
            nodeStructureVersion = std::numeric_limits<std::size_t>::max();
            topologyVersion = std::numeric_limits<std::size_t>::max();
            properPartVersion = std::numeric_limits<std::size_t>::max();
        }
    };
    /** @brief Stores reusable node-support comparison metadata. */
    mutable NodeSupportMetadataCache nodeSupportMetadataCache_;

    // Version counters for iterator invalidation.
    /** @brief Node structure version used to detect stale derived state. */
    std::size_t nodeStructureVersion_ = 0;
    /** @brief Topology version used to detect stale derived state. */
    std::size_t topologyVersion_ = 0;
    /** @brief Proper part version used to detect stale derived state. */
    std::size_t properPartVersion_ = 0;
    /** @brief Mutation version used to detect stale derived state. */
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
            properHead_[slotId] = InvalidPixel;
            properTail_[slotId] = InvalidPixel;
            properPartCardinalityByNode_[slotId] = 0;
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
        properHead_.push_back(InvalidPixel);
        properTail_.push_back(InvalidPixel);
        properPartCardinalityByNode_.push_back(0);
        return slotId;
    }

    /**
     * @brief Resets the proper-part linked-list storage to the requested size.
     *
     * @param numPixels Number of pixels in the finite domain.
     */
    inline void initializeProperPartStorage(size_t numPixels) {
        nextProperPart_.assign(numPixels, InvalidPixel);
        prevProperPart_.assign(numPixels, InvalidPixel);
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
     * @brief Clears the current topology and keeps only an empty pixel domain.
     *
     * @param numPixels Number of pixels in the finite domain.
     */
    inline void initializeEmptyStorage(size_t numPixels) {
        rootNodeId_ = InvalidNode;
        numNodes_ = 0;
        smallestNodeMap_.assign(numPixels, InvalidNode);

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
        properPartCardinalityByNode_.clear();
        initializeProperPartStorage(numPixels);
        invalidateHigraNodeIdSpace();

        dfsIntervalCache_.entryIndex.clear();
        dfsIntervalCache_.exitIndex.clear();
        invalidateDfsIntervalCache();
        nodeSupportMetadataCache_.cardinalityByNode.clear();
        nodeSupportMetadataCache_.smallestPixelByNode.clear();
        nodeSupportMetadataCache_.invalidate();
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
        properHead_[slotId] = InvalidPixel;
        properTail_[slotId] = InvalidPixel;
        properPartCardinalityByNode_[slotId] = 0;
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
     * @param nodeId Dense internal node identifier.
     * @param context Operation context or diagnostic label.
     */
    inline void requireAliveNode(NodeId nodeId, const char* context) const {
        MMCFILTERS_CONTRACT_REQUIRE(isAlive(nodeId), throw std::invalid_argument(std::string(context) + " requires a live internal NodeId."));
    }

    /**
     * @brief Rejects invalid, released, or root node ids for non-root-only edits.
     *
     * @param nodeId Dense internal node identifier.
     * @param context Operation context or diagnostic label.
     */
    inline void requireAliveNonRootNode(NodeId nodeId, const char* context) const {
        requireAliveNode(nodeId, context);
        MMCFILTERS_CONTRACT_REQUIRE(!isRoot(nodeId), throw std::invalid_argument(std::string(context) + " cannot target the root node."));
    }

    /**
     * @brief Reconstructs direct proper-part linked lists from `smallestNodeMap_`.
     */
    inline void rebuildProperPartLinksFromSmallestNodeMap() {
        properHead_.assign(nodeParent_.size(), InvalidPixel);
        properTail_.assign(nodeParent_.size(), InvalidPixel);
        properPartCardinalityByNode_.assign(nodeParent_.size(), 0);
        initializeProperPartStorage(smallestNodeMap_.size());

        for (PixelId pixel = 0; pixel < static_cast<PixelId>(smallestNodeMap_.size()); ++pixel) {
            const NodeId smallestNodeSlotId = smallestNodeMap_[static_cast<std::size_t>(pixel)];
            if (smallestNodeSlotId == InvalidNode || smallestNodeSlotId >= static_cast<NodeId>(nodeParent_.size()) || isFreeSlot(smallestNodeSlotId)) {
                continue;
            }

            const PixelId tailProperPartId = properTail_[smallestNodeSlotId];
            if (tailProperPartId == InvalidPixel) {
                properHead_[smallestNodeSlotId] = pixel;
                properTail_[smallestNodeSlotId] = pixel;
            } else {
                nextProperPart_[tailProperPartId] = pixel;
                prevProperPart_[pixel] = tailProperPartId;
                properTail_[smallestNodeSlotId] = pixel;
            }
            properPartCardinalityByNode_[smallestNodeSlotId]++;
        }
    }

    /**
     * @brief Removes one pixel from the proper part of its current smallest node.
     *
     * @param smallestNodeSlotId Dense slot of the pixel's current smallest node.
     * @param pixel Pixel identifier.
     */
    inline void unlinkPixelFromProperPart(NodeId smallestNodeSlotId, PixelId pixel) noexcept {
        const PixelId prev = prevProperPart_[static_cast<size_t>(pixel)];
        const PixelId next = nextProperPart_[static_cast<size_t>(pixel)];

        if (prev == InvalidPixel) {
            properHead_[static_cast<size_t>(smallestNodeSlotId)] = next;
        } else {
            nextProperPart_[static_cast<size_t>(prev)] = next;
        }

        if (next == InvalidPixel) {
            properTail_[static_cast<size_t>(smallestNodeSlotId)] = prev;
        } else {
            prevProperPart_[static_cast<size_t>(next)] = prev;
        }

        nextProperPart_[static_cast<size_t>(pixel)] = InvalidPixel;
        prevProperPart_[static_cast<size_t>(pixel)] = InvalidPixel;
        --properPartCardinalityByNode_[static_cast<size_t>(smallestNodeSlotId)];
    }

    /**
     * @brief Appends one detached pixel to the target node's proper part.
     *
     * @param targetSlotId Destination.
     * @param pixel Pixel identifier.
     */
    inline void appendDetachedProperPart(NodeId targetSlotId, PixelId pixel) noexcept {
        smallestNodeMap_[static_cast<size_t>(pixel)] = targetSlotId;
        const PixelId tail = properTail_[static_cast<size_t>(targetSlotId)];
        if (tail == InvalidPixel) {
            properHead_[static_cast<size_t>(targetSlotId)] = pixel;
            properTail_[static_cast<size_t>(targetSlotId)] = pixel;
        } else {
            nextProperPart_[static_cast<size_t>(tail)] = pixel;
            prevProperPart_[static_cast<size_t>(pixel)] = tail;
            properTail_[static_cast<size_t>(targetSlotId)] = pixel;
        }
        ++properPartCardinalityByNode_[static_cast<size_t>(targetSlotId)];
    }

    /**
     * @brief Splices all direct proper parts from `sourceSlotId` to the tail of `targetSlotId`.
     *
     * @param targetSlotId Destination.
     * @param sourceSlotId Input.
     */
    inline void spliceProperPartsSlots(NodeId targetSlotId, NodeId sourceSlotId) noexcept {
        const PixelId sourceHead = properHead_[static_cast<size_t>(sourceSlotId)];
        if (sourceHead == InvalidPixel) {
            return;
        }

        for (PixelId pixel = sourceHead; pixel != InvalidPixel; pixel = nextProperPart_[static_cast<size_t>(pixel)]) {
            smallestNodeMap_[static_cast<size_t>(pixel)] = targetSlotId;
        }

        const PixelId sourceTail = properTail_[static_cast<size_t>(sourceSlotId)];
        const int sourceCount = properPartCardinalityByNode_[static_cast<size_t>(sourceSlotId)];
        const PixelId targetTail = properTail_[static_cast<size_t>(targetSlotId)];

        if (targetTail == InvalidPixel) {
            properHead_[static_cast<size_t>(targetSlotId)] = sourceHead;
            properTail_[static_cast<size_t>(targetSlotId)] = sourceTail;
        } else {
            nextProperPart_[static_cast<size_t>(targetTail)] = sourceHead;
            prevProperPart_[static_cast<size_t>(sourceHead)] = targetTail;
            properTail_[static_cast<size_t>(targetSlotId)] = sourceTail;
        }

        properPartCardinalityByNode_[static_cast<size_t>(targetSlotId)] += sourceCount;
        properHead_[static_cast<size_t>(sourceSlotId)] = InvalidPixel;
        properTail_[static_cast<size_t>(sourceSlotId)] = InvalidPixel;
        properPartCardinalityByNode_[static_cast<size_t>(sourceSlotId)] = 0;
    }

    /**
     * @brief Releases a live detached slot and invalidates dependent caches.
     *
     * @param slotNodeId Node identifier.
     */
    inline void releaseSlotNode(NodeId slotNodeId) {
        releaseSlotStorage(slotNodeId);
        numNodes_--;
        invalidateDfsIntervalCache();
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
        nodeSupportMetadataCache_.invalidate();
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
        invalidateDfsIntervalCache();
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
     * @brief Invalidates DFS entry/exit intervals used by topology queries.
     */
    inline void invalidateDfsIntervalCache() const noexcept { dfsIntervalCache_.invalidate(); }

    /**
     * @brief Recomputes interleaved zero-based DFS entry/exit event indices.
     */
    inline void recomputeDfsIntervalCache() const {
        dfsIntervalCache_.entryIndex.assign(nodeParent_.size(), -1);
        dfsIntervalCache_.exitIndex.assign(nodeParent_.size(), -1);

        if (rootNodeId_ == InvalidNode) {
            dfsIntervalCache_.valid = true;
            return;
        }

        int nextDfsEventIndex = 0;
        computeIncrementalAttributes(
            const_cast<MorphologicalTree*>(this), root(),
            [&](NodeId nodeId) -> void { dfsIntervalCache_.entryIndex[nodeId] = nextDfsEventIndex++; },
            [](NodeId, NodeId) -> void {},
            [&](NodeId nodeId) -> void { dfsIntervalCache_.exitIndex[nodeId] = nextDfsEventIndex++; });

        dfsIntervalCache_.valid = true;
    }

    /**
     * @brief Ensures that DFS entry/exit intervals are available.
     */
    inline void ensureDfsIntervalCache() const {
        if (!dfsIntervalCache_.valid) {
            recomputeDfsIntervalCache();
        }
    }

    /**
     * @brief Lazily builds the Euler-tour LCA cache on first use.
     *
     * @return Mutable reference to the updated object.
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
     * @param tree Tree topology.
     * @param rootNodeId Identifier of the traversal root.
     * @param preProcessing Callback invoked before visiting a node children.
     * @param mergeProcessing Callback invoked after completing one child.
     * @param postProcessing Callback invoked after all children are merged.
     */
    template <class PreProcessing, class MergeProcessing, class PostProcessing>
    static void computeIncrementalAttributes(MorphologicalTree* tree, NodeId rootNodeId, PreProcessing&& preProcessing, MergeProcessing&& mergeProcessing,
                                             PostProcessing&& postProcessing) {
        preProcessing(rootNodeId);
        for (NodeId childNodeId : tree->children(rootNodeId)) {
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
     * @param parentSlotId Parent-node value.
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
     * @param parentSlotId Parent-node value.
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
     * @param nodeId Dense internal node identifier.
     */
    inline void releaseNode(NodeId nodeId) {
        if (!isNode(nodeId) || !isAlive(nodeId) || isRoot(nodeId)) {
            return;
        }
        const NodeId nodeSlot = nodeId;
        if (nodeParent_[nodeSlot] != nodeSlot) {
            return;
        }
        if (numChildrenByNode_[nodeSlot] != 0 || properPartCardinalityByNode_[nodeSlot] != 0) {
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
        invalidateDfsIntervalCache();
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
        invalidateDfsIntervalCache();
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
            properPartCardinalityByNode_[static_cast<std::size_t>(childSlotId)] == 0) {
            freeNodeIds_.reserve(freeNodeIds_.size() + 1);
        }
        unlinkChildSlot(parentSlotId, childSlotId);
        nodeParent_[childSlotId] = childSlotId;
        invalidateDfsIntervalCache();
        if (releaseNodeFlag) {
            releaseNode(childId);
        }
    }

    /**
     * @brief Attaches `nodeId` as the last child of `parentNodeId`.
     *
     * @param parentNodeId Identifier of the parent node.
     * @param nodeId Dense internal node identifier.
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
     * @param nodeId Dense internal node identifier.
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
     * @param nodeId Dense internal node identifier.
     * @param newParentId Parent-node value.
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
     * @param sourceId Input.
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
     * @param targetNodeId Node identifier.
     * @param sourceNodeId Node identifier.
     * @param pixel Proper-part identifier.
     */
    inline void movePixelToProperPart(NodeId targetNodeId, NodeId sourceNodeId, PixelId pixel) {
        if (!isAlive(targetNodeId) || !isAlive(sourceNodeId) || !isPixel(pixel) || targetNodeId == sourceNodeId) {
            return;
        }
        const NodeId sourceSlotId = sourceNodeId;
        if (smallestNodeMap_[pixel] != sourceSlotId) {
            return;
        }
        unlinkPixelFromProperPart(sourceSlotId, pixel);
        appendDetachedProperPart(targetNodeId, pixel);
        bumpProperPartVersion();
    }

    /**
     * @brief Transfers all direct proper parts from `sourceNodeId` to `targetNodeId`.
     *
     * @param targetNodeId Node identifier.
     * @param sourceNodeId Node identifier.
     */
    inline void mergeProperParts(NodeId targetNodeId, NodeId sourceNodeId) {
        if (!isAlive(targetNodeId) || !isAlive(sourceNodeId) || targetNodeId == sourceNodeId) {
            return;
        }
        if (properHead_[static_cast<size_t>(sourceNodeId)] == InvalidPixel) {
            return;
        }
        spliceProperPartsSlots(targetNodeId, sourceNodeId);
        bumpProperPartVersion();
    }

    /**
     * @brief Promotes `nodeId` to be the connected root.
     *
     * @param nodeId Dense internal node identifier.
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
     * @brief Imports native node-parent and smallest-node-map buffers.
     *
     * This shared materialization path supports any builder that emits
     * independent node and pixel domains, including live nodes with an empty
     * proper part.
     *
     * @param nodeParent Parent-node value.
     * @param smallestNodeMap Pixel-indexed inclusion-smallest node map.
     * @param root Root node of the traversal.
     * @param gridDomain2D Optional regular-grid pixel domain.
     * @param semantics Hierarchy semantics validated by the operation.
     */
    void initializeNativeTopologyStorage(std::vector<NodeId> nodeParent, std::vector<NodeId> smallestNodeMap, NodeId root,
                                         std::optional<GridDomain2D> gridDomain2D, MorphologicalTreeSemantics semantics) {
        if (gridDomain2D && smallestNodeMap.size() != gridDomain2D->size("Native topology 2D pixel domain")) {
            throw std::invalid_argument("Native topology pixel domain must match the attached 2D grid.");
        }
        if (nodeParent.size() > static_cast<std::size_t>(std::numeric_limits<NodeId>::max())) {
            throw std::invalid_argument("Native topology internal-node domain exceeds NodeId range.");
        }
        if (smallestNodeMap.size() > static_cast<std::size_t>(std::numeric_limits<PixelId>::max())) {
            throw std::invalid_argument("Native topology pixel domain exceeds PixelId range.");
        }
        switch (semantics.kind) {
        case MorphologicalTreeKind::Generic:
        case MorphologicalTreeKind::MaxTree:
        case MorphologicalTreeKind::MinTree:
        case MorphologicalTreeKind::TreeOfShapes:
        case MorphologicalTreeKind::UnrestrictedResidualTree:
        case MorphologicalTreeKind::SaturatedResidualTree:
            break;
        default:
            throw std::invalid_argument("Native topology kind is not supported.");
        }
        const int numNodeSlots = static_cast<int>(nodeParent.size());
        const std::size_t numPixels = smallestNodeMap.size();
        if (numNodeSlots <= 0) {
            throw std::invalid_argument("Native topology import requires at least one internal node.");
        }
        if (numPixels == 0) {
            throw std::invalid_argument("Native topology import requires at least one pixel.");
        }
        if (root < 0 || root >= numNodeSlots) {
            throw std::invalid_argument("Native topology import requires a valid root node id.");
        }
        validateMorphologicalTreeSemantics(semantics);
        semantics_ = std::move(semantics);
        gridDomain2D_ = gridDomain2D;
        const auto requireMatchingGrid = [this](const RegularGridAdjacency2D& adjacency, const char* context) {
            if (!gridDomain2D_) {
                throw std::invalid_argument(std::string(context) + " requires an attached 2D grid domain.");
            }
            if (adjacency.getNumRows() != gridDomain2D_->rows || adjacency.getNumColumns() != gridDomain2D_->columns) {
                throw std::invalid_argument(std::string(context) + " must match the attached 2D grid.");
            }
        };
        if (const auto* context = std::get_if<SharedAdjacencyContext>(&semantics_.constructionContext)) {
            requireMatchingGrid(context->adjacency, "SharedAdjacencyContext adjacency");
        } else if (const auto* context = std::get_if<SaturatedResidualContext>(&semantics_.constructionContext)) {
            requireMatchingGrid(context->adjacency, "SaturatedResidualContext adjacency");
            if (context->infinityPixel < 0 || static_cast<std::size_t>(context->infinityPixel) >= gridDomain2D_->size("Saturated residual domain")) {
                throw std::invalid_argument("SaturatedResidualContext infinity pixel must belong to the attached 2D grid.");
            }
        } else if (const auto* convention = std::get_if<TopographicConvention>(&semantics_.constructionContext)) {
            if (convention->infinityPixel < 0) {
                throw std::invalid_argument("TopographicConvention infinity pixel must be non-negative.");
            }
            if (!gridDomain2D_) {
                throw std::invalid_argument("TopographicConvention requires an attached 2D grid domain.");
            }
            const std::int64_t extension = convention->domainExtension == TopographicDomainExtension::ExteriorRing ? 1 : -1;
            const std::int64_t activeRows = 2 * static_cast<std::int64_t>(gridDomain2D_->rows) + extension;
            const std::int64_t activeColumns = 2 * static_cast<std::int64_t>(gridDomain2D_->columns) + extension;
            if (activeRows <= 0 || activeColumns <= 0 || static_cast<std::int64_t>(convention->infinityPixel) >= activeRows * activeColumns) {
                throw std::invalid_argument("TopographicConvention infinity pixel must belong to the active topographic domain.");
            }
            if (const auto* immersion = std::get_if<ComplementaryGridImmersion>(&convention->immersion)) {
                requireMatchingGrid(immersion->complementaryAdjacencies.minAdjacency, "Topographic minimum adjacency");
                requireMatchingGrid(immersion->complementaryAdjacencies.maxAdjacency, "Topographic maximum adjacency");
            }
        }

        initializeEmptyStorage(numPixels);
        nodeParent_ = std::move(nodeParent);
        firstChild_.assign(static_cast<size_t>(numNodeSlots), InvalidNode);
        nextSibling_.assign(static_cast<size_t>(numNodeSlots), InvalidNode);
        prevSibling_.assign(static_cast<size_t>(numNodeSlots), InvalidNode);
        lastChild_.assign(static_cast<size_t>(numNodeSlots), InvalidNode);
        numChildrenByNode_.assign(static_cast<size_t>(numNodeSlots), 0);
        alive_.assign(static_cast<size_t>(numNodeSlots), 1);
        freeNodeIds_.clear();
        properHead_.assign(static_cast<size_t>(numNodeSlots), InvalidPixel);
        properTail_.assign(static_cast<size_t>(numNodeSlots), InvalidPixel);
        properPartCardinalityByNode_.assign(static_cast<size_t>(numNodeSlots), 0);
        smallestNodeMap_ = std::move(smallestNodeMap);
        initializeProperPartStorage(numPixels);

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

        rebuildProperPartLinksFromSmallestNodeMap();
        invalidateDfsIntervalCache();
        invalidateAllIterators();
        preservedExternalNodeIdOffset_.reset();
    }

    /**
     * @brief Copies and fully validates a native topology input.
     *
     * @param nodeParent Parent-node value.
     * @param smallestNodeMap Pixel-indexed inclusion-smallest node map.
     * @param root Root node of the traversal.
     * @param gridDomain2D Optional regular-grid pixel domain.
     * @param semantics Hierarchy semantics validated by the operation.
     */
    void initializeNativeTopology(std::span<const NodeId> nodeParent, std::span<const NodeId> smallestNodeMap, NodeId root,
                                  std::optional<GridDomain2D> gridDomain2D, MorphologicalTreeSemantics semantics) {
        if (nodeParent.size() > static_cast<std::size_t>(std::numeric_limits<NodeId>::max())) {
            throw std::invalid_argument("Native topology internal-node domain exceeds NodeId range.");
        }
        const int numNodeSlots = static_cast<int>(nodeParent.size());
        for (NodeId smallestNodeId : smallestNodeMap) {
            if (smallestNodeId < 0 || smallestNodeId >= numNodeSlots) {
                throw std::invalid_argument("Native topology import found a smallest-node-map entry outside the internal-node domain.");
            }
        }
        initializeNativeTopologyStorage(std::vector<NodeId>(nodeParent.begin(), nodeParent.end()),
                                        std::vector<NodeId>(smallestNodeMap.begin(), smallestNodeMap.end()), root, gridDomain2D, std::move(semantics));
        validateConnectedRootedTree();
    }

    /**
     * @brief Consumes producer-owned buffers paired with generic structural
     * evidence.
     *
     * @param nodeParent Parent-node value.
     * @param smallestNodeMap Pixel-indexed inclusion-smallest node map.
     * @param root Root node of the traversal.
     * @param gridDomain2D Optional regular-grid pixel domain.
     * @param semantics Hierarchy semantics validated by the operation.
     * @param topologyProof Proof that the native topology invariants hold.
     */
    void initializeValidatedNativeTopology(std::vector<NodeId>&& nodeParent, std::vector<NodeId>&& smallestNodeMap, NodeId root,
                                           std::optional<GridDomain2D> gridDomain2D, MorphologicalTreeSemantics semantics,
                                           detail::NativeTopologyProof&& topologyProof) {
        topologyProof.requireMatches(nodeParent.size(), smallestNodeMap.size(), root);
        initializeNativeTopologyStorage(std::move(nodeParent), std::move(smallestNodeMap), root, gridDomain2D, std::move(semantics));
#ifndef NDEBUG
        validateConnectedRootedTree();
#endif
    }

    /**
     * @brief Factory-only materialization of a proven owning native topology.
     *
     * @param nodeParent Parent-node value.
     * @param smallestNodeMap Pixel-indexed inclusion-smallest node map.
     * @param root Root node of the traversal.
     * @param gridDomain2D Optional regular-grid pixel domain.
     * @param semantics Hierarchy semantics validated by the operation.
     * @param topologyProof Proof that the native topology invariants hold.
     */
    MorphologicalTree(detail::MorphologicalTreeConstructionTag, std::vector<NodeId>&& nodeParent, std::vector<NodeId>&& smallestNodeMap, NodeId root,
                      std::optional<GridDomain2D> gridDomain2D, MorphologicalTreeSemantics semantics, detail::NativeTopologyProof&& topologyProof) {
        initializeValidatedNativeTopology(std::move(nodeParent), std::move(smallestNodeMap), root, gridDomain2D, std::move(semantics),
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
        other.semantics_ = MorphologicalTreeSemantics{};
        gridDomain2D_ = std::move(other.gridDomain2D_);
        other.gridDomain2D_.reset();
        numNodes_ = std::exchange(other.numNodes_, 0);
        preservedExternalNodeIdOffset_ = std::move(other.preservedExternalNodeIdOffset_);
        other.preservedExternalNodeIdOffset_.reset();
        editSessionOpen_ = false;
        editValidationStatistics_ = other.editValidationStatistics_;

        smallestNodeMap_ = std::move(other.smallestNodeMap_);
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
        properPartCardinalityByNode_ = std::move(other.properPartCardinalityByNode_);
        nextProperPart_ = std::move(other.nextProperPart_);
        prevProperPart_ = std::move(other.prevProperPart_);

        dfsIntervalCache_ = std::move(other.dfsIntervalCache_);
        other.dfsIntervalCache_ = {};
        lcaCache_.reset();
        other.lcaCache_.reset();
        nodeSupportMetadataCache_ = std::move(other.nodeSupportMetadataCache_);
        other.nodeSupportMetadataCache_ = {};

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
     * @brief Copying is disabled to keep topology storage explicit.
     */
    MorphologicalTree(const MorphologicalTree&) = delete;

    /**
     * @brief Copy assignment is disabled to keep topology storage explicit.
     *
     * @return Mutable reference to the updated object.
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
     * @return Mutable reference to the updated object.
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
     * @brief Tag-protected import from native mmcfilters topology buffers.
     *
     * This path is used by builders that already materialize internal-node
     * parent links and a row-major smallest-node map. Nodes are not required to
     * have a non-empty proper part, so the representation supports the general
     * connected-subset tree model.
     * The buffers are copied into canonical tree storage and validated as one
     * connected rooted hierarchy in which every node has non-empty subtree
     * support.
     *
     * @param nodeParent Parent-node value.
     * @param smallestNodeMap Pixel-indexed inclusion-smallest node map.
     * @param root Root node of the traversal.
     * @param rows Number of rows in the domain.
     * @param columns Number of columns in the domain.
     * @param semantics Hierarchy semantics validated by the operation.
     */
    MorphologicalTree(detail::MorphologicalTreeConstructionTag, std::span<const NodeId> nodeParent, std::span<const NodeId> smallestNodeMap, NodeId root,
                      int rows, int columns, MorphologicalTreeSemantics semantics) {
        initializeNativeTopology(nodeParent, smallestNodeMap, root, GridDomain2D{rows, columns}, std::move(semantics));
    }

    /**
     * @brief Imports a native hierarchy over an abstract finite pixel set.
     *
     * No row/column interpretation is attached. Regular-grid reconstruction
     * and geometry-dependent algorithms consequently reject this tree
     * explicitly, while purely topological and support-based algorithms remain
     * available.
     *
     * @param nodeParent Parent-node value.
     * @param smallestNodeMap Pixel-indexed inclusion-smallest node map.
     * @param root Root node of the traversal.
     * @param semantics Hierarchy semantics validated by the operation.
     */
    MorphologicalTree(detail::MorphologicalTreeConstructionTag, std::span<const NodeId> nodeParent, std::span<const NodeId> smallestNodeMap, NodeId root,
                      MorphologicalTreeSemantics semantics) {
        initializeNativeTopology(nodeParent, smallestNodeMap, root, std::nullopt, std::move(semantics));
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
        cloned.smallestNodeMap_ = smallestNodeMap_;
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
        cloned.properPartCardinalityByNode_ = properPartCardinalityByNode_;
        cloned.nextProperPart_ = nextProperPart_;
        cloned.prevProperPart_ = prevProperPart_;
        cloned.dfsIntervalCache_ = dfsIntervalCache_;
        cloned.lcaCache_.reset();
        cloned.nodeSupportMetadataCache_ = nodeSupportMetadataCache_;
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
    class ProperPartIterator;
    class ProperPartRange;
    class NodeSupportIterator;
    class NodeSupportRange;
    class PostOrderNodeIterator;
    class PostOrderNodeRange;
    class BreadthFirstNodeIterator;
    class BreadthFirstNodeRange;
    class AncestorNodeIterator;
    class AncestorNodeRange;
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
        MMCFILTERS_CONTRACT_REQUIRE(mutationVersion_ == expectedVersion,
                                    throw std::logic_error(std::string(context) + " cannot be used after the referenced tree topology has changed."));
    }

    /**
     * @brief Returns the size of the dense internal-node id domain.
     *
     * Some slots may currently be free and therefore not correspond to live nodes.
     *
     * @return The size of the dense internal-node id domain.
     */
    inline int numInternalNodeSlots() const { return static_cast<int>(nodeParent_.size()); }

    /**
     * @brief Returns the cardinality of the pixel domain.
     *
     * @return The number of pixels in the tree domain.
     */
    inline int numPixels() const { return static_cast<int>(smallestNodeMap_.size()); }

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
        return *preservedExternalNodeIdOffset_ + numInternalNodeSlots();
    }

    /**
     * @brief Returns the size of the requested node-id domain.
     *
     * `NodeIdSpace::Higra` means the preserved imported Higra domain, not the
     * compact domain that would be generated by exporting the current tree.
     *
     * @param outputSpace Node-id domain used to index the output.
     * @return The size of the requested node-id domain.
     */
    inline int getNodeIdSpaceSize(NodeIdSpace outputSpace) const {
        switch (outputSpace) {
        case NodeIdSpace::MorphologicalTree:
            return numInternalNodeSlots();
        case NodeIdSpace::Higra:
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
     * @param nodeId Dense internal node identifier.
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
    inline NodeId root() const { return rootNodeId_; }

    /**
     * @brief Tests whether `nodeId` belongs to the internal-node id domain.
     *
     * @param nodeId Dense internal node identifier.
     * @return True if nodeId belongs to the internal-node id domain; otherwise false.
     */
    inline bool isNode(NodeId nodeId) const noexcept { return nodeId >= 0 && nodeId < static_cast<int>(nodeParent_.size()); }

    /**
     * @brief Tests whether `pixel` belongs to the pixel domain.
     *
     * @param pixel Pixel identifier.
     * @return True if pixel belongs to the pixel domain; otherwise false.
     */
    inline bool isPixel(PixelId pixel) const noexcept { return pixel >= 0 && pixel < static_cast<int>(smallestNodeMap_.size()); }

    /**
     * @brief Tests whether a node slot currently represents a live node.
     *
     * @param nodeId Dense internal node identifier.
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
     * @param nodeId Dense internal node identifier.
     * @return True if nodeId is the current root; otherwise false.
     */
    inline bool isRoot(NodeId nodeId) const { return nodeId == root(); }

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
    inline int numLeafNodes() const {
        int count = 0;
        for (NodeId nodeId : aliveNodeIds()) {
            if (isLeaf(nodeId)) {
                ++count;
            }
        }
        return count;
    }

    /**
     * @brief Returns the number of direct children of `nodeId`.
     *
     * @param nodeId Dense internal node identifier.
     * @return The number of direct children of nodeId.
     */
    inline int numChildren(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::numChildren");
        return numChildrenByNode_[nodeId];
    }

    /**
     * @brief Returns the number of internal descendants of `nodeId`.
     *
     * @param nodeId Dense internal node identifier.
     * @return The number of internal descendants of nodeId.
     */
    inline int numDescendants(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::numDescendants");
        ensureDfsIntervalCache();
        return (dfsIntervalCache_.exitIndex[nodeId] - dfsIntervalCache_.entryIndex[nodeId] - 1) / 2;
    }

    /**
     * @brief Returns the number of siblings of `nodeId`.
     *
     * @param nodeId Dense internal node identifier.
     * @return The number of siblings of nodeId.
     */
    inline int numSiblings(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::numSiblings");
        if (isRoot(nodeId)) {
            return 0;
        }
        const NodeId parentNodeId = parent(nodeId);
        return std::max(0, numChildren(parentNodeId) - 1);
    }

    /**
     * @brief Returns the zero-based DFS entry-event index of `nodeId`.
     *
     * @param nodeId Dense internal node identifier.
     * @return The index at which DFS enters nodeId in the interleaved event sequence.
     */
    inline int dfsEntryIndex(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::dfsEntryIndex");
        ensureDfsIntervalCache();
        return dfsIntervalCache_.entryIndex[nodeId];
    }

    /**
     * @brief Returns the zero-based DFS exit-event index of `nodeId`.
     *
     * @param nodeId Dense internal node identifier.
     * @return The index at which DFS exits nodeId in the interleaved event sequence.
     */
    inline int dfsExitIndex(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::dfsExitIndex");
        ensureDfsIntervalCache();
        return dfsIntervalCache_.exitIndex[nodeId];
    }

    /**
     * @brief Returns the first direct child of `nodeId`, or `InvalidNode`.
     *
     * @param nodeId Dense internal node identifier.
     * @return The first direct child of nodeId, or InvalidNode.
     */
    inline NodeId getFirstChild(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getFirstChild");
        return firstChild_[nodeId];
    }

    /**
     * @brief Returns the next sibling of `nodeId`, or `InvalidNode`.
     *
     * @param nodeId Dense internal node identifier.
     * @return The next sibling of nodeId, or InvalidNode.
     */
    inline NodeId getNextSibling(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::getNextSibling");
        return nextSibling_[nodeId];
    }

    /**
     * @brief Tests whether `nodeId` has no direct children.
     *
     * @param nodeId Dense internal node identifier.
     * @return True if nodeId has no direct children; otherwise false.
     */
    inline bool isLeaf(NodeId nodeId) const { return getFirstChild(nodeId) == InvalidNode; }

    /**
     * @brief Returns the cardinality of the proper part of `nodeId`.
     *
     * @param nodeId Dense internal node identifier.
     * @return The number of pixels in the proper part of nodeId.
     */
    inline int properPartCardinality(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::properPartCardinality");
        return properPartCardinalityByNode_[nodeId];
    }

    /**
     * @brief Tests whether `nodeId` has an empty proper part.
     *
     * @param nodeId Dense internal node identifier.
     * @return True if the proper part of nodeId is empty; otherwise false.
     */
    inline bool hasEmptyProperPart(NodeId nodeId) const { return properPartCardinality(nodeId) == 0; }

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
        return parent(childId) == parentNodeId;
    }

    /**
     * @brief Returns the direct parent of `nodeId`.
     *
     * The root and detached nodes report themselves as parent.
     *
     * @param nodeId Dense internal node identifier.
     * @return The direct parent of nodeId.
     */
    inline NodeId parent(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::parent");
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
     * @brief Returns the inclusion-smallest node containing `pixel`.
     *
     * @param pixel Pixel in the finite tree domain.
     * @return The smallest node whose support contains pixel, or InvalidNode for an invalid pixel in checked builds.
     */
    inline NodeId smallestNode(PixelId pixel) const {
        if constexpr (contract::validationsEnabled) {
            return isPixel(pixel) ? smallestNodeMap_[static_cast<size_t>(pixel)] : InvalidNode;
        }
        return smallestNodeMap_[static_cast<size_t>(pixel)];
    }

    /**
     * @brief Returns the pixel-indexed smallest-node map.
     *
     * @return A read-only view `P` satisfying `pixel` in `properPart(P[pixel])`.
     */
    [[nodiscard]] inline std::span<const NodeId> smallestNodeMap() const noexcept { return smallestNodeMap_; }

    /**
     * @brief Returns all live leaf nodes in the current hierarchy.
     *
     * @return All live leaf nodes in the current hierarchy.
     */
    std::vector<NodeId> leaves() const {
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
     * Algorithms must use explicit capabilities such as `nodeAltitudeOrder()`
     * and the typed construction context instead of dispatching on this label.
     *
     * @return The optional descriptive hierarchy-family label.
     */
    inline MorphologicalTreeKind kind() const noexcept { return semantics_.kind; }

    /**
     * @brief Returns all generic semantic capabilities of this hierarchy.
     *
     * @return All generic semantic capabilities of this hierarchy.
     */
    inline const MorphologicalTreeSemantics& semantics() const noexcept { return semantics_; }

    /**
     * @brief Returns the global parent-to-child altitude ordering constraint.
     *
     * @return The global parent-to-child altitude ordering constraint.
     */
    inline NodeAltitudeOrder nodeAltitudeOrder() const noexcept { return semantics_.nodeAltitudeOrder; }

    /**
     * @brief Returns the typed construction context retained by this tree.
     *
     * @return The typed construction context retained by this tree.
     */
    inline const MorphologicalTreeConstructionContext& constructionContext() const noexcept { return semantics_.constructionContext; }

    /** @brief Returns the shared-adjacency context, or `nullptr`. @return Stored context when its variant is active. */
    inline const SharedAdjacencyContext* sharedAdjacencyContext() const noexcept {
        return std::get_if<SharedAdjacencyContext>(&semantics_.constructionContext);
    }

    /** @brief Returns the saturated-residual context, or `nullptr`. @return Stored context when its variant is active. */
    inline const SaturatedResidualContext* saturatedResidualContext() const noexcept {
        return std::get_if<SaturatedResidualContext>(&semantics_.constructionContext);
    }

    /** @brief Returns the topographic convention, or `nullptr`. @return Stored convention when its variant is active. */
    inline const TopographicConvention* topographicConvention() const noexcept {
        return std::get_if<TopographicConvention>(&semantics_.constructionContext);
    }

    /**
     * @brief Returns the number of currently live nodes.
     *
     * @return The number of currently live nodes.
     */
    inline int numNodes() const noexcept { return numNodes_; }

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
        MMCFILTERS_CONTRACT_REQUIRE(!isEditing(),
                                    throw std::logic_error(std::string(context) + " requires a committed MorphologicalTree; an edit session is still open."));
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
     * @brief Tests whether pixel ids have an attached row/column layout.
     *
     * @return True if pixel ids have an attached row/column layout; otherwise false.
     */
    [[nodiscard]] inline bool hasGridDomain2D() const noexcept { return gridDomain2D_.has_value(); }

    /**
     * @brief Returns the optional regular 2D pixel domain.
     *
     * @return The optional regular 2D pixel domain.
     */
    [[nodiscard]] inline const std::optional<GridDomain2D>& gridDomain2D() const noexcept { return gridDomain2D_; }

    /**
     * @brief Returns the regular 2D domain or rejects a geometry-dependent call.
     *
     * @param context Operation context or diagnostic label.
     * @return The regular 2D domain or rejects a geometry-dependent call.
     */
    inline const GridDomain2D& requireGridDomain2D(const char* context) const {
        if (!gridDomain2D_) {
            throw std::invalid_argument(std::string(context) + " requires a regular 2D pixel domain.");
        }
        return *gridDomain2D_;
    }

    /**
     * @brief Returns the number of rows in the regular 2D pixel domain.
     *
     * @return The number of rows in the regular 2D pixel domain.
     */
    inline int numRows() const { return requireGridDomain2D("MorphologicalTree::numRows").rows; }

    /**
     * @brief Returns the number of columns in the regular 2D pixel domain.
     *
     * @return The number of columns in the regular 2D pixel domain.
     */
    inline int numColumns() const { return requireGridDomain2D("MorphologicalTree::numColumns").columns; }

    /**
     * @brief Validates one rooted tree with a complete smallest-node map and non-empty supports.
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

        const int numSlots = numInternalNodeSlots();
        const int domainPixelCount = numPixels();
        std::vector<int> expectedChildrenByNode(static_cast<size_t>(numSlots), 0);
        std::vector<int> expectedProperPartCardinalityByNode(static_cast<size_t>(numSlots), 0);
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
            throw std::runtime_error("Connected-tree validation found numNodes() out of sync with the alive node slots.");
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

        for (PixelId pixel = 0; pixel < domainPixelCount; ++pixel) {
            const NodeId smallestNodeId = smallestNodeMap_[static_cast<size_t>(pixel)];
            if (!isAlive(smallestNodeId)) {
                throw std::runtime_error("Connected-tree validation found a pixel without a live smallest node.");
            }
            expectedProperPartCardinalityByNode[static_cast<size_t>(smallestNodeId)] += 1;
        }

        std::vector<uint8_t> seenProperPart(static_cast<size_t>(domainPixelCount), 0);
        for (NodeId nodeId = 0; nodeId < numSlots; ++nodeId) {
            if (!isAlive(nodeId)) {
                continue;
            }

            int actualProperPartCardinality = 0;
            PixelId previousProperPartId = InvalidPixel;
            for (PixelId pixel = properHead_[static_cast<size_t>(nodeId)]; pixel != InvalidPixel;
                 pixel = nextProperPart_[static_cast<size_t>(pixel)]) {
                if (!isPixel(pixel)) {
                    throw std::runtime_error("Connected-tree validation found an invalid pixel id in a proper-part list.");
                }
                if (smallestNodeMap_[static_cast<size_t>(pixel)] != nodeId) {
                    throw std::runtime_error("Connected-tree validation found a proper-part list that disagrees with the smallest-node map.");
                }
                if (prevProperPart_[static_cast<size_t>(pixel)] != previousProperPartId) {
                    throw std::runtime_error("Connected-tree validation found broken previous-proper-part links.");
                }
                if (seenProperPart[static_cast<size_t>(pixel)] != 0) {
                    throw std::runtime_error("Connected-tree validation found a proper part referenced multiple times.");
                }
                seenProperPart[static_cast<size_t>(pixel)] = 1;
                previousProperPartId = pixel;
                ++actualProperPartCardinality;
                if (actualProperPartCardinality > domainPixelCount) {
                    throw std::runtime_error("Connected-tree validation found a cycle in a proper-part list.");
                }
            }

            if (properHead_[static_cast<size_t>(nodeId)] == InvalidPixel) {
                if (properTail_[static_cast<size_t>(nodeId)] != InvalidPixel) {
                    throw std::runtime_error("Connected-tree validation found an empty proper-part list with a non-empty tail pointer.");
                }
            } else {
                if (properTail_[static_cast<size_t>(nodeId)] != previousProperPartId) {
                    throw std::runtime_error("Connected-tree validation found an incorrect proper-part tail pointer.");
                }
                if (nextProperPart_[static_cast<size_t>(previousProperPartId)] != InvalidPixel) {
                    throw std::runtime_error("Connected-tree validation found a proper-part list whose tail still points forward.");
                }
            }

            if (actualProperPartCardinality != properPartCardinalityByNode_[static_cast<size_t>(nodeId)]) {
                throw std::runtime_error("Connected-tree validation found an incorrect direct proper-part count cache.");
            }
            if (actualProperPartCardinality != expectedProperPartCardinalityByNode[static_cast<size_t>(nodeId)]) {
                throw std::runtime_error("Connected-tree validation found proper-part lists out of sync with the smallest-node map.");
            }
        }

        for (PixelId pixel = 0; pixel < domainPixelCount; ++pixel) {
            if (seenProperPart[static_cast<size_t>(pixel)] == 0) {
                throw std::runtime_error("Connected-tree validation found a pixel missing from the proper-part lists.");
            }
        }

        std::vector<int> subtreeSupport = expectedProperPartCardinalityByNode;
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
     * @brief Tests whether every live node has a non-empty proper part.
     *
     * @return True exactly when the connected-subset tree is a morphological
     * tree of partial partitions.
     */
    [[nodiscard]] bool isTreeOfPartialPartitions() const {
        if (numNodes_ <= 0) {
            return false;
        }
        for (NodeId nodeId = 0; nodeId < numInternalNodeSlots(); ++nodeId) {
            if (isAlive(nodeId) && properPartCardinalityByNode_[static_cast<std::size_t>(nodeId)] == 0) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Validates the connected-subset invariants and the non-empty proper-part specialization.
     *
     * @throws std::runtime_error if the tree is invalid or any live node has an empty proper part.
     */
    void validateTreeOfPartialPartitions() const {
        validateConnectedRootedTree();
        for (NodeId nodeId = 0; nodeId < numInternalNodeSlots(); ++nodeId) {
            if (isAlive(nodeId) && properPartCardinalityByNode_[static_cast<std::size_t>(nodeId)] == 0) {
                throw std::runtime_error("Tree-of-partial-partitions validation found a node with an empty proper part: " + std::to_string(nodeId) + ".");
            }
        }
    }

    /**
     * @brief Tests whether `u` is an ancestor of `v`.
     *
     * @param u First endpoint or element.
     * @param v Second endpoint or element.
     * @return True if u is an ancestor of v; otherwise false.
     */
    inline bool isAncestor(NodeId u, NodeId v) const {
        requireAliveNode(u, "MorphologicalTree::isAncestor");
        requireAliveNode(v, "MorphologicalTree::isAncestor");
        ensureDfsIntervalCache();
        const NodeId slotU = u;
        const NodeId slotV = v;
        return dfsIntervalCache_.entryIndex[slotU] <= dfsIntervalCache_.entryIndex[slotV] &&
               dfsIntervalCache_.exitIndex[slotU] >= dfsIntervalCache_.exitIndex[slotV];
    }

    /**
     * @brief Tests whether `u` is a descendant of `v`.
     *
     * @param u First endpoint or element.
     * @param v Second endpoint or element.
     * @return True if u is a descendant of v; otherwise false.
     */
    inline bool isDescendant(NodeId u, NodeId v) const {
        requireAliveNode(u, "MorphologicalTree::isDescendant");
        requireAliveNode(v, "MorphologicalTree::isDescendant");
        ensureDfsIntervalCache();
        const NodeId slotU = u;
        const NodeId slotV = v;
        return dfsIntervalCache_.entryIndex[slotV] <= dfsIntervalCache_.entryIndex[slotU] &&
               dfsIntervalCache_.exitIndex[slotV] >= dfsIntervalCache_.exitIndex[slotU];
    }
    /**
     * @brief Tests whether `u` and `v` are comparable in the ancestry order.
     *
     * @param u First endpoint or element.
     * @param v Second endpoint or element.
     * @return True if u and v are comparable in the ancestry order; otherwise false.
     */
    inline bool isComparable(NodeId u, NodeId v) const { return isAncestor(u, v) || isAncestor(v, u); }

    /**
     * @brief Tests whether `u` is a strict ancestor of `v`.
     *
     * @param u First endpoint or element.
     * @param v Second endpoint or element.
     * @return True if u is a strict ancestor of v; otherwise false.
     */
    inline bool isStrictAncestor(NodeId u, NodeId v) const { return u != v && isAncestor(u, v); }

    /**
     * @brief Tests whether `u` is a strict descendant of `v`.
     *
     * @param u First endpoint or element.
     * @param v Second endpoint or element.
     * @return True if u is a strict descendant of v; otherwise false.
     */
    inline bool isStrictDescendant(NodeId u, NodeId v) const { return u != v && isDescendant(u, v); }

    /**
     * @brief Tests whether `u` and `v` are strictly comparable in the ancestry order.
     *
     * @param u First endpoint or element.
     * @param v Second endpoint or element.
     * @return True if u and v are strictly comparable in the ancestry order; otherwise false.
     */
    inline bool isStrictComparable(NodeId u, NodeId v) const { return isStrictAncestor(u, v) || isStrictAncestor(v, u); }

    /**
     * @brief Returns the lowest common ancestor of `u` and `v`.
     *
     * @param u First endpoint or element.
     * @param v Second endpoint or element.
     * @return The lowest common ancestor of u and v.
     */
    inline NodeId lowestCommonAncestor(NodeId u, NodeId v) const {
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
        if (parent(u) == u || parent(v) == v) {
            return InvalidNode;
        }
        return ensureLcaCache().findLowestCommonAncestor(u, v);
    }

    /**
     * @brief Returns a fail-fast range over all live node ids.
     *
     * @return A fail-fast range over all live node ids.
     */
    inline AliveNodeRange aliveNodeIds() const { return AliveNodeRange(this, 0, numInternalNodeSlots(), nodeStructureVersion_); }

    /**
     * @brief Returns a fail-fast range over the direct children of `nodeId`.
     *
     * @param nodeId Dense internal node identifier.
     * @return A fail-fast range over the direct children of nodeId.
     */
    inline ChildrenRange children(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::children");
        return ChildrenRange(this, firstChild_[nodeId], topologyVersion_);
    }

    /**
     * @brief Returns a fail-fast range over the pixels in the proper part of `nodeId`.
     *
     * @param nodeId Dense internal node identifier.
     * @return A fail-fast range over the pixels in the proper part of nodeId.
     */
    inline ProperPartRange properPart(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::properPart");
        return ProperPartRange(this, properHead_[nodeId], properPartVersion_);
    }

    /**
     * @brief Returns a fail-fast range over the pixels in the support of `nodeId`.
     *
     * The range walks the subtree rooted at `nodeId` and yields every pixel in
     * the proper parts of those nodes, without materialising a vector.
     *
     * @param nodeId Dense internal node identifier.
     * @return A fail-fast range over the pixels in the support represented by nodeId.
     */
    inline NodeSupportRange nodeSupport(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::nodeSupport");
        return NodeSupportRange(this, nodeId, topologyVersion_, properPartVersion_);
    }

    /**
     * @brief Returns a post-order traversal range rooted at the connected root.
     *
     * @return A post-order traversal range rooted at the connected root.
     */
    inline PostOrderNodeRange postOrder() const { return PostOrderNodeRange(this, root(), topologyVersion_); }

    /**
     * @brief Returns a post-order traversal range rooted at `rootNodeId`.
     *
     * @param rootNodeId Identifier of the traversal root.
     * @return A post-order traversal range rooted at rootNodeId.
     */
    inline PostOrderNodeRange postOrder(NodeId rootNodeId) const {
        requireAliveNode(rootNodeId, "MorphologicalTree::postOrder");
        return PostOrderNodeRange(this, rootNodeId, topologyVersion_);
    }

    /**
     * @brief Returns a breadth-first traversal range rooted at the connected root.
     *
     * @return A breadth-first traversal range rooted at the connected root.
     */
    inline BreadthFirstNodeRange breadthFirstTraversal() const { return BreadthFirstNodeRange(this, root(), topologyVersion_); }

    /**
     * @brief Returns a breadth-first traversal range rooted at `rootNodeId`.
     *
     * @param rootNodeId Identifier of the traversal root.
     * @return A breadth-first traversal range rooted at rootNodeId.
     */
    inline BreadthFirstNodeRange breadthFirstTraversal(NodeId rootNodeId) const {
        requireAliveNode(rootNodeId, "MorphologicalTree::breadthFirstTraversal");
        return BreadthFirstNodeRange(this, rootNodeId, topologyVersion_);
    }

    /**
     * @brief Returns `nodeId` and its proper ancestors through the connected root.
     *
     * The queried node is the first element, so this range is inclusive of
     * self. The connected root is the last element.
     *
     * @param nodeId Dense internal node identifier.
     * @return The inclusive ancestor chain from nodeId to the connected root.
     */
    inline AncestorNodeRange ancestors(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::ancestors");
        return AncestorNodeRange(this, nodeId, topologyVersion_);
    }

    /**
     * @brief Returns the path that connects `sourceNodeId` and `targetNodeId`.
     *
     * @param sourceNodeId Node identifier.
     * @param targetNodeId Node identifier.
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
     * @param nodeId Dense internal node identifier.
     * @return A pre-order traversal range over the subtree of nodeId.
     */
    inline SubtreeNodeRange subtreeNodes(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::subtreeNodes");
        return SubtreeNodeRange(this, nodeId, topologyVersion_);
    }

    /**
     * @brief Returns a range over all proper descendants of `nodeId`.
     *
     * @param nodeId Dense internal node identifier.
     * @return A range over all proper descendants of nodeId.
     */
    inline DescendantNodeRange descendants(NodeId nodeId) const {
        requireAliveNode(nodeId, "MorphologicalTree::descendants");
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
     * @param nodeId Dense internal node identifier.
     */
    inline void pruneNode(NodeId nodeId) {
        requireNotEditing("MorphologicalTree::pruneNode");
        requireAliveNonRootNode(nodeId, "MorphologicalTree::pruneNode");
        const NodeId parentNodeId = parent(nodeId);
        if (parentNodeId == InvalidNode || parentNodeId == nodeId) {
            throw std::invalid_argument("MorphologicalTree::pruneNode requires an attached non-root node.");
        }
        const NodeId parentSlotId = parentNodeId;
        std::size_t releaseCount = 0;
        for ([[maybe_unused]] NodeId subtreeNode : subtreeNodes(nodeId)) {
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
            mergeProperParts(parentSlotId, currentId);
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
     * @param nodeId Dense internal node identifier.
     */
    inline void mergeNodeIntoParent(NodeId nodeId) {
        requireNotEditing("MorphologicalTree::mergeNodeIntoParent");
        requireAliveNonRootNode(nodeId, "MorphologicalTree::mergeNodeIntoParent");
        const NodeId parentNodeId = parent(nodeId);
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
        /** @brief Dense node identifier of the euler. */
        std::vector<NodeId> euler_;
        /** @brief Depth buffer. */
        std::vector<int> depth_;
        /** @brief First occurrence buffer. */
        std::vector<int> firstOccurrence_;
        /** @brief Log2 buffer. */
        std::vector<int> log2_;
        /** @brief Sparse table buffer. */
        std::vector<int> sparseTable_;
        /** @brief Sparse table stride. */
        int sparseTableStride_ = 0;
        /** @brief References the tree used by the component. */
        const MorphologicalTree* tree_ = nullptr;

        /**
         * @brief Appends one DFS step to the Euler tour and recurses into children.
         *
         * @param nodeId Dense internal node identifier.
         * @param currentDepth Depth of the current traversal node.
         */
        void depthFirstTraversal(NodeId nodeId, int currentDepth) {
            const NodeId slotId = nodeId;
            if (firstOccurrence_[static_cast<size_t>(slotId)] == -1) {
                firstOccurrence_[static_cast<size_t>(slotId)] = static_cast<int>(euler_.size());
            }

            euler_.push_back(nodeId);
            depth_.push_back(currentDepth);

            for (NodeId childNodeId : tree_->children(nodeId)) {
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
         * @param column Zero-based column coordinate.
         * @return The mapped sparse-table row/column pair to the flat storage index.
         */
        std::size_t sparseTableIndex(int row, int column) const {
            return static_cast<std::size_t>(row) * static_cast<std::size_t>(sparseTableStride_) + static_cast<std::size_t>(column);
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
         * @param tree Tree topology.
         */
        explicit LCAEulerRMQ(const MorphologicalTree* tree) : tree_(tree) {
            if (tree_ == nullptr || tree_->root() == InvalidNode) {
                return;
            }

            euler_.reserve(static_cast<size_t>(std::max(0, 2 * tree_->numNodes() - 1)));
            depth_.reserve(static_cast<size_t>(std::max(0, 2 * tree_->numNodes() - 1)));
            firstOccurrence_.assign(static_cast<size_t>(tree_->numInternalNodeSlots()), -1);

            depthFirstTraversal(tree_->root(), 0);
            buildSparseTable();
        }

        /**
         * @brief Returns the lowest common ancestor of `u` and `v`.
         *
         * @param u First endpoint or element.
         * @param v Second endpoint or element.
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
        /** @brief Dense node identifier of the current. */
        NodeId current_ = InvalidNode;
        /** @brief Dense node identifier of the end. */
        NodeId end_ = InvalidNode;
        /** @brief Expected version used to detect stale derived state. */
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
         * @param tree Tree topology.
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
         * @return Mutable reference to the updated object.
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
        /** @brief Dense node identifier of the begin. */
        NodeId begin_ = InvalidNode;
        /** @brief Dense node identifier of the end. */
        NodeId end_ = InvalidNode;
        /** @brief Expected version used to detect stale derived state. */
        std::size_t expectedVersion_ = 0;

      public:
        /**
         * @brief Creates an empty live-node range.
         */
        AliveNodeRange() = default;

        /**
         * @brief Creates a fail-fast live-node range over dense slots.
         *
         * @param tree Tree topology.
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
        /** @brief Dense node identifier of the current local. */
        NodeId currentLocal_ = InvalidNode;
        /** @brief Expected version used to detect stale derived state. */
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
         * @param tree Tree topology.
         * @param currentLocal Current local iterator position.
         * @param expectedVersion Mutation version required by the operation.
         */
        ChildrenIterator(const MorphologicalTree* tree, NodeId currentLocal, std::size_t expectedVersion)
            : T_(tree), currentLocal_(currentLocal), expectedVersion_(expectedVersion) {}

        /**
         * @brief Advances to the next sibling in the child list.
         *
         * @return Mutable reference to the updated object.
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
        /** @brief Dense node identifier of the first local. */
        NodeId firstLocal_ = InvalidNode;
        /** @brief Expected version used to detect stale derived state. */
        std::size_t expectedVersion_ = 0;

      public:
        /**
         * @brief Creates an empty child range.
         */
        ChildrenRange() = default;

        /**
         * @brief Creates a range over a linked list of direct children.
         *
         * @param tree Tree topology.
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
     * @brief Iterator over the pixels in one node's proper part.
     */
    class ProperPartIterator {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Pixel identifier of the current proper part. */
        PixelId currentProperPart_ = InvalidPixel;
        /** @brief Expected version used to detect stale derived state. */
        std::size_t expectedVersion_ = 0;

      public:
        /// Standard iterator category exposed for STL compatibility.
        using iterator_category = std::forward_iterator_tag;
        /// Proper-part id yielded by the iterator.
        using value_type = PixelId;
        /// Signed distance type exposed for STL compatibility.
        using difference_type = std::ptrdiff_t;
        /// Pointer type exposed for STL compatibility.
        using pointer = const PixelId*;
        /// Reference type exposed for STL compatibility.
        using reference = const PixelId&;

        /**
         * @brief Creates the end/sentinel proper-part iterator.
         */
        ProperPartIterator() = default;

        /**
         * @brief Creates a proper-part iterator starting from a linked-list entry.
         *
         * @param tree Tree topology.
         * @param currentProperPart Proper-part data.
         * @param expectedVersion Mutation version required by the operation.
         */
        ProperPartIterator(const MorphologicalTree* tree, PixelId currentProperPart, std::size_t expectedVersion)
            : T_(tree), currentProperPart_(currentProperPart), expectedVersion_(expectedVersion) {}

        /**
         * @brief Advances to the next pixel in the proper part.
         *
         * @return Mutable reference to the updated object.
         */
        ProperPartIterator& operator++() {
            T_->checkProperPartIteratorVersion(expectedVersion_);
            if (T_ && currentProperPart_ != InvalidPixel) {
                currentProperPart_ = T_->nextProperPart_[currentProperPart_];
            }
            return *this;
        }

        /**
         * @brief Returns the current proper-part id.
         *
         * @return The current proper-part id.
         */
        PixelId operator*() const {
            T_->checkProperPartIteratorVersion(expectedVersion_);
            return currentProperPart_;
        }

        /**
         * @brief Compares proper-part iterator positions.
         *
         * @param other Object to compare with or transfer from.
         * @return True when the documented condition holds; otherwise false.
         */
        bool operator==(const ProperPartIterator& other) const { return currentProperPart_ == other.currentProperPart_; }

        /**
         * @brief Compares proper-part iterator positions for inequality.
         *
         * @param other Object to compare with or transfer from.
         * @return True when the documented condition holds; otherwise false.
         */
        bool operator!=(const ProperPartIterator& other) const { return !(*this == other); }
    };

    /**
     * @brief Range wrapper for direct proper-part iteration.
     */
    class ProperPartRange {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Pixel identifier of the first proper part. */
        PixelId firstProperPart_ = InvalidPixel;
        /** @brief Expected version used to detect stale derived state. */
        std::size_t expectedVersion_ = 0;

      public:
        /**
         * @brief Creates an empty direct proper-part range.
         */
        ProperPartRange() = default;

        /**
         * @brief Creates a range over one node's direct proper-part list.
         *
         * @param tree Tree topology.
         * @param firstProperPart Proper-part data.
         * @param expectedVersion Mutation version required by the operation.
         */
        ProperPartRange(const MorphologicalTree* tree, PixelId firstProperPart, std::size_t expectedVersion)
            : T_(tree), firstProperPart_(firstProperPart), expectedVersion_(expectedVersion) {}

        /**
         * @brief Returns an iterator at the first direct proper part.
         *
         * @return An iterator at the first direct proper part.
         */
        ProperPartIterator begin() const { return ProperPartIterator(T_, firstProperPart_, expectedVersion_); }

        /**
         * @brief Returns the direct proper-part range sentinel.
         *
         * @return The direct proper-part range sentinel.
         */
        ProperPartIterator end() const { return ProperPartIterator(T_, InvalidPixel, expectedVersion_); }
    };

    /**
     * @brief Iterator over all pixels in one node support.
     */
    class NodeSupportIterator {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Dense node identifier of the node stack. */
        std::vector<NodeId> nodeStack_;
        /** @brief Pixel identifier of the current proper part. */
        PixelId currentProperPart_ = InvalidPixel;
        /** @brief Expected topology version used to detect stale derived state. */
        std::size_t expectedTopologyVersion_ = 0;
        /** @brief Expected proper part version used to detect stale derived state. */
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
         * @param nodeId Dense internal node identifier.
         */
        void pushChildren(NodeId nodeId) {
            std::vector<NodeId> children;
            for (NodeId childId : T_->children(nodeId)) {
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
            while (currentProperPart_ == InvalidPixel && !nodeStack_.empty()) {
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
        using value_type = PixelId;
        /// Signed distance type exposed for STL compatibility.
        using difference_type = std::ptrdiff_t;
        /// Pointer type exposed for STL compatibility.
        using pointer = const PixelId*;
        /// Reference type exposed for STL compatibility.
        using reference = const PixelId&;

        /**
         * @brief Creates the end/sentinel node-support iterator.
         */
        NodeSupportIterator() = default;

        /**
         * @brief Creates an iterator over every proper part in a rooted subtree.
         *
         * @param tree Tree topology.
         * @param rootNodeId Identifier of the traversal root.
         * @param expectedTopologyVersion Topology version captured by the iterator.
         * @param expectedProperPartVersion Proper-part data.
         */
        NodeSupportIterator(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedTopologyVersion, std::size_t expectedProperPartVersion)
            : T_(tree), expectedTopologyVersion_(expectedTopologyVersion), expectedProperPartVersion_(expectedProperPartVersion) {
            if (T_ && rootNodeId != InvalidNode) {
                nodeStack_.push_back(rootNodeId);
                settle();
            }
        }

        /**
         * @brief Advances to the next pixel in the node support.
         *
         * @return Mutable reference to the updated object.
         */
        NodeSupportIterator& operator++() {
            checkVersions();
            if (currentProperPart_ != InvalidPixel) {
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
        PixelId operator*() const {
            checkVersions();
            return currentProperPart_;
        }

        /**
         * @brief Compares node-support iterator positions.
         *
         * @param other Object to compare with or transfer from.
         * @return True when the documented condition holds; otherwise false.
         */
        bool operator==(const NodeSupportIterator& other) const { return currentProperPart_ == other.currentProperPart_; }

        /**
         * @brief Compares node-support iterator positions for inequality.
         *
         * @param other Object to compare with or transfer from.
         * @return True when the documented condition holds; otherwise false.
         */
        bool operator!=(const NodeSupportIterator& other) const { return !(*this == other); }
    };

    /**
     * @brief Range wrapper for node-support pixel iteration.
     */
    class NodeSupportRange {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Dense node identifier of the root node identifier. */
        NodeId rootNodeId_ = InvalidNode;
        /** @brief Expected topology version used to detect stale derived state. */
        std::size_t expectedTopologyVersion_ = 0;
        /** @brief Expected proper part version used to detect stale derived state. */
        std::size_t expectedProperPartVersion_ = 0;

      public:
        /**
         * @brief Creates an empty node-support range.
         */
        NodeSupportRange() = default;

        /**
         * @brief Creates a range over all proper parts in the subtree of `rootNodeId`.
         *
         * @param tree Tree topology.
         * @param rootNodeId Identifier of the traversal root.
         * @param expectedTopologyVersion Topology version captured by the iterator.
         * @param expectedProperPartVersion Proper-part data.
         */
        NodeSupportRange(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedTopologyVersion, std::size_t expectedProperPartVersion)
            : T_(tree), rootNodeId_(rootNodeId), expectedTopologyVersion_(expectedTopologyVersion), expectedProperPartVersion_(expectedProperPartVersion) {}

        /**
         * @brief Returns an iterator at the first pixel in the node support.
         *
         * @return An iterator at the first pixel in the node support.
         */
        NodeSupportIterator begin() const { return NodeSupportIterator(T_, rootNodeId_, expectedTopologyVersion_, expectedProperPartVersion_); }

        /**
         * @brief Returns the node-support range sentinel.
         *
         * @return The node-support range sentinel.
         */
        NodeSupportIterator end() const { return NodeSupportIterator(); }
    };

    /**
     * @brief Post-order iterator over one subtree.
     */
    class PostOrderNodeIterator {
      private:
        /** @brief Stores one traversal stack entry and its expansion state. */
        struct Item {
            /** @brief Dense node identifier of the identifier. */
            NodeId id;
            /** @brief Indicates whether the node's children were already expanded. */
            bool expanded;
        };
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Stack buffer. */
        std::vector<Item> stack_;
        /** @brief Dense node identifier of the current. */
        NodeId current_ = InvalidNode;
        /** @brief Expected version used to detect stale derived state. */
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
                    for (NodeId childId : T_->children(top.id)) {
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
         * @param tree Tree topology.
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
         * @return Mutable reference to the updated object.
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
        /** @brief Dense node identifier of the root node identifier. */
        NodeId rootNodeId_ = InvalidNode;
        /** @brief Expected version used to detect stale derived state. */
        std::size_t expectedVersion_ = 0;

      public:
        /**
         * @brief Creates an empty post-order range.
         */
        PostOrderNodeRange() = default;

        /**
         * @brief Creates a post-order range rooted at `rootNodeId`.
         *
         * @param tree Tree topology.
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
        /** @brief Dense node identifier of the queue. */
        FastQueue<NodeId> queue_;
        /** @brief Expected version used to detect stale derived state. */
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
         * @param tree Tree topology.
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
         * @return Mutable reference to the updated object.
         */
        BreadthFirstNodeIterator& operator++() {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            if (!queue_.empty()) {
                NodeId current = queue_.pop();
                for (NodeId childId : T_->children(current)) {
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
        /** @brief Dense node identifier of the root node identifier. */
        NodeId rootNodeId_ = InvalidNode;
        /** @brief Expected version used to detect stale derived state. */
        std::size_t expectedVersion_ = 0;

      public:
        /**
         * @brief Creates an empty breadth-first range.
         */
        BreadthFirstNodeRange() = default;

        /**
         * @brief Creates a breadth-first range rooted at `rootNodeId`.
         *
         * @param tree Tree topology.
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
    class AncestorNodeIterator {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Dense node identifier of the current. */
        NodeId current_ = InvalidNode;
        /** @brief Expected version used to detect stale derived state. */
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
        AncestorNodeIterator() = default;

        /**
         * @brief Creates an iterator starting at one node and walking to the root.
         *
         * @param tree Tree topology.
         * @param current Current iterator or traversal position.
         * @param expectedVersion Mutation version required by the operation.
         */
        AncestorNodeIterator(const MorphologicalTree* tree, NodeId current, std::size_t expectedVersion)
            : T_(tree), current_(current), expectedVersion_(expectedVersion) {}

        /**
         * @brief Advances one step toward the root.
         *
         * @return Mutable reference to the updated object.
         */
        AncestorNodeIterator& operator++() {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            if (T_ && current_ != InvalidNode) {
                if (T_->isRoot(current_)) {
                    current_ = InvalidNode;
                } else {
                    current_ = T_->parent(current_);
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
        bool operator==(const AncestorNodeIterator& other) const { return current_ == other.current_; }

        /**
         * @brief Compares path-to-root iterator positions for inequality.
         *
         * @param other Object to compare with or transfer from.
         * @return True when the documented condition holds; otherwise false.
         */
        bool operator!=(const AncestorNodeIterator& other) const { return !(*this == other); }
    };

    /**
     * @brief Range wrapper for rootward path traversal.
     */
    class AncestorNodeRange {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Dense node identifier of the start. */
        NodeId start_ = InvalidNode;
        /** @brief Expected version used to detect stale derived state. */
        std::size_t expectedVersion_ = 0;

      public:
        /**
         * @brief Creates an empty rootward-path range.
         */
        AncestorNodeRange() = default;

        /**
         * @brief Creates a range from `start` to the connected root.
         *
         * @param tree Tree topology.
         * @param start Inclusive start position.
         * @param expectedVersion Mutation version required by the operation.
         */
        AncestorNodeRange(const MorphologicalTree* tree, NodeId start, std::size_t expectedVersion)
            : T_(tree), start_(start), expectedVersion_(expectedVersion) {}

        /**
         * @brief Returns an iterator at the path start node.
         *
         * @return An iterator at the path start node.
         */
        AncestorNodeIterator begin() const { return AncestorNodeIterator(T_, start_, expectedVersion_); }

        /**
         * @brief Returns the rootward-path range sentinel.
         *
         * @return The rootward-path range sentinel.
         */
        AncestorNodeIterator end() const { return AncestorNodeIterator(); }
    };

    /**
     * @brief Iterator over a materialised path between two nodes.
     */
    class PathBetweenNodesIterator {
      private:
        /** @brief References the t used by the component. */
        const MorphologicalTree* T_ = nullptr;
        /** @brief Dense node identifier of the path. */
        const std::vector<NodeId>* path_ = nullptr;
        /** @brief Index. */
        std::size_t index_ = 0;
        /** @brief Expected version used to detect stale derived state. */
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
         * @param tree Tree topology.
         * @param path Node path updated by the edit.
         * @param index Zero-based index.
         * @param expectedVersion Mutation version required by the operation.
         */
        PathBetweenNodesIterator(const MorphologicalTree* tree, const std::vector<NodeId>* path, std::size_t index, std::size_t expectedVersion)
            : T_(tree), path_(path), index_(index), expectedVersion_(expectedVersion) {}

        /**
         * @brief Advances to the next node in the materialised path.
         *
         * @return Mutable reference to the updated object.
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
        /** @brief Dense node identifier of the path. */
        std::vector<NodeId> path_;
        /** @brief Expected version used to detect stale derived state. */
        std::size_t expectedVersion_ = 0;

        /**
         * @brief Materialises the path joining two nodes in the same connected component.
         *
         * @param tree Tree topology.
         * @param sourceNodeId Node identifier.
         * @param targetNodeId Node identifier.
         * @return Values produced by the operation.
         */
        static std::vector<NodeId> buildPath(const MorphologicalTree* tree, NodeId sourceNodeId, NodeId targetNodeId) {
            if (tree == nullptr || !tree->isAlive(sourceNodeId) || !tree->isAlive(targetNodeId)) {
                return {};
            }

            auto componentAnchor = [tree](NodeId nodeId) {
                NodeId currentNodeId = nodeId;
                while (true) {
                    const NodeId parentNodeId = tree->parent(currentNodeId);
                    if (parentNodeId == InvalidNode || parentNodeId == currentNodeId) {
                        return currentNodeId;
                    }
                    currentNodeId = parentNodeId;
                }
            };

            if (componentAnchor(sourceNodeId) != componentAnchor(targetNodeId)) {
                return {};
            }

            const NodeId lcaNodeId = tree->lowestCommonAncestor(sourceNodeId, targetNodeId);
            if (lcaNodeId == InvalidNode) {
                return {};
            }

            std::vector<NodeId> path;
            for (NodeId currentNodeId = sourceNodeId;; currentNodeId = tree->parent(currentNodeId)) {
                path.push_back(currentNodeId);
                if (currentNodeId == lcaNodeId) {
                    break;
                }

                const NodeId parentNodeId = tree->parent(currentNodeId);
                if (parentNodeId == InvalidNode || parentNodeId == currentNodeId) {
                    return {};
                }
            }

            std::vector<NodeId> descendingTail;
            for (NodeId currentNodeId = targetNodeId; currentNodeId != lcaNodeId; currentNodeId = tree->parent(currentNodeId)) {
                descendingTail.push_back(currentNodeId);

                const NodeId parentNodeId = tree->parent(currentNodeId);
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
         * @param tree Tree topology.
         * @param sourceNodeId Node identifier.
         * @param targetNodeId Node identifier.
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
        /** @brief Dense node identifier of the stack. */
        std::vector<NodeId> stack_;
        /** @brief Expected version used to detect stale derived state. */
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
         * @param tree Tree topology.
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
         * @return Mutable reference to the updated object.
         */
        SubtreeNodeIterator& operator++() {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            if (!stack_.empty()) {
                NodeId current = stack_.back();
                stack_.pop_back();
                std::vector<NodeId> children;
                for (NodeId childId : T_->children(current)) {
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
        /** @brief Dense node identifier of the root node identifier. */
        NodeId rootNodeId_ = InvalidNode;
        /** @brief Expected version used to detect stale derived state. */
        std::size_t expectedVersion_ = 0;

      public:
        /**
         * @brief Creates an empty subtree range.
         */
        SubtreeNodeRange() = default;

        /**
         * @brief Creates a pre-order range over the subtree rooted at `rootNodeId`.
         *
         * @param tree Tree topology.
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
        /** @brief Dense node identifier of the root node identifier. */
        NodeId rootNodeId_ = InvalidNode;
        /** @brief Expected version used to detect stale derived state. */
        std::size_t expectedVersion_ = 0;

      public:
        /**
         * @brief Creates an empty descendant range.
         */
        DescendantNodeRange() = default;

        /**
         * @brief Creates a range over proper descendants of `rootNodeId`.
         *
         * @param tree Tree topology.
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
