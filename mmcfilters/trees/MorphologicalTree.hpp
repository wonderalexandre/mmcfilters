#pragma once

#include "../utils/AdjacencyRelation.hpp"
#include "../utils/Common.hpp"
#include "BuilderMorphologicalTreeByUnionFind.hpp"
#include <cassert>
#include <memory>
#include <stdexcept>
#include <utility>

namespace mmcfilters {

/**
 * @brief Selects the node-id domain used by attribute buffers exposed to callers.
 */
enum class NodeIdSpace {
    MORPHOLOGICAL_TREE,
    HIGRA
};

/**
 * @brief Selects the interpolation mode used to build a tree of shapes.
 */
enum class ToSInterpolation {
    SelfDual,
    Min4cMax8c
};

// Forward declarations for helper wrappers.
class MorphologicalTree;
class MorphologicalTreePybind;
class TreeEditor;
class WeightedMorphologicalTree;

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
 *   `WeightedMorphologicalTree`.
 *
 * The class intentionally keeps only the canonical structural state and a small
 * number of derived caches. Higher-level attribute computation is delegated to
 * the incremental attribute computers in `mmcfilters/attributes`.
 */
class MorphologicalTree {
private:
    // ========================= Private attributes ========================= //
    friend class MorphologicalTreePybind;
    friend class MorphologicalTreeContext;
    friend class TreeEditor;
    friend class WeightedMorphologicalTree;

    struct PrePostOrderCache {
        std::vector<int> timePreOrder;
        std::vector<int> timePostOrder;
        bool valid = false;

        /**
         * @brief Marks the cached traversal timestamps as stale.
         */
        void invalidate() noexcept {
            valid = false;
        }
    };

    class LCAEulerRMQ;

    NodeId rootNodeId_ = InvalidNode;
    int treeType_ = 0; //0-mintree, 1-maxtree, 2-tree of shapes
    int numRows_ = 0;
    int numCols_ = 0;
    std::optional<AdjacencyRelation> adj_; // optional adjacency context; not part of the structural tree contract
    int numNodes_ = 0;
    bool preservesHigraNodeIdSpace_ = false;

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
    mutable PrePostOrderCache prePostOrderCache_;
    mutable std::unique_ptr<LCAEulerRMQ> lcaCache_;


    // Version counters for iterator invalidation.
    std::size_t nodeStructureVersion_ = 0;
    std::size_t topologyVersion_ = 0;
    std::size_t properPartVersion_ = 0;

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

    inline void invalidateHigraNodeIdSpace() noexcept { preservesHigraNodeIdSpace_ = false; }

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
     * @brief Builds a component-tree or tree-of-shapes topology from an image builder.
     */
    inline void build(const ImageUInt8Ptr& imgPtr, IMorphologicalTreeBuilder& builderUF) {
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
     * @brief Rebuilds the tree topology from an imported static Higra hierarchy.
     *
     * The input follows the Higra convention where leaves occupy
     * `[0, numProperParts)` and internal nodes occupy
     * `[numProperParts, parent.size())`. On success, the imported Higra
     * node-id domain is preserved until the tree topology or proper-part
     * ownership is edited.
     */
    inline void resetFromHigraTopology(std::span<const NodeId> parent, NodeId numProperParts) {
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
     * @brief Counts direct proper parts owned by one node slot.
     */
    inline int countOwnedProperPartsBySlot(NodeId slotNodeId) const noexcept {
        int count = 0;
        for (NodeId ownerSlotId : properPartOwner_) {
            if (ownerSlotId == slotNodeId) {
                ++count;
            }
        }
        return count;
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
        invalidateLcaCache();
        invalidateHigraNodeIdSpace();
    }

    /**
     * @brief Marks topology-derived traversals and LCA queries as stale.
     */
    inline void bumpTopologyVersion() noexcept {
        ++topologyVersion_;
        invalidateLcaCache();
        invalidatePrePostOrderCache();
        invalidateHigraNodeIdSpace();
    }

    /**
     * @brief Marks direct proper-part iterators as stale.
     */
    inline void bumpProperPartVersion() noexcept {
        ++properPartVersion_;
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
     * @brief Returns whether linking `candidateParentId` under `forbiddenAncestorId`
     * would close a parent-chain cycle.
     *
     * This helper intentionally walks the current parent chain directly instead
     * of using ancestry caches so it stays cheap and valid while topology
     * caches are stale during local mutations.
     */
    inline bool wouldCreateParentCycle(NodeId forbiddenAncestorId, NodeId candidateParentId) const noexcept {
        if (forbiddenAncestorId == InvalidNode || candidateParentId == InvalidNode) {
            return false;
        }

        NodeId currentNodeId = candidateParentId;
        while (currentNodeId != InvalidNode) {
            if (currentNodeId == forbiddenAncestorId) {
                return true;
            }

            if (!isNode(currentNodeId)) {
                return false;
            }

            const NodeId parentNodeId = nodeParent_[currentNodeId];
            if (parentNodeId == InvalidNode || parentNodeId == currentNodeId) {
                return false;
            }
            currentNodeId = parentNodeId;
        }

        return false;
    }

    /**
     * @brief Debug-only precondition guard against creating a cycle while reparenting.
     *
     * This helper intentionally checks only the candidate parent chain before
     * the mutation. Post-mutation structural validation belongs in the unit
     * tests and debug invariant sweeps, not in the hot path of every mutator.
     */
    inline void assertNoCycleFromReparent([[maybe_unused]] NodeId movingNodeId, [[maybe_unused]] NodeId newParentId) const noexcept {
        assert(!wouldCreateParentCycle(movingNodeId, newParentId) && "MorphologicalTree mutation would create a cycle.");
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
     * @brief Computes subtree areas bottom-up from direct proper-part counts.
     */
    inline std::vector<int32_t> computeAreasIncrementally() const {
        std::vector<int32_t> areaByNode(nodeParent_.size(), 0);

        if (rootNodeId_ == InvalidNode) {
            return areaByNode;
        }

        computeIncrementalAttributes(
            const_cast<MorphologicalTree*>(this),
            getRoot(),
            [&](NodeId nodeId) -> void {
                areaByNode[nodeId] = getNumProperParts(nodeId);
            },
            [&](NodeId parentNodeId, NodeId childNodeId) -> void {
                areaByNode[parentNodeId] += areaByNode[childNodeId];
            },
            [&](NodeId) -> void {}
        );

        return areaByNode;
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
     * @brief Keeps the largest-area descendant associated with one ascendant slot.
     */
    void maxAreaDescendants(const std::vector<int32_t>& areaByNode, std::vector<NodeId>& descendants, NodeId ascendantNodeId, NodeId descendantNodeId) {
        const NodeId nodeAscSlot = ascendantNodeId;
        if (descendants[nodeAscSlot] == InvalidNode) {
            descendants[nodeAscSlot] = descendantNodeId;
        }

        if (areaByNode[descendants[nodeAscSlot]] < areaByNode[descendantNodeId]) {
            descendants[nodeAscSlot] = descendantNodeId;
        }
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

        // Se o filho já tem pai, desconecta antes de ler P/C
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
     * @brief Moves all children of `fromId` to the tail of the child list of `toId`.
     */
    inline void spliceChildrenSlots(NodeId toId, NodeId fromId) {
        if (toId < 0 || fromId < 0 || toId == fromId) return;

        NodeId firstFrom = firstChild_[fromId];
        if (firstFrom == InvalidNode) return; // nada para mover

        // 1) todos os filhos de 'fromId' passam a ter pai 'toId'
        for (NodeId c = firstChild_[fromId]; c != InvalidNode; c = nextSibling_[c]) {
            nodeParent_[c] = toId;
        }

        // 2) concatena a lista de filhos de 'fromId' no fim da lista de 'toId'
        if (firstChild_[toId] == InvalidNode) {
            // 'toId' não tinha filhos — vira exatamente a lista de 'fromId'
            firstChild_[toId] = firstChild_[fromId];
            lastChild_[toId]  = lastChild_[fromId];
            // o primeiro filho já tem prevSibling_ == InvalidNode porque era o primeiro de 'fromId'
        } else {
            // 'toId' já tinha filhos — encadeia no final
            nextSibling_[ lastChild_[toId] ] = firstChild_[fromId];
            prevSibling_[ firstChild_[fromId] ] = lastChild_[toId];
            lastChild_[toId] = lastChild_[fromId];
        }

        // 3) atualiza contadores
        numChildrenByNode_[toId] += numChildrenByNode_[fromId];

        // 4) zera a lista de 'fromId'
        firstChild_[fromId] = InvalidNode;
        lastChild_[fromId]  = InvalidNode;
        numChildrenByNode_[fromId]   = 0;
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
        if (numChildrenByNode_[nodeSlot] != 0 || countOwnedProperPartsBySlot(nodeSlot) != 0) {
            return;
        }
        releaseSlotNode(nodeSlot);
    }

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
     *
     * In debug builds we assert only the reparenting precondition that prevents
     * an immediate cycle. Full post-mutation invariant checking is covered by
     * the dedicated regression tests.
     */
    inline void attachNode(NodeId parentNodeId, NodeId nodeId) {
        if (!isAlive(parentNodeId) || !isAlive(nodeId) || isRoot(nodeId) || parentNodeId == nodeId) {
            return;
        }
        assertNoCycleFromReparent(nodeId, parentNodeId);
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
     *
     * The debug assertion is a precondition check only; deeper invariant
     * validation is intentionally delegated to the test suite.
     */
    inline void moveNode(NodeId nodeId, NodeId newParentId) {
        if (!isAlive(nodeId) || !isAlive(newParentId) || isRoot(nodeId) || nodeId == newParentId) {
            return;
        }
        assertNoCycleFromReparent(nodeId, newParentId);
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
     *
     * As with the other low-level mutators, the runtime debug check is limited
     * to the precondition needed to avoid closing a cycle immediately.
     */
    inline void moveChildren(NodeId parentNodeId, NodeId sourceId) {
        if (!isAlive(parentNodeId) || !isAlive(sourceId) || parentNodeId == sourceId) {
            return;
        }
        assertNoCycleFromReparent(sourceId, parentNodeId);
        spliceChildrenSlots(parentNodeId, sourceId);
    }

    /**
     * @brief Transfers one direct proper part from `sourceNodeId` to `targetNodeId`.
     */
    inline void moveProperPart(NodeId targetNodeId, NodeId sourceNodeId, NodeId pixelId) {
        if (!isAlive(targetNodeId) || !isAlive(sourceNodeId) || !isProperPart(pixelId)) {
            return;
        }
        const NodeId sourceSlotId = sourceNodeId;
        if (properPartOwner_[pixelId] != sourceSlotId) {
            return;
        }
        properPartOwner_[pixelId] = targetNodeId;
        rebuildProperPartLinksFromOwnership();
        bumpProperPartVersion();
    }

    /**
     * @brief Transfers all direct proper parts from `sourceNodeId` to `targetNodeId`.
     */
    inline void moveProperParts(NodeId targetNodeId, NodeId sourceNodeId) {
        if (!isAlive(targetNodeId) || !isAlive(sourceNodeId) || targetNodeId == sourceNodeId) {
            return;
        }
        std::vector<NodeId> sourcePixels;
        for (NodeId pixelId : getProperParts(sourceNodeId)) {
            sourcePixels.push_back(pixelId);
        }
        if (sourcePixels.empty()) {
            return;
        }
        const NodeId targetSlotId = targetNodeId;
        for (NodeId pixelId : sourcePixels) {
            properPartOwner_[pixelId] = targetSlotId;
        }
        rebuildProperPartLinksFromOwnership();
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


public:
    // ========================= Public attributes ========================= //

    static const int MAX_TREE = 0;
    static const int MIN_TREE = 1;
    static const int TREE_OF_SHAPES = 2;

    // ========================= Public constructors ========================= //
    /**
     * @brief Constructs an empty tree shell.
     */
    MorphologicalTree() = default;

    MorphologicalTree(const MorphologicalTree&) = delete;
    MorphologicalTree& operator=(const MorphologicalTree&) = delete;
    MorphologicalTree(MorphologicalTree&&) noexcept = default;
    MorphologicalTree& operator=(MorphologicalTree&&) noexcept = default;

    virtual ~MorphologicalTree() = default;

    // Forward declarations for nested iterator/range types whose definitions stay at the end of the class.
    class AliveNodeIterator;
    class AliveNodeRange;
    class ChildrenIterator;
    class ChildrenRange;
    class ProperPartsIterator;
    class ProperPartsRange;
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

    // ========================= Public methods ========================= //

    /**
     * @brief Creates a max-tree or min-tree from an image.
     */
    static MorphologicalTree createComponentTree(ImageUInt8Ptr img, bool isMaxtree, double radius = 1.5) {
        MorphologicalTree tree;
        tree.treeType_ = isMaxtree ? MAX_TREE : MIN_TREE;
        tree.numRows_ = img->getNumRows();
        tree.numCols_ = img->getNumCols();
        tree.adj_.emplace(tree.numRows_, tree.numCols_, radius);
        tree.numNodes_ = 0;
        tree.properPartOwner_.resize(img->getSize(), InvalidNode);

        BuilderComponentTree builderUF(&*tree.adj_, isMaxtree);
        tree.build(img, builderUF);
        return tree;
    }

    /**
     * @brief Creates a tree of shapes from an image.
     */
    static MorphologicalTree createTreeOfShapes(ImageUInt8Ptr img, ToSInterpolation interpolation = ToSInterpolation::SelfDual) {
        MorphologicalTree tree;
        tree.treeType_ = TREE_OF_SHAPES;
        tree.numRows_ = img->getNumRows();
        tree.numCols_ = img->getNumCols();
        tree.adj_ = std::nullopt;
        tree.numNodes_ = 0;
        tree.properPartOwner_.resize(img->getSize(), InvalidNode);

        BuilderTreeOfShape builderUF(interpolation == ToSInterpolation::Min4cMax8c);
        tree.build(img, builderUF);
        return tree;
    }

    /**
     * @brief Creates a tree from an imported static Higra parent representation.
     *
     * Leaves/proper parts occupy `[0, rows * cols)` and internal nodes occupy
     * `[rows * cols, parent.size())`.
     *
     * The imported Higra node-id space is preserved only while the topology
     * and proper-part ownership remain unmodified. After any edit,
     * `exportHigraHierarchy()` should be used to create a new compact Higra
     * representation.
     */
    static MorphologicalTree createFromHigraParent(
        std::span<const NodeId> parent,
        int rows,
        int cols,
        int treeType,
        std::optional<AdjacencyRelation> adjacency = std::nullopt) {
        if (treeType != MAX_TREE && treeType != MIN_TREE && treeType != TREE_OF_SHAPES) {
            throw std::invalid_argument("Unsupported morphological tree type.");
        }

        MorphologicalTree tree;
        tree.treeType_ = treeType;
        tree.numRows_ = rows;
        tree.numCols_ = cols;
        if (adjacency) {
            tree.adj_ = std::move(adjacency);
        } else if (treeType == TREE_OF_SHAPES) {
            tree.adj_ = std::nullopt;
        } else {
            throw std::invalid_argument("Higra import of max/min trees requires an explicit adjacency relation.");
        }
        tree.resetFromHigraTopology(parent, static_cast<NodeId>(rows * cols));
        return tree;
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
        return numChildrenByNode_[nodeId];
    }

    /**
     * @brief Returns the number of internal descendants of `nodeId`.
     */
    inline int getNodeNumDescendants(NodeId nodeId) const {
        ensurePrePostOrderCache();
        const NodeId localId = nodeId;
        return (prePostOrderCache_.timePostOrder[localId] - prePostOrderCache_.timePreOrder[localId] - 1) / 2;
    }

    /**
     * @brief Returns the number of siblings of `nodeId`.
     */
    inline int getNodeNumSiblings(NodeId nodeId) const noexcept {
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
        ensurePrePostOrderCache();
        return prePostOrderCache_.timePreOrder[nodeId];
    }

    /**
     * @brief Returns the postorder time of `nodeId` in the cached DFS traversal.
     */
    inline int getNodeTimePostOrder(NodeId nodeId) const {
        ensurePrePostOrderCache();
        return prePostOrderCache_.timePostOrder[nodeId];
    }

    /**
     * @brief Returns the first direct child of `nodeId`, or `InvalidNode`.
     */
    inline NodeId getFirstChild(NodeId nodeId) const {
        return firstChild_[nodeId];
    }

    /**
     * @brief Returns the next sibling of `nodeId`, or `InvalidNode`.
     */
    inline NodeId getNextSibling(NodeId nodeId) const {
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
        return numProperPartsByNode_[nodeId];
    }

    /**
     * @brief Tests whether `childId` is a direct child of `parentNodeId`.
     */
    inline bool hasChild(NodeId parentNodeId, NodeId childId) const { return getNodeParent(childId) == parentNodeId; }

    /**
     * @brief Returns the direct parent of `nodeId`.
     *
     * The root and detached nodes report themselves as parent.
     */
    inline NodeId getNodeParent(NodeId nodeId) const {
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
     * @brief Returns the live node that directly owns `pixelId`.
     */
    inline NodeId getSmallestComponent(int pixelId) const {
        return (pixelId >= 0 && pixelId < static_cast<int>(properPartOwner_.size()))
            ? properPartOwner_[pixelId]
            : InvalidNode;
    }

    /**
     * @brief Returns all live leaf nodes in the current hierarchy.
     */
    std::vector<NodeId> getLeaves() {
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
     * @brief Returns the current tree type (`MAX_TREE`, `MIN_TREE`, or `TREE_OF_SHAPES`).
     */
    inline int getTreeType() const noexcept{return treeType_;}

    /**
     * @brief Tests whether the instance represents a max-tree.
     */
    inline bool isMaxtree()const noexcept{ return treeType_ == MAX_TREE;}

    /**
     * @brief Returns the number of currently live nodes.
     */
    inline int getNumNodes()const noexcept{ return numNodes_; }

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
     * @brief Tests whether `u` is an ancestor of `v`.
     */
    inline bool isAncestor(NodeId u, NodeId v) const {
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
        return ChildrenRange(this, firstChild_[nodeId], topologyVersion_);
    }

    /**
     * @brief Returns a fail-fast range over the direct proper parts of `nodeId`.
     */
    inline ProperPartsRange getProperParts(NodeId nodeId) const {
        return ProperPartsRange(this, properHead_[nodeId], properPartVersion_);
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
        return BreadthFirstNodeRange(this, rootNodeId, topologyVersion_);
    }

    /**
     * @brief Returns the path from `nodeId` to the connected root.
     */
    inline PathToRootRange getPathToRootNodes(NodeId nodeId) const {
        return PathToRootRange(this, nodeId, topologyVersion_);
    }

    /**
     * @brief Returns the path that connects `sourceNodeId` and `targetNodeId`.
     */
    inline PathBetweenNodesRange getPathBetweenNodes(NodeId sourceNodeId, NodeId targetNodeId) const {
        return PathBetweenNodesRange(this, sourceNodeId, targetNodeId, topologyVersion_);
    }

    /**
     * @brief Returns a pre-order traversal range over the subtree of `nodeId`.
     */
    inline SubtreeNodeRange getNodeSubtree(NodeId nodeId) const {
        return SubtreeNodeRange(this, nodeId, topologyVersion_);
    }

    /**
     * @brief Returns a range over all proper descendants of `nodeId`.
     */
    inline DescendantNodeRange getDescendants(NodeId nodeId) const {
        return DescendantNodeRange(this, nodeId, topologyVersion_);
    }

    /**
     * @brief Prunes the subtree of `nodeId`, moving all its support to the parent.
     */
    inline void pruneNode(NodeId nodeId) {
        if (!isAlive(nodeId) || isRoot(nodeId)) {
            return;
        }
        const NodeId parentNodeId = getNodeParent(nodeId);
        const NodeId parentSlotId = parentNodeId;

        std::vector<NodeId> subtreeNodes;
        std::vector<NodeId> movedPixels;
        for (NodeId subtreeNodeId : getPostOrderNodes(nodeId)) {
            subtreeNodes.push_back(subtreeNodeId);
            for (NodeId pixelId : getProperParts(subtreeNodeId)) {
                movedPixels.push_back(pixelId);
            }
        }
        for (NodeId pixelId : movedPixels) {
            properPartOwner_[pixelId] = parentSlotId;
        }
        rebuildProperPartLinksFromOwnership();
        bumpProperPartVersion();

        for (NodeId subtreeNodeId : subtreeNodes) {
            const NodeId subtreeSlotId = subtreeNodeId;
            const NodeId currentParentSlotId = nodeParent_[subtreeSlotId];
            if (currentParentSlotId != InvalidNode && currentParentSlotId != subtreeSlotId) {
                unlinkChildSlot(currentParentSlotId, subtreeSlotId, false);
            }
            nodeParent_[subtreeSlotId] = subtreeSlotId;
            releaseNode(subtreeNodeId);
        }
    }

    /**
     * @brief Merges `nodeId` into its parent and releases the emptied slot.
     */
    inline void mergeNodeIntoParent(NodeId nodeId) {
        if (!isAlive(nodeId) || isRoot(nodeId)) {
            return;
        }
        const NodeId parentNodeId = getNodeParent(nodeId);
        const NodeId parentSlotId = parentNodeId;
        const NodeId nodeSlotId = nodeId;

        std::vector<NodeId> movedPixels;
        for (NodeId pixelId : getProperParts(nodeId)) {
            movedPixels.push_back(pixelId);
        }
        for (NodeId pixelId : movedPixels) {
            properPartOwner_[pixelId] = parentSlotId;
        }
        unlinkChildSlot(parentSlotId, nodeSlotId, false);
        spliceChildrenSlots(parentSlotId, nodeSlotId);
        nodeParent_[nodeSlotId] = nodeSlotId;
        rebuildProperPartLinksFromOwnership();
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
        std::vector<std::vector<int>> sparseTable_;
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
                sparseTable_.clear();
                return;
            }

            int maxLog = 1;
            while ((1 << maxLog) <= n) {
                ++maxLog;
            }

            sparseTable_.assign(static_cast<size_t>(n), std::vector<int>(static_cast<size_t>(maxLog), 0));

            for (int i = 0; i < n; ++i) {
                sparseTable_[static_cast<size_t>(i)][0] = i;
            }

            for (int j = 1; j < maxLog; ++j) {
                const int blockSize = 1 << j;
                const int halfBlock = blockSize >> 1;
                for (int i = 0; i + blockSize <= n; ++i) {
                    const int leftIndex = sparseTable_[static_cast<size_t>(i)][static_cast<size_t>(j - 1)];
                    const int rightIndex = sparseTable_[static_cast<size_t>(i + halfBlock)][static_cast<size_t>(j - 1)];
                    sparseTable_[static_cast<size_t>(i)][static_cast<size_t>(j)] =
                        depth_[static_cast<size_t>(leftIndex)] <= depth_[static_cast<size_t>(rightIndex)] ? leftIndex : rightIndex;
                }
            }
        }

        /**
         * @brief Returns the Euler-tour position of minimum depth on `[left, right]`.
         */
        int rmq(int left, int right) const {
            const int length = right - left + 1;
            int logLength = 0;
            while ((1 << (logLength + 1)) <= length) {
                ++logLength;
            }

            const int leftIndex = sparseTable_[static_cast<size_t>(left)][static_cast<size_t>(logLength)];
            const int rightIndex = sparseTable_[static_cast<size_t>(right - (1 << logLength) + 1)][static_cast<size_t>(logLength)];
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
        using iterator_category = std::input_iterator_tag;
        using value_type = NodeId;
        using difference_type = std::ptrdiff_t;
        using pointer = const NodeId*;
        using reference = const NodeId&;

        AliveNodeIterator() = default;
        AliveNodeIterator(const MorphologicalTree* tree, NodeId current, NodeId end, std::size_t expectedVersion)
            : T_(tree), current_(current), end_(end), expectedVersion_(expectedVersion) {
            settle_();
        }

        AliveNodeIterator& operator++() {
            T_->checkNodeIteratorVersion(expectedVersion_);
            if (current_ != InvalidNode) {
                ++current_;
                settle_();
            }
            return *this;
        }

        NodeId operator*() const {
            T_->checkNodeIteratorVersion(expectedVersion_);
            return current_;
        }
        bool operator==(const AliveNodeIterator& other) const { return current_ == other.current_; }
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
        AliveNodeRange() = default;
        AliveNodeRange(const MorphologicalTree* tree, NodeId begin, NodeId end, std::size_t expectedVersion)
            : T_(tree), begin_(begin), end_(end), expectedVersion_(expectedVersion) {}

        AliveNodeIterator begin() const { return AliveNodeIterator(T_, begin_, end_, expectedVersion_); }
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
        using iterator_category = std::forward_iterator_tag;
        using value_type = NodeId;
        using difference_type = std::ptrdiff_t;
        using pointer = const NodeId*;
        using reference = const NodeId&;

        ChildrenIterator() = default;
        ChildrenIterator(const MorphologicalTree* tree, NodeId currentLocal, std::size_t expectedVersion)
            : T_(tree), currentLocal_(currentLocal), expectedVersion_(expectedVersion) {}

        ChildrenIterator& operator++() {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            if (T_ && currentLocal_ != InvalidNode) {
                currentLocal_ = T_->nextSibling_[currentLocal_];
            }
            return *this;
        }

        NodeId operator*() const {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            return currentLocal_;
        }
        bool operator==(const ChildrenIterator& other) const { return currentLocal_ == other.currentLocal_; }
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
        ChildrenRange() = default;
        ChildrenRange(const MorphologicalTree* tree, NodeId firstLocal, std::size_t expectedVersion)
            : T_(tree), firstLocal_(firstLocal), expectedVersion_(expectedVersion) {}

        ChildrenIterator begin() const { return ChildrenIterator(T_, firstLocal_, expectedVersion_); }
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
        using iterator_category = std::forward_iterator_tag;
        using value_type = NodeId;
        using difference_type = std::ptrdiff_t;
        using pointer = const NodeId*;
        using reference = const NodeId&;

        ProperPartsIterator() = default;
        ProperPartsIterator(const MorphologicalTree* tree, NodeId currentPixel, std::size_t expectedVersion)
            : T_(tree), currentPixel_(currentPixel), expectedVersion_(expectedVersion) {}

        ProperPartsIterator& operator++() {
            T_->checkProperPartIteratorVersion(expectedVersion_);
            if (T_ && currentPixel_ != InvalidNode) {
                currentPixel_ = T_->nextProperPart_[currentPixel_];
            }
            return *this;
        }

        NodeId operator*() const {
            T_->checkProperPartIteratorVersion(expectedVersion_);
            return currentPixel_;
        }
        bool operator==(const ProperPartsIterator& other) const { return currentPixel_ == other.currentPixel_; }
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
        ProperPartsRange() = default;
        ProperPartsRange(const MorphologicalTree* tree, NodeId firstPixel, std::size_t expectedVersion)
            : T_(tree), firstPixel_(firstPixel), expectedVersion_(expectedVersion) {}

        ProperPartsIterator begin() const { return ProperPartsIterator(T_, firstPixel_, expectedVersion_); }
        ProperPartsIterator end() const { return ProperPartsIterator(T_, InvalidNode, expectedVersion_); }
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
        using iterator_category = std::input_iterator_tag;
        using value_type = NodeId;
        using difference_type = std::ptrdiff_t;
        using pointer = const NodeId*;
        using reference = const NodeId&;

        PostOrderNodeIterator() = default;
        PostOrderNodeIterator(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedVersion) : T_(tree), expectedVersion_(expectedVersion) {
            if (T_ && rootNodeId != InvalidNode) {
                stack_.push_back({rootNodeId, false});
                settle_();
            }
        }

        PostOrderNodeIterator& operator++() {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            if (!stack_.empty()) {
                stack_.pop_back();
                settle_();
            }
            return *this;
        }

        NodeId operator*() const {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            return current_;
        }
        bool operator==(const PostOrderNodeIterator& other) const { return current_ == other.current_; }
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
        PostOrderNodeRange() = default;
        PostOrderNodeRange(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedVersion)
            : T_(tree), rootNodeId_(rootNodeId), expectedVersion_(expectedVersion) {}

        PostOrderNodeIterator begin() const { return PostOrderNodeIterator(T_, rootNodeId_, expectedVersion_); }
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
        using iterator_category = std::input_iterator_tag;
        using value_type = NodeId;
        using difference_type = std::ptrdiff_t;
        using pointer = const NodeId*;
        using reference = const NodeId&;

        BreadthFirstNodeIterator() = default;
        BreadthFirstNodeIterator(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedVersion) : T_(tree), expectedVersion_(expectedVersion) {
            if (T_ && rootNodeId != InvalidNode) {
                queue_.push(rootNodeId);
            }
        }

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

        NodeId operator*() const {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            return queue_.front();
        }
        bool operator==(const BreadthFirstNodeIterator& other) const { return queue_.empty() == other.queue_.empty(); }
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
        BreadthFirstNodeRange() = default;
        BreadthFirstNodeRange(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedVersion)
            : T_(tree), rootNodeId_(rootNodeId), expectedVersion_(expectedVersion) {}

        BreadthFirstNodeIterator begin() const { return BreadthFirstNodeIterator(T_, rootNodeId_, expectedVersion_); }
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
        using iterator_category = std::input_iterator_tag;
        using value_type = NodeId;
        using difference_type = std::ptrdiff_t;
        using pointer = const NodeId*;
        using reference = const NodeId&;

        PathToRootIterator() = default;
        PathToRootIterator(const MorphologicalTree* tree, NodeId current, std::size_t expectedVersion)
            : T_(tree), current_(current), expectedVersion_(expectedVersion) {}

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

        NodeId operator*() const {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            return current_;
        }
        bool operator==(const PathToRootIterator& other) const { return current_ == other.current_; }
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
        PathToRootRange() = default;
        PathToRootRange(const MorphologicalTree* tree, NodeId start, std::size_t expectedVersion)
            : T_(tree), start_(start), expectedVersion_(expectedVersion) {}

        PathToRootIterator begin() const { return PathToRootIterator(T_, start_, expectedVersion_); }
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
        using iterator_category = std::input_iterator_tag;
        using value_type = NodeId;
        using difference_type = std::ptrdiff_t;
        using pointer = const NodeId*;
        using reference = const NodeId&;

        PathBetweenNodesIterator() = default;
        PathBetweenNodesIterator(
            const MorphologicalTree* tree,
            const std::vector<NodeId>* path,
            std::size_t index,
            std::size_t expectedVersion)
            : T_(tree), path_(path), index_(index), expectedVersion_(expectedVersion) {}

        PathBetweenNodesIterator& operator++() {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            if (path_ && index_ < path_->size()) {
                ++index_;
            }
            return *this;
        }

        NodeId operator*() const {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            return (*path_)[index_];
        }

        bool operator==(const PathBetweenNodesIterator& other) const {
            return path_ == other.path_ && index_ == other.index_;
        }
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
        PathBetweenNodesRange() = default;
        PathBetweenNodesRange(
            const MorphologicalTree* tree,
            NodeId sourceNodeId,
            NodeId targetNodeId,
            std::size_t expectedVersion)
            : T_(tree),
              path_(buildPath(tree, sourceNodeId, targetNodeId)),
              expectedVersion_(expectedVersion) {}

        PathBetweenNodesIterator begin() const {
            return PathBetweenNodesIterator(T_, &path_, 0, expectedVersion_);
        }
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
        using iterator_category = std::input_iterator_tag;
        using value_type = NodeId;
        using difference_type = std::ptrdiff_t;
        using pointer = const NodeId*;
        using reference = const NodeId&;

        SubtreeNodeIterator() = default;
        SubtreeNodeIterator(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedVersion) : T_(tree), expectedVersion_(expectedVersion) {
            if (T_ && rootNodeId != InvalidNode) {
                stack_.push_back(rootNodeId);
            }
        }

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

        NodeId operator*() const {
            T_->checkTopologyIteratorVersion(expectedVersion_);
            return stack_.back();
        }
        bool operator==(const SubtreeNodeIterator& other) const { return stack_.empty() == other.stack_.empty(); }
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
        SubtreeNodeRange() = default;
        SubtreeNodeRange(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedVersion)
            : T_(tree), rootNodeId_(rootNodeId), expectedVersion_(expectedVersion) {}

        SubtreeNodeIterator begin() const { return SubtreeNodeIterator(T_, rootNodeId_, expectedVersion_); }
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
        DescendantNodeRange() = default;
        DescendantNodeRange(const MorphologicalTree* tree, NodeId rootNodeId, std::size_t expectedVersion)
            : T_(tree), rootNodeId_(rootNodeId), expectedVersion_(expectedVersion) {}

        SubtreeNodeIterator begin() const {
            auto it = SubtreeNodeIterator(T_, rootNodeId_, expectedVersion_);
            ++it;
            return it;
        }
        SubtreeNodeIterator end() const { return SubtreeNodeIterator(); }
    };



};

} // namespace mmcfilters
