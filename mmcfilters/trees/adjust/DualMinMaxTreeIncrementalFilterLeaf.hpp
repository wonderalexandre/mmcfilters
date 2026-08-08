#pragma once

#include "../WeightedMorphologicalTree.hpp"
#include "../../utils/GenerationStampSet.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

#ifndef MMCFILTERS_COMPONENT_TREE_ADJUSTMENT_DENSE_MAX_BITS
#define MMCFILTERS_COMPONENT_TREE_ADJUSTMENT_DENSE_MAX_BITS 8
#endif

namespace mmcfilters::adjust {

/**
 * @brief Leaf-by-leaf incremental updater for a paired min-tree / max-tree state.
 *
 * This is the leaf-granularity counterpart of
 * `DualMinMaxTreeIncrementalFilter`. It is adapted from the MorphoTreeAdjust
 * `DualMinMaxTreeIncrementalFilterLeaf`, but targets this repository's
 * `WeightedMorphologicalTree<T>` and `WeightedTreeEditor<T>` representation.
 *
 * The final min/max semantics are the same as the subtree updater: one tree is
 * pruned and the dual tree is updated in place so both trees keep representing
 * the same image. The difference is the elementary operation:
 *
 * - here, one live non-root leaf is removed from the primal tree;
 * - the dual tree is updated immediately after that leaf removal;
 * - a larger threshold or subtree operation is expressed as a sequence of these
 *   smaller leaf updates.
 *
 * Consequences of this granularity:
 *
 * - each update touches a smaller local support;
 * - the total number of update events can be much larger;
 * - SDRT construction can record one residual event per consumed leaf support;
 * - the class exposes the dual-tree nodes touched by the last update so callers
 *   can refresh only affected agenda candidates.
 *
 * Algorithmic intuition:
 *
 * - the proper parts of the removed leaf define the set `C`;
 * - the dual tree is queried to find `nodeCa`, the extremal component that
 *   contains `C` at altitude `a`;
 * - graph neighbors of `C` seed ancestor paths in the dual tree;
 * - those paths are bucketed by altitude and swept from the outside of the
 *   interval toward `nodeCa`;
 * - nodes at the same level are merged, frontiers above/below the interval are
 *   reattached, and empty nodes are contracted.
 *
 * Operational invariants expected by the implementation:
 *
 * - `mintree_` and `maxtree_` are externally owned and outlive the adjuster;
 * - both trees share the same image domain and adjacency relation;
 * - the selected leaf is alive, non-root, and has at least one direct proper
 *   part at the time `updateTree` is called;
 * - this is a mutable weighted-tree boundary, not a `WeightedTreeView`
 *   boundary: one update changes topology/proper-part ownership and the
 *   owned weighted state that subsequent steps observe;
 * - the `altitude_t` template parameter is the owned tree altitude type and
 *   also controls merge-bucket keys/backends;
 * - all topology mutations go through one `WeightedTreeEditor<T>` edit session
 *   per dual update and publish with the generic move-only proof protocol; the
 *   update establishes its invariants during the existing algorithmic passes,
 *   and assertion-enabled builds retain the complete oracle.
 *
 * Complexity of one leaf step is local. With `P_L` as the number of proper
 * parts in the leaf, `A_L` as the number of adjacent seed nodes, `M_L` as the
 * number of nodes bucketed for merging, and `K_L` as the number of nodes edited
 * in the dual tree, the practical cost is approximately
 * `O(P_L + A_L + M_L + K_L)`.
 */
template <AltitudeValue altitude_t> class DualMinMaxTreeIncrementalFilterLeaf {
  private:
    /** @brief Defines the `tree_t` alias used by the component. */
    using tree_t = WeightedMorphologicalTree<altitude_t>;
    /** @brief Defines the `editor_t` alias used by the component. */
    using editor_t = WeightedTreeEditor<altitude_t>;

    /**
     * @brief Decides whether altitude buckets can be represented by a dense array.
     *
     * MorphoTreeAdjust assumes an 8-bit non-negative altitude domain and uses one
     * vector per possible level. In MAF, `std::uint8_t` can vary, so this helper
     * keeps that dense backend only for small integral domains and falls back to a
     * sparse map for wider or non-integral altitude types.
     *
     * @return True when the documented condition holds; otherwise false.
     */
    static constexpr bool usesDenseLevels() {
        if constexpr (std::is_integral_v<altitude_t> && !std::is_same_v<std::remove_cv_t<altitude_t>, bool>) {
            using unsigned_altitude_t = std::make_unsigned_t<altitude_t>;
            return std::numeric_limits<unsigned_altitude_t>::digits <= MMCFILTERS_COMPONENT_TREE_ADJUSTMENT_DENSE_MAX_BITS;
        } else {
            return false;
        }
    }

    /** @brief Indicates whether altitude levels use the dense lookup representation. */
    static constexpr bool use_dense_levels = usesDenseLevels();

    /**
     * @brief Per-update collection of dual-tree nodes grouped by altitude.
     *
     * This mirrors the role of MorphoTreeAdjust's `MergedNodesCollection`: during
     * one leaf update, adjacent seed components are climbed toward `nodeCa`, and
     * the nodes that must be merged are placed in altitude buckets. The sweep then
     * consumes those buckets in descending order for max-tree updates and ascending
     * order for min-tree updates.
     *
     * The collection also stores frontier nodes outside the `[a, b]` interval.
     * Those frontier nodes are not merged, but must be reattached to the union node
     * created at level `b` so the dual tree remains connected.
     */
    class MergedNodesCollection {
      private:
        /**
         * @brief Selects dense or sparse bucket storage for the configured altitude type.
         */
        template <bool dense, typename altitude_type> struct StorageSelector;

        /**
         * @brief Dense storage: one bucket for each representable altitude value.
         */
        template <typename altitude_type> struct StorageSelector<true, altitude_type> {
            /** @brief Stores the domain size. */
            static constexpr std::size_t domain_size = static_cast<std::size_t>(static_cast<long long>(std::numeric_limits<altitude_type>::max()) -
                                                                                static_cast<long long>(std::numeric_limits<altitude_type>::lowest()) + 1);
            /** @brief Defines the `type` alias used by the component. */
            using type = std::array<std::vector<NodeId>, domain_size>;
        };

        /**
         * @brief Sparse storage used when the altitude domain is too large to index directly.
         */
        template <typename altitude_type> struct StorageSelector<false, altitude_type> {
            /** @brief Defines the `type` alias used by the component. */
            using type = std::map<altitude_type, std::vector<NodeId>>;
        };

        /** @brief Defines the `storage_t` alias used by the component. */
        using storage_t = typename StorageSelector<use_dense_levels, altitude_t>::type;

        /** @brief Merge buckets keyed by altitude level. */
        storage_t mergeNodesByLevelStorage_;

        /** @brief Levels that currently have non-empty merge buckets in this update. */
        std::vector<altitude_t> mergeLevels_;

        /** @brief Distinct dual-tree components adjacent to the removed leaf support. */
        std::vector<NodeId> adjacentNodes_;

        /** @brief Boundary nodes outside the merge interval that must be reattached at `b`. */
        std::vector<NodeId> frontierNodesAboveB_;

        /** @brief Marks all nodes already collected as either merge nodes or frontiers. */
        GenerationStampSet collectedNodeMarks_;

        /** @brief Marks only nodes inserted in merge buckets. */
        GenerationStampSet mergeBucketNodeMarks_;

        /** @brief Marks dual components already used as adjacency seeds. */
        GenerationStampSet adjacentSeedMarks_;

        /** @brief Largest bucket size seen in the current step, useful for scratch sizing. */
        std::size_t maxBucketSize_ = 0;

        /** @brief Cursor into `mergeLevels_` while the level sweep is running. */
        int currentMergeLevelIndex_ = 0;

        /** @brief Direction flag: max-tree sweeps high-to-low, min-tree sweeps low-to-high. */
        bool isMaxtree_ = false;

      public:
        /**
         * @brief Converts an altitude value into an array index for dense storage.
         *
         * @param level Altitude level used by the operation.
         * @return The converted altitude value into an array index for dense storage.
         */
        static std::size_t denseBucketIndex(altitude_t level) {
            if constexpr (std::is_signed_v<altitude_t>) {
                return static_cast<std::size_t>(static_cast<long long>(level) - static_cast<long long>(std::numeric_limits<altitude_t>::lowest()));
            } else {
                return static_cast<std::size_t>(level);
            }
        }

        /**
         * @brief Converts a dense bucket index back to the corresponding altitude value.
         *
         * @param index Zero-based index used by the operation.
         * @return The converted dense bucket index back to the corresponding altitude value.
         */
        static altitude_t levelFromDenseBucketIndex(std::size_t index) {
            if constexpr (std::is_signed_v<altitude_t>) {
                return static_cast<altitude_t>(static_cast<long long>(index) + static_cast<long long>(std::numeric_limits<altitude_t>::lowest()));
            } else {
                return static_cast<altitude_t>(index);
            }
        }

        /**
         * @brief Creates collection state sized for the fixed node-id domain.
         *
         * Mark sets are generation-stamped, so resetting the collection between
         * leaf updates is normally O(1) for marks and O(number of touched levels)
         * for bucket cleanup.
         *
         * @param maxNodes Max-tree nodes associated with the leaf event.
         */
        explicit MergedNodesCollection(int maxNodes = 0)
            : collectedNodeMarks_(static_cast<std::size_t>(std::max(0, maxNodes))), mergeBucketNodeMarks_(static_cast<std::size_t>(std::max(0, maxNodes))),
              adjacentSeedMarks_(static_cast<std::size_t>(std::max(0, maxNodes))) {}

        /**
         * @brief Clears transient state and sets the sweep direction for the next step.
         *
         * @param isMaxtree Flag controlling is maxtree.
         */
        void resetCollection(bool isMaxtree) {
            isMaxtree_ = isMaxtree;
            for (altitude_t level : mergeLevels_) {
                if constexpr (use_dense_levels) {
                    mergeNodesByLevelStorage_[denseBucketIndex(level)].clear();
                } else {
                    auto it = mergeNodesByLevelStorage_.find(level);
                    if (it != mergeNodesByLevelStorage_.end()) {
                        it->second.clear();
                    }
                }
            }
            mergeLevels_.clear();
            adjacentNodes_.clear();
            frontierNodesAboveB_.clear();
            collectedNodeMarks_.resetAll();
            mergeBucketNodeMarks_.resetAll();
            adjacentSeedMarks_.resetAll();
            maxBucketSize_ = 0;
            currentMergeLevelIndex_ = 0;
        }

        /**
         * @brief Returns the merge bucket for `level`, creating it for sparse storage.
         *
         * @param level Altitude level used by the operation.
         * @return The merge bucket for level, creating it for sparse storage.
         */
        std::vector<NodeId>& getMergedNodes(const altitude_t& level) {
            if constexpr (use_dense_levels) {
                return mergeNodesByLevelStorage_[denseBucketIndex(level)];
            } else {
                return mergeNodesByLevelStorage_[level];
            }
        }

        /**
         * @brief Exposes the distinct adjacent seed nodes collected in this step.
         *
         * @return Reference to the resulting object.
         */
        std::vector<NodeId>& getAdjacentNodes() { return adjacentNodes_; }

        /**
         * @brief Exposes frontier nodes outside the merge interval.
         *
         * @return Reference to the resulting object.
         */
        std::vector<NodeId>& getFrontierNodesAboveB() { return frontierNodesAboveB_; }

        /**
         * @brief Adds `nodeId` as a frontier node if it was not collected before.
         *
         * @param nodeId Identifier of the node used by the operation.
         */
        void addFrontierNodeAboveB(NodeId nodeId) {
            if (nodeId == InvalidNode || collectedNodeMarks_.isMarked(static_cast<std::size_t>(nodeId))) {
                return;
            }
            frontierNodesAboveB_.push_back(nodeId);
            collectedNodeMarks_.mark(static_cast<std::size_t>(nodeId));
        }

        /**
         * @brief Adds `nodeId` to the merge bucket matching its current altitude.
         *
         * @param tree Tree topology used by the operation.
         * @param nodeId Identifier of the node used by the operation.
         */
        void addMergeNode(const tree_t& tree, NodeId nodeId) {
            if (nodeId == InvalidNode || collectedNodeMarks_.isMarked(static_cast<std::size_t>(nodeId))) {
                return;
            }

            auto& bucket = getMergedNodes(tree.getAltitude(nodeId));
            bucket.push_back(nodeId);
            maxBucketSize_ = std::max(maxBucketSize_, bucket.size());
            collectedNodeMarks_.mark(static_cast<std::size_t>(nodeId));
            mergeBucketNodeMarks_.mark(static_cast<std::size_t>(nodeId));
        }

        /**
         * @brief Returns whether `nodeId` was bucketed as a merge node in this step.
         *
         * @param nodeId Identifier of the node used by the operation.
         * @return Whether nodeId was bucketed as a merge node in this step.
         */
        bool isMergeNode(NodeId nodeId) const { return nodeId != InvalidNode && mergeBucketNodeMarks_.isMarked(static_cast<std::size_t>(nodeId)); }

        /**
         * @brief Registers a dual-tree adjacency seed once and records it for inspection.
         *
         * @param nodeId Identifier of the node used by the operation.
         * @return True when the documented condition holds; otherwise false.
         */
        bool markAdjacentSeed(NodeId nodeId) {
            if (nodeId == InvalidNode || adjacentSeedMarks_.isMarked(static_cast<std::size_t>(nodeId))) {
                return false;
            }

            adjacentNodes_.push_back(nodeId);
            adjacentSeedMarks_.mark(static_cast<std::size_t>(nodeId));
            return true;
        }

        /**
         * @brief Returns the largest merge bucket size observed in this update.
         *
         * @return The largest merge bucket size observed in this update.
         */
        std::size_t getMaxBucketSize() const { return maxBucketSize_; }

        /**
         * @brief Materializes non-empty levels and returns the first level to sweep.
         *
         * @return The materialized non-empty levels and returns the first level to sweep.
         */
        altitude_t firstMergeLevel() {
            mergeLevels_.clear();
            if constexpr (use_dense_levels) {
                for (std::size_t i = 0; i < mergeNodesByLevelStorage_.size(); ++i) {
                    if (!mergeNodesByLevelStorage_[i].empty()) {
                        mergeLevels_.push_back(levelFromDenseBucketIndex(i));
                    }
                }
            } else {
                mergeLevels_.reserve(mergeNodesByLevelStorage_.size());
                for (const auto& entry : mergeNodesByLevelStorage_) {
                    if (!entry.second.empty()) {
                        mergeLevels_.push_back(entry.first);
                    }
                }
            }

            if (mergeLevels_.empty()) {
                return altitude_t{};
            }

            currentMergeLevelIndex_ = isMaxtree_ ? static_cast<int>(mergeLevels_.size()) - 1 : 0;
            return mergeLevels_[static_cast<std::size_t>(currentMergeLevelIndex_)];
        }

        /**
         * @brief Tests whether the level sweep cursor still points to a valid level.
         *
         * @return True if the level sweep cursor still points to a valid level; otherwise false.
         */
        bool hasMergeLevel() const {
            return !mergeLevels_.empty() && currentMergeLevelIndex_ >= 0 && currentMergeLevelIndex_ < static_cast<int>(mergeLevels_.size());
        }

        /**
         * @brief Advances the sweep cursor according to the tree polarity.
         *
         * @return Next merge level in the polarity-ordered sweep.
         */
        altitude_t nextMergeLevel() {
            currentMergeLevelIndex_ = isMaxtree_ ? currentMergeLevelIndex_ - 1 : currentMergeLevelIndex_ + 1;
            if (!hasMergeLevel()) {
                return altitude_t{};
            }
            return mergeLevels_[static_cast<std::size_t>(currentMergeLevelIndex_)];
        }
    };

    /** @brief Mutable min-tree owned by the caller. */
    tree_t* mintree_ = nullptr;

    /** @brief Mutable max-tree owned by the caller. */
    tree_t* maxtree_ = nullptr;

    /** @brief Shared pixel adjacency used to find dual neighbors of the removed leaf support. */
    const RegularGridAdjacency2D* graph_ = nullptr;

    /** @brief Merge buckets, adjacency seeds, and frontiers for the current leaf step. */
    MergedNodesCollection mergeNodesByLevel_;

    /** @brief Marks empty dual-tree nodes that should be contracted after the sweep. */
    GenerationStampSet removedMarks_;

    /** @brief Empty marked nodes whose absorption is deferred until local reattachments finish. */
    std::vector<NodeId> removedNodesPendingAbsorption_;

    /** @brief Marks pixels that belong to the primal leaf support `C`. */
    GenerationStampSet pixelsInLeafMarks_;

    /** @brief Marks dual-tree nodes whose ancestor path has already been climbed. */
    GenerationStampSet climbedNodeMarks_;

    /** @brief Dual-tree nodes that became candidate-relevant in the last update. */
    std::vector<NodeId> lastCandidateNodes_;

    /** @brief Deduplicates roots whose current parent relation changed. */
    GenerationStampSet topologyChangedMarks_;

    /**
     * @brief Nodes whose parent edge changed or whose slot died in the last update.
     *
     * Marking a superset is intentional: dynamic ancestry indexes may invalidate
     * the corresponding snapshot subtrees without inspecting the editor journal.
     */
    std::vector<NodeId> lastTopologyChangedNodes_;

    /** @brief Whether the last update transferred the selected support in bulk. */
    bool lastBulkProperPartTransfer_ = false;

    /**
     * @brief Exact external dual owners supplied by a preceding leaf certificate.
     *
     * Min/max residual construction already scans the selected leaf boundary
     * before committing it. These fields let the immediately following dual
     * update reuse those owners instead of traversing the same pixel boundary
     * a second time. Ordinary callers leave the cache empty and retain the
     * original self-contained update path.
     */
    tree_t* primedDualTree_ = nullptr;
    /** @brief Stores the primed leaf identifier. */
    NodeId primedLeafId_ = InvalidNode;
    /** @brief Stores the primed adjacent nodes. */
    std::vector<NodeId> primedAdjacentNodes_;

    /**
     * @brief Contiguous support used when an ordinary caller did not certify it.
     *
     * The saturated builder supplies its already collected support directly.
     * Other callers retain the self-contained API and pay one local copy, after
     * which all phases of the update reuse the same contiguous pixel sequence.
     */
    std::vector<NodeId> selectedSupportScratch_;

    /**
     * @brief Returns the read-only topology view associated with a weighted tree.
     *
     * @param tree Tree topology used by the operation.
     * @return The read-only topology view associated with a weighted tree.
     */
    static const MorphologicalTree& topologyOf(const tree_t* tree) {
        assert(tree != nullptr);
        return tree->topology();
    }

    /**
     * @brief Returns the current altitude of one tree node.
     *
     * @param tree Tree topology used by the operation.
     * @param nodeId Identifier of the node used by the operation.
     * @return The current altitude of one tree node.
     */
    altitude_t nodeAltitude(const tree_t* tree, NodeId nodeId) const {
        assert(tree != nullptr);
        return tree->getAltitude(nodeId);
    }

    /**
     * @brief Returns the primal tree for an update whose dual tree has polarity `isMaxtree`.
     *
     * @param isMaxtree Flag controlling is maxtree.
     * @return The primal tree for an update whose dual tree has polarity isMaxtree.
     */
    tree_t* getPrimalTree(bool isMaxtree) { return isMaxtree ? mintree_ : maxtree_; }

    /**
     * @brief Records one root whose snapshot subtree must be invalidated.
     *
     * @param nodeId Identifier of the node processed by the operation.
     */
    void noteTopologyChanged(NodeId nodeId) {
        if (nodeId == InvalidNode || nodeId < 0 || static_cast<std::size_t>(nodeId) >= topologyChangedMarks_.n) {
            return;
        }
        const std::size_t index = static_cast<std::size_t>(nodeId);
        if (topologyChangedMarks_.isMarked(index)) {
            return;
        }
        topologyChangedMarks_.mark(index);
        lastTopologyChangedNodes_.push_back(nodeId);
    }

    /**
     * @brief Clears all per-step marks and scratch buffers before a leaf update.
     *
     * @param isMaxtree Flag controlling is maxtree.
     */
    void initializeUpdateStepState(bool isMaxtree) {
        mergeNodesByLevel_.resetCollection(isMaxtree);
        removedMarks_.resetAll();
        removedNodesPendingAbsorption_.clear();
        pixelsInLeafMarks_.resetAll();
        climbedNodeMarks_.resetAll();
        lastCandidateNodes_.clear();
        topologyChangedMarks_.resetAll();
        lastTopologyChangedNodes_.clear();
        lastBulkProperPartTransfer_ = false;
    }

    /**
     * @brief Detaches a live non-root node from its parent without releasing it.
     *
     * This mirrors the efficient editor usage in `DualMinMaxTreeIncrementalFilter`:
     * when the caller only needs a temporary detached subtree, `editor.detach`
     * avoids the parent/child relation check required by `removeChild`.
     *
     * @param topology Tree topology used by the operation.
     * @param editor Active tree editor used by the operation.
     * @param nodeId Identifier of the node used by the operation.
     */
    void detachIfNonRoot(const MorphologicalTree& topology, editor_t& editor, NodeId nodeId) {
        if (nodeId == InvalidNode || !topology.isNode(nodeId) || !topology.isAlive(nodeId) || topology.isRoot(nodeId)) {
            return;
        }
        noteTopologyChanged(nodeId);
        editor.detach(nodeId);
    }

    /**
     * @brief Attaches one subtree and records its changed parent edge.
     *
     * @param editor Mutable editor updated by the operation.
     * @param parentId Identifier of the parent node.
     * @param childId Identifier of the child node.
     */
    void attachChanged(editor_t& editor, NodeId parentId, NodeId childId) {
        noteTopologyChanged(childId);
        editor.attach(parentId, childId);
    }

    /**
     * @brief Detaches a node from its current parent and optionally releases it.
     *
     * Non-release detaches use the lighter `TreeEditor::detach` path. Release
     * operations keep using `removeChild(..., true)` because they must return the
     * now-empty node slot to the tree's free list.
     *
     * @param tree Tree topology used by the operation.
     * @param editor Active tree editor used by the operation.
     * @param nodeId Identifier of the node used by the operation.
     * @param releaseNode Node identifier represented by `releaseNode`.
     */
    void disconnect(tree_t* tree, editor_t& editor, NodeId nodeId, bool releaseNode) {
        assert(tree != nullptr);
        const MorphologicalTree& topology = topologyOf(tree);
        if (!releaseNode) {
            detachIfNonRoot(topology, editor, nodeId);
            return;
        }
        if (topology.isRoot(nodeId)) {
            return;
        }

        const NodeId parentId = topology.getNodeParent(nodeId);
        if (parentId == InvalidNode || parentId == nodeId) {
            return;
        }

        noteTopologyChanged(nodeId);
        editor.removeChild(parentId, nodeId, releaseNode);
    }

    /**
     * @brief Moves every child and every direct proper part from `childId` to `parentId`.
     *
     * After this operation, `childId` is structurally empty but still attached
     * until the caller explicitly disconnects or releases it.
     *
     * @param tree Tree topology used by the operation.
     * @param editor Active tree editor used by the operation.
     * @param parentId Identifier of the parent node.
     * @param childId Identifier of the child node.
     */
    void mergedParentAndChildren(tree_t* tree, editor_t& editor, NodeId parentId, NodeId childId) {
        assert(tree != nullptr);
        const MorphologicalTree& topology = topologyOf(tree);
        for (NodeId child : topology.getChildren(childId)) {
            noteTopologyChanged(child);
        }
        editor.moveChildren(parentId, childId);
        editor.moveProperParts(parentId, childId);
    }

    /**
     * @brief Moves the selected support `C` to the current union node.
     *
     * The support is represented by the direct proper parts of `leafId` in the
     * primal tree. Each pixel is looked up in the dual tree and moved from its
     * current owner to `unionNode`. Donor nodes that become empty are marked for
     * deferred absorption; if the emptied donor is `nodeCa`, the caller must handle
     * the special replacement/root cases after the sweep.
     *
     * @param dualTree Dual tree used by the operation.
     * @param editor Active tree editor used by the operation.
     * @param unionNode Node identifier represented by `unionNode`.
     * @param support Direct proper-part support of the selected primal leaf.
     * @param nodeCa Identifier of the common-ancestor node.
     * @param nodeCaRemoved Whether the common-ancestor node was removed.
     * @param certifiedWholeSupportOwner Optional certified owner of the complete
     *        support, enabling one bulk proper-part transfer.
     */
    void moveSelectedProperPartsToNode(tree_t* dualTree, editor_t& editor, NodeId unionNode, std::span<const NodeId> support, NodeId nodeCa,
                                       bool& nodeCaRemoved, NodeId certifiedWholeSupportOwner) {
        const MorphologicalTree& dualTopology = topologyOf(dualTree);
        if (certifiedWholeSupportOwner != InvalidNode && certifiedWholeSupportOwner != unionNode && dualTopology.isNode(certifiedWholeSupportOwner) &&
            dualTopology.isAlive(certifiedWholeSupportOwner) &&
            dualTopology.getNumProperParts(certifiedWholeSupportOwner) == static_cast<int>(support.size())) {
#ifndef NDEBUG
            for (NodeId pixelId : support) {
                assert(dualTopology.getProperPartOwner(pixelId) == certifiedWholeSupportOwner);
            }
#endif
            editor.moveProperParts(unionNode, certifiedWholeSupportOwner);
            lastBulkProperPartTransfer_ = true;
            if (dualTopology.isAlive(certifiedWholeSupportOwner) && dualTopology.getNumProperParts(certifiedWholeSupportOwner) == 0 &&
                !removedMarks_.isMarked(static_cast<std::size_t>(certifiedWholeSupportOwner))) {
                removedMarks_.mark(static_cast<std::size_t>(certifiedWholeSupportOwner));
                removedNodesPendingAbsorption_.push_back(certifiedWholeSupportOwner);
                if (certifiedWholeSupportOwner == nodeCa) {
                    nodeCaRemoved = true;
                }
            }
            return;
        }
        for (NodeId pixelId : support) {
            const NodeId ownerId = dualTopology.getProperPartOwner(pixelId);
            if (ownerId == InvalidNode || ownerId == unionNode) {
                continue;
            }

            editor.moveProperPart(unionNode, ownerId, pixelId);
            if (dualTopology.isAlive(ownerId) && dualTopology.getNumProperParts(ownerId) == 0 && !removedMarks_.isMarked(static_cast<std::size_t>(ownerId))) {
                removedMarks_.mark(static_cast<std::size_t>(ownerId));
                removedNodesPendingAbsorption_.push_back(ownerId);
                if (ownerId == nodeCa) {
                    nodeCaRemoved = true;
                }
            }
        }
    }

    /**
     * @brief Tests whether a marked empty node is still alive and safe to absorb.
     *
     * @param tree Tree topology used by the operation.
     * @param nodeId Identifier of the node used by the operation.
     * @return True if a marked empty node is still alive and safe to absorb; otherwise false.
     */
    bool canAbsorbRemovedNode(tree_t* tree, NodeId nodeId) const {
        return tree != nullptr && nodeId != InvalidNode && topologyOf(tree).isNode(nodeId) && topologyOf(tree).isAlive(nodeId) &&
               removedMarks_.isMarked(static_cast<std::size_t>(nodeId)) && topologyOf(tree).getNumProperParts(nodeId) == 0;
    }

    /**
     * @brief Contracts an empty root by promoting its altitude-compatible child.
     *
     * If a removed node is currently the root, its children cannot simply be moved
     * to a parent. The child with the extremal altitude for the dual tree polarity
     * becomes the new root and all other children are reattached below it.
     *
     * @param tree Tree topology used by the operation.
     * @param editor Active tree editor used by the operation.
     * @param removedNodeId Node identifier represented by `removedNodeId`.
     */
    void absorbRemovedRootNode(tree_t* tree, editor_t& editor, NodeId removedNodeId) {
        const bool isMaxtree = tree == maxtree_;
        const MorphologicalTree& topology = topologyOf(tree);
        const NodeId firstChild = topology.getFirstChild(removedNodeId);
        if (firstChild == InvalidNode) {
            return;
        }

        NodeId newRoot = firstChild;
        for (NodeId childId = firstChild; childId != InvalidNode; childId = topology.getNextSibling(childId)) {
            if ((isMaxtree && nodeAltitude(tree, childId) < nodeAltitude(tree, newRoot)) ||
                (!isMaxtree && nodeAltitude(tree, childId) > nodeAltitude(tree, newRoot))) {
                newRoot = childId;
            }
        }

        for (NodeId childId = firstChild; childId != InvalidNode;) {
            const NodeId nextChildId = topology.getNextSibling(childId);
            if (childId != newRoot && !topology.hasChild(newRoot, childId)) {
                detachIfNonRoot(topology, editor, childId);
                attachChanged(editor, newRoot, childId);
            }
            childId = nextChildId;
        }

        noteTopologyChanged(newRoot);
        noteTopologyChanged(removedNodeId);
        editor.setRoot(newRoot);
        editor.releaseNode(removedNodeId);
    }

    /**
     * @brief Contracts one empty non-root node into its current parent.
     *
     * The node's children and direct proper parts are first moved to the parent,
     * then the empty node is released. If the parent becomes empty as a result of
     * the merge, the parent is marked so the post-order absorption pass can keep
     * contracting upward.
     *
     * @param tree Tree topology used by the operation.
     * @param editor Active tree editor used by the operation.
     * @param removedNodeId Node identifier represented by `removedNodeId`.
     * @return The surviving parent empty non-root node into its current parent.
     */
    NodeId absorbRemovedNonRootNode(tree_t* tree, editor_t& editor, NodeId removedNodeId) {
        const MorphologicalTree& topology = topologyOf(tree);
        const NodeId parentId = topology.getNodeParent(removedNodeId);
        if (parentId == InvalidNode || parentId == removedNodeId || !topology.isAlive(parentId)) {
            return InvalidNode;
        }

        mergedParentAndChildren(tree, editor, parentId, removedNodeId);
        disconnect(tree, editor, removedNodeId, true);

        if (topology.isAlive(parentId) && topology.getNumProperParts(parentId) == 0) {
            removedMarks_.mark(static_cast<std::size_t>(parentId));
            return parentId;
        }

        return InvalidNode;
    }

    /**
     * @brief Contracts all empty nodes accumulated during the current update.
     *
     * The pass is iterative post-order. Children are considered before the marked
     * node itself so moving children/proper parts never leaves dangling links in
     * the connected component. Newly emptied parents are pushed back on the stack.
     *
     * @param tree Tree topology used by the operation.
     * @param editor Active tree editor used by the operation.
     * @param removedNodeIds Node identifier represented by `removedNodeIds`.
     */
    void absorbRemovedNodes(tree_t* tree, editor_t& editor, const std::vector<NodeId>& removedNodeIds) {
        struct Frame {
            NodeId nodeId = InvalidNode;
            NodeId nextChildId = InvalidNode;
        };

        if (tree == nullptr || removedNodeIds.empty()) {
            return;
        }

        std::vector<Frame> stack;
        stack.reserve(std::max<std::size_t>(64, removedNodeIds.size()));
        const auto makeFrame = [tree](NodeId nodeId) {
            if (nodeId == InvalidNode || !topologyOf(tree).isNode(nodeId) || !topologyOf(tree).isAlive(nodeId)) {
                return Frame{nodeId, InvalidNode};
            }
            return Frame{nodeId, topologyOf(tree).getFirstChild(nodeId)};
        };

        for (auto it = removedNodeIds.rbegin(); it != removedNodeIds.rend(); ++it) {
            if (*it != InvalidNode) {
                stack.push_back(makeFrame(*it));
            }
        }

        while (!stack.empty()) {
            Frame& frame = stack.back();
            NodeId currentNodeId = InvalidNode;

            if (!canAbsorbRemovedNode(tree, frame.nodeId)) {
                stack.pop_back();
            } else if (frame.nextChildId != InvalidNode) {
                const NodeId childId = frame.nextChildId;
                frame.nextChildId = topologyOf(tree).getNextSibling(childId);
                stack.push_back(makeFrame(childId));
            } else {
                currentNodeId = frame.nodeId;
                stack.pop_back();
            }

            if (canAbsorbRemovedNode(tree, currentNodeId)) {
                if (topologyOf(tree).isRoot(currentNodeId)) {
                    absorbRemovedRootNode(tree, editor, currentNodeId);
                } else {
                    const NodeId parentId = absorbRemovedNonRootNode(tree, editor, currentNodeId);
                    if (parentId != InvalidNode) {
                        stack.push_back(makeFrame(parentId));
                    }
                }
            }
        }
    }

    /**
     * @brief Collects a subtree in post-order so leaf pruning remains the elementary step.
     *
     * @param tree Tree topology used by the operation.
     * @param nodeId Identifier of the node used by the operation.
     * @param out Destination receiving the result.
     */
    static void collectPostOrderNodes(const tree_t& tree, NodeId nodeId, std::vector<NodeId>& out) {
        for (NodeId childId : tree.topology().getChildren(nodeId)) {
            collectPostOrderNodes(tree, childId, out);
        }
        out.push_back(nodeId);
    }

  public:
    /**
     * @brief Reports whether this instantiation uses dense per-level buckets.
     *
     * @return True if this instantiation uses dense per-level buckets; otherwise false.
     */
    static constexpr bool usesDenseLevelBackend() { return use_dense_levels; }

    /**
     * @brief Returns the bit-width threshold used to choose the dense bucket backend.
     *
     * @return The bit-width threshold used to choose the dense bucket backend.
     */
    static constexpr int denseLevelBackendMaxBits() { return MMCFILTERS_COMPONENT_TREE_ADJUSTMENT_DENSE_MAX_BITS; }

    /**
     * @brief Builds the leaf adjuster over a fixed weighted min-tree / max-tree pair.
     *
     * The adjuster does not own the trees or the adjacency relation. Public prune
     * methods always operate on this same pair; they do not rebind trees per call.
     * Both trees must represent the same image domain and must use node ids from
     * this repository's internal dense node-id space.
     *
     * @param minTree Min-tree used by the operation.
     * @param maxTree Max-tree used by the operation.
     * @param graph Region or adjacency graph used by the operation.
     */
    DualMinMaxTreeIncrementalFilterLeaf(tree_t* minTree, tree_t* maxTree, const RegularGridAdjacency2D& graph)
        : mintree_(minTree), maxtree_(maxTree), graph_(&graph), mergeNodesByLevel_(std::max(minTree ? minTree->topology().getNumInternalNodeSlots() : 0,
                                                                                            maxTree ? maxTree->topology().getNumInternalNodeSlots() : 0)),
          removedMarks_(static_cast<std::size_t>(
              std::max(minTree ? minTree->topology().getNumInternalNodeSlots() : 0, maxTree ? maxTree->topology().getNumInternalNodeSlots() : 0))),
          pixelsInLeafMarks_(static_cast<std::size_t>(
              std::max(minTree ? minTree->topology().getNumTotalProperParts() : 0, maxTree ? maxTree->topology().getNumTotalProperParts() : 0))),
          climbedNodeMarks_(static_cast<std::size_t>(
              std::max(minTree ? minTree->topology().getNumInternalNodeSlots() : 0, maxTree ? maxTree->topology().getNumInternalNodeSlots() : 0))),
          topologyChangedMarks_(static_cast<std::size_t>(
              std::max(minTree ? minTree->topology().getNumInternalNodeSlots() : 0, maxTree ? maxTree->topology().getNumInternalNodeSlots() : 0))) {
        assert(mintree_ != nullptr);
        assert(maxtree_ != nullptr);
        assert(graph_ != nullptr);
        assert(mintree_->topology().getNumRowsOfGridDomain2D() == maxtree_->topology().getNumRowsOfGridDomain2D());
        assert(mintree_->topology().getNumColsOfGridDomain2D() == maxtree_->topology().getNumColsOfGridDomain2D());
    }

    /**
     * @brief Exposes adjacent seed nodes collected for the current update.
     *
     * @return Reference to the resulting object.
     */
    std::vector<NodeId>& getAdjacentNodes() { return mergeNodesByLevel_.getAdjacentNodes(); }

    /**
     * @brief Exposes frontier nodes outside the merge interval for the current update.
     *
     * @return Reference to the resulting object.
     */
    std::vector<NodeId>& getFrontierNodesAboveB() { return mergeNodesByLevel_.getFrontierNodesAboveB(); }

    /**
     * @brief Exposes the merge bucket associated with an altitude level.
     *
     * @param level Altitude level used by the operation.
     * @return Reference to the resulting object.
     */
    std::vector<NodeId>& getMergedNodes(const altitude_t& level) { return mergeNodesByLevel_.getMergedNodes(level); }

    /**
     * @brief Returns dual-tree nodes that may affect the caller's candidate agenda.
     *
     * SDRT construction uses this after each leaf prune to refresh only affected
     * dual candidates. The vector is cleared at the beginning of each update and
     * receives the union node created at every consumed merge level, matching the
     * MorphoTreeAdjust leaf-builder behavior.
     *
     * @return Dual-tree nodes that may affect the caller's candidate agenda.
     */
    const std::vector<NodeId>& getLastCandidateNodes() const { return lastCandidateNodes_; }

    /**
     * @brief Returns a safe superset of changed parent-edge roots.
     *
     * Each reported node either changed its parent relation, was released, or
     * roots a subtree whose direct attachment changed in the last dual update.
     *
     * @return The last topology changed nodes.
     */
    const std::vector<NodeId>& getLastTopologyChangedNodes() const { return lastTopologyChangedNodes_; }

    /**
     * @brief Reports whether the last selected proper part used a bulk transfer.
     *
     * @return `true` when the last selected proper part used a bulk transfer.
     */
    [[nodiscard]] bool getLastBulkProperPartTransfer() const noexcept { return lastBulkProperPartTransfer_; }

    /**
     * @brief Reuses a certified max-leaf boundary in the next min-tree update.
     *
     * `adjacentNodes` must contain every distinct current min-tree owner
     * touching the external pixel boundary of `leafId`. The cache is consumed
     * by the next matching update and ignored by all other calls.
     *
     * @param leafId Leaf identifier used by the operation.
     * @param adjacentNodes Adjacent nodes used by the operation.
     */
    void primeMaxLeafBoundary(NodeId leafId, std::span<const NodeId> adjacentNodes) {
        primedDualTree_ = mintree_;
        primedLeafId_ = leafId;
        primedAdjacentNodes_.assign(adjacentNodes.begin(), adjacentNodes.end());
    }

    /**
     * @brief Reuses a certified min-leaf boundary in the next max-tree update.
     *
     * `adjacentNodes` must contain every distinct current max-tree owner
     * touching the external pixel boundary of `leafId`.
     *
     * @param leafId Leaf identifier used by the operation.
     * @param adjacentNodes Adjacent nodes used by the operation.
     */
    void primeMinLeafBoundary(NodeId leafId, std::span<const NodeId> adjacentNodes) {
        primedDualTree_ = maxtree_;
        primedLeafId_ = leafId;
        primedAdjacentNodes_.assign(adjacentNodes.begin(), adjacentNodes.end());
    }

    /**
     * @brief Builds merge buckets and frontiers for one removed leaf.
     *
     * Starting from the leaf support `C` in `primalTree`, the method marks all
     * pixels in `C`, scans their graph neighbors, maps each external neighbor to
     * its owner in `dualTree`, and climbs ancestor paths that belong to the
     * altitude interval governed by `nodeCa` and `b`.
     *
     * Nodes whose level is inside the merge side of the interval are inserted in
     * `mergeNodesByLevel_`. Boundary nodes outside the merge side are recorded as
     * frontiers and later attached to the union node created at level `b`.
     *
     * @param dualTree Dual tree used by the operation.
     * @param support Direct proper-part support of the selected primal leaf.
     * @param certifiedAdjacentNodes Exact distinct dual owners on the external
     *        boundary. When non-empty, the method does not rescan or mark the
     *        support to rediscover that boundary.
     * @param nodeCa Identifier of the common-ancestor node.
     * @param b Branch or side selector used by the operation.
     * @param isMaxtree Flag controlling is maxtree.
     */
    void buildMergedAndNestedCollections(tree_t* dualTree, std::span<const NodeId> support, std::span<const NodeId> certifiedAdjacentNodes, NodeId nodeCa,
                                         altitude_t b, bool isMaxtree) {
        assert(dualTree != nullptr);
        const MorphologicalTree& dualTopology = topologyOf(dualTree);
        mergeNodesByLevel_.resetCollection(isMaxtree);
        assert(nodeCa != InvalidNode);
        const altitude_t altitudeCa = nodeAltitude(dualTree, nodeCa);
        pixelsInLeafMarks_.resetAll();
        climbedNodeMarks_.resetAll();

        const auto addAdjacentNode = [&](NodeId nodeQ) {
            if (nodeQ == InvalidNode || !dualTopology.isAlive(nodeQ)) {
                return;
            }
            const altitude_t altitudeQ = nodeAltitude(dualTree, nodeQ);
            const bool validSeed = (isMaxtree && altitudeQ >= altitudeCa) || (!isMaxtree && altitudeQ <= altitudeCa);
            if (!validSeed || !mergeNodesByLevel_.markAdjacentSeed(nodeQ)) {
                return;
            }

            NodeId nodeSubtree = nodeQ;
            NodeId n = nodeQ;
            while (n != InvalidNode && dualTopology.isAlive(n) && !climbedNodeMarks_.isMarked(static_cast<std::size_t>(n))) {
                const altitude_t levelCurrent = nodeAltitude(dualTree, n);
                if (!((isMaxtree && levelCurrent >= altitudeCa) || (!isMaxtree && levelCurrent <= altitudeCa))) {
                    break;
                }

                climbedNodeMarks_.mark(static_cast<std::size_t>(n));
                nodeSubtree = n;

                if ((isMaxtree && levelCurrent <= b) || (!isMaxtree && levelCurrent >= b)) {
                    mergeNodesByLevel_.addMergeNode(*dualTree, nodeSubtree);
                } else {
                    NodeId parentId = dualTopology.getNodeParent(nodeSubtree);
                    if (parentId == nodeSubtree) {
                        parentId = InvalidNode;
                    }
                    if (!(parentId != InvalidNode &&
                          ((isMaxtree && nodeAltitude(dualTree, parentId) > b) || (!isMaxtree && nodeAltitude(dualTree, parentId) < b)))) {
                        mergeNodesByLevel_.addFrontierNodeAboveB(nodeSubtree);
                    }
                }

                const NodeId parentId = dualTopology.getNodeParent(n);
                if (parentId == n) {
                    break;
                }
                n = parentId;
            }
        };

        if (!certifiedAdjacentNodes.empty()) {
            for (NodeId nodeQ : certifiedAdjacentNodes) {
                addAdjacentNode(nodeQ);
            }
        } else {
            for (NodeId pixelId : support) {
                pixelsInLeafMarks_.mark(static_cast<std::size_t>(pixelId));
            }
            for (NodeId p : support) {
                for (NodeId q : graph_->getNeighborIndices(p)) {
                    if (pixelsInLeafMarks_.isMarked(static_cast<std::size_t>(q))) {
                        continue;
                    }
                    addAdjacentNode(dualTopology.getProperPartOwner(q));
                }
            }
        }
    }

    /**
     * @brief Updates the dual tree after one leaf is removed from the primal tree.
     *
     * This is the algorithmic core of leaf mode. The method:
     *
     * - identifies the parent level `b` of the primal leaf;
     * - finds `nodeCa`, the extremal dual component containing the leaf support;
     * - builds merge buckets/frontiers for the affected dual region;
     * - sweeps altitude levels toward `nodeCa`, merging same-level nodes;
     * - moves the leaf support to the union node at level `b`;
     * - reattaches frontiers and previous-level union nodes;
     * - handles the cases where `nodeCa` becomes empty or the root changes;
     * - contracts all empty donor nodes collected during the update.
     *
     * The caller is responsible for pruning `leafId` from the primal tree after
     * this dual update completes.
     *
     * @param dualTree Dual tree used by the operation.
     * @param leafId Identifier of the leaf node.
     * @param certifiedSupport Optional support already collected by a preceding
     *        certificate.
     * @param certifiedNodeCa Optional exact extremal dual owner of that support.
     * @param certifiedAdjacentNodes Optional exact distinct owners touching the
     *        external boundary.
     * @param certifiedWholeSupportOwner Optional certified owner of the complete
     *        support, enabling one bulk proper-part transfer.
     */
    void updateTree(tree_t* dualTree, NodeId leafId, std::span<const NodeId> certifiedSupport = {}, NodeId certifiedNodeCa = InvalidNode,
                    std::span<const NodeId> certifiedAdjacentNodes = {}, NodeId certifiedWholeSupportOwner = InvalidNode) {
        assert(dualTree != nullptr);
        assert(leafId != InvalidNode);

        const bool isMaxtree = dualTree == maxtree_;
        initializeUpdateStepState(isMaxtree);
        tree_t* primalTree = getPrimalTree(isMaxtree);
        assert(primalTree != nullptr);
        const MorphologicalTree& dualTopology = topologyOf(dualTree);
        const MorphologicalTree& primalTopology = topologyOf(primalTree);

        assert(primalTopology.isNode(leafId) && primalTopology.isAlive(leafId));
        assert(primalTopology.isLeaf(leafId));

        const NodeId leafParent = primalTopology.getNodeParent(leafId);
        assert(leafParent != InvalidNode && leafParent != leafId);
        const altitude_t b = nodeAltitude(primalTree, leafParent);

        assert(primalTopology.getNumProperParts(leafId) > 0);
        std::span<const NodeId> support = certifiedSupport;
        if (support.empty()) {
            const auto properParts = primalTopology.getProperParts(leafId);
            selectedSupportScratch_.assign(properParts.begin(), properParts.end());
            support = selectedSupportScratch_;
        }
        assert(support.size() == static_cast<std::size_t>(primalTopology.getNumProperParts(leafId)));

        NodeId nodeCa = certifiedNodeCa;
        altitude_t altitudeCa = altitude_t{};
        if (nodeCa != InvalidNode) {
            assert(dualTopology.isNode(nodeCa) && dualTopology.isAlive(nodeCa));
            altitudeCa = nodeAltitude(dualTree, nodeCa);
        } else {
            for (NodeId pixelId : support) {
                const NodeId ownerId = dualTopology.getProperPartOwner(pixelId);
                if (ownerId == InvalidNode) {
                    continue;
                }

                const altitude_t ownerAltitude = nodeAltitude(dualTree, ownerId);
                if (nodeCa == InvalidNode || (isMaxtree && ownerAltitude < altitudeCa) || (!isMaxtree && ownerAltitude > altitudeCa)) {
                    nodeCa = ownerId;
                    altitudeCa = ownerAltitude;
                }
            }
        }
        assert(nodeCa != InvalidNode);

        std::span<const NodeId> adjacentNodes = certifiedAdjacentNodes;
        if (adjacentNodes.empty() && primedDualTree_ == dualTree && primedLeafId_ == leafId) {
            adjacentNodes = primedAdjacentNodes_;
        }
        buildMergedAndNestedCollections(dualTree, support, adjacentNodes, nodeCa, b, isMaxtree);
        primedDualTree_ = nullptr;
        primedLeafId_ = InvalidNode;
        primedAdjacentNodes_.clear();

        editor_t editor = mmcfilters::detail::beginEstablishedWeightedEdit(*dualTree);
        altitude_t mergeLevel = mergeNodesByLevel_.firstMergeLevel();
        NodeId unionNode = InvalidNode;
        NodeId previousUnionNode = InvalidNode;
        bool nodeCaRemoved = false;

        while (mergeNodesByLevel_.hasMergeLevel() && ((isMaxtree && mergeLevel > altitudeCa) || (!isMaxtree && mergeLevel < altitudeCa))) {
            auto& mergeNodesAtLevel = mergeNodesByLevel_.getMergedNodes(mergeLevel);
            unionNode = InvalidNode;

            for (auto nodeId : mergeNodesAtLevel) {
                if (!dualTopology.isAlive(nodeId)) {
                    continue;
                }
                if (unionNode == InvalidNode || nodeId < unionNode) {
                    unionNode = nodeId;
                }
            }

            for (auto nodeId : mergeNodesAtLevel) {
                if (!dualTopology.isAlive(nodeId)) {
                    continue;
                }

                if (nodeId == unionNode) {
                    disconnect(dualTree, editor, unionNode, false);
                    continue;
                }

                if (unionNode == InvalidNode) {
                    continue;
                }

                mergedParentAndChildren(dualTree, editor, unionNode, nodeId);
                disconnect(dualTree, editor, nodeId, true);
            }

            if (unionNode == InvalidNode) {
                mergeLevel = mergeNodesByLevel_.nextMergeLevel();
                unionNode = previousUnionNode;
                continue;
            }

            lastCandidateNodes_.push_back(unionNode);

            if (mergeLevel == b) {
                moveSelectedProperPartsToNode(dualTree, editor, unionNode, support, nodeCa, nodeCaRemoved, certifiedWholeSupportOwner);

                for (auto nodeId : mergeNodesByLevel_.getFrontierNodesAboveB()) {
                    if (!dualTopology.isAlive(nodeId)) {
                        continue;
                    }
                    disconnect(dualTree, editor, nodeId, false);
                    attachChanged(editor, unionNode, nodeId);
                }
            }

            if (previousUnionNode != InvalidNode && dualTopology.isAlive(previousUnionNode) && !dualTopology.hasChild(unionNode, previousUnionNode)) {
                if (!dualTopology.isRoot(previousUnionNode)) {
                    disconnect(dualTree, editor, previousUnionNode, false);
                }
                attachChanged(editor, unionNode, previousUnionNode);
            }

            previousUnionNode = unionNode;
            mergeLevel = mergeNodesByLevel_.nextMergeLevel();
        }

        const NodeId finalUnionNode = previousUnionNode;
        if (finalUnionNode == InvalidNode || !dualTopology.isAlive(finalUnionNode)) {
            auto proof = editor.proveIncremental();
            editor.commit(std::move(proof));
            return;
        }

        if (nodeCaRemoved) {
            if (!dualTopology.isRoot(nodeCa)) {
                const NodeId parentIdNodeCa = dualTopology.getNodeParent(nodeCa);
                if (!dualTopology.isRoot(finalUnionNode)) {
                    disconnect(dualTree, editor, finalUnionNode, false);
                }
                attachChanged(editor, parentIdNodeCa, finalUnionNode);

                for (NodeId n = dualTopology.getFirstChild(nodeCa); n != InvalidNode;) {
                    const NodeId next = dualTopology.getNextSibling(n);
                    if (n != finalUnionNode && !dualTopology.hasChild(finalUnionNode, n)) {
                        detachIfNonRoot(dualTopology, editor, n);
                        attachChanged(editor, finalUnionNode, n);
                    }
                    n = next;
                }

                disconnect(dualTree, editor, nodeCa, true);
            } else {
                NodeId newRoot = finalUnionNode;
                for (NodeId n = dualTopology.getFirstChild(nodeCa); n != InvalidNode; n = dualTopology.getNextSibling(n)) {
                    if ((isMaxtree && nodeAltitude(dualTree, n) < nodeAltitude(dualTree, newRoot)) ||
                        (!isMaxtree && nodeAltitude(dualTree, n) > nodeAltitude(dualTree, newRoot))) {
                        newRoot = n;
                    }
                }

                if (newRoot != finalUnionNode) {
                    if (!dualTopology.isRoot(finalUnionNode)) {
                        disconnect(dualTree, editor, finalUnionNode, false);
                    }
                    attachChanged(editor, newRoot, finalUnionNode);
                }

                for (NodeId n = dualTopology.getFirstChild(nodeCa); n != InvalidNode;) {
                    const NodeId next = dualTopology.getNextSibling(n);
                    if (n != newRoot && !dualTopology.hasChild(newRoot, n)) {
                        detachIfNonRoot(dualTopology, editor, n);
                        attachChanged(editor, newRoot, n);
                    }
                    n = next;
                }

                noteTopologyChanged(newRoot);
                noteTopologyChanged(nodeCa);
                editor.setRoot(newRoot);
                editor.releaseNode(nodeCa);
            }
        } else if (finalUnionNode != nodeCa) {
            if (!dualTopology.isRoot(finalUnionNode)) {
                disconnect(dualTree, editor, finalUnionNode, false);
            }
            attachChanged(editor, nodeCa, finalUnionNode);
        }

        absorbRemovedNodes(dualTree, editor, removedNodesPendingAbsorption_);
        auto proof = editor.proveIncremental();
        editor.commit(std::move(proof));
    }

    /**
     * @brief Prunes max-tree subtrees and updates this adjuster's min-tree leaf by leaf.
     *
     * Each requested subtree is first converted to post-order so children are
     * processed before their parent. Nodes that have already disappeared during
     * previous leaf removals are skipped.
     *
     * @param nodesToPrune Node identifiers selected for pruning.
     */
    void pruneMaxTreeAndUpdateMinTree(std::vector<NodeId>& nodesToPrune) {
        std::vector<NodeId> postOrderNodes;
        for (NodeId rootSubtree : nodesToPrune) {
            if (rootSubtree == InvalidNode || rootSubtree == maxtree_->topology().getRoot() || !maxtree_->topology().isNode(rootSubtree) ||
                !maxtree_->topology().isAlive(rootSubtree)) {
                continue;
            }

            postOrderNodes.clear();
            collectPostOrderNodes(*maxtree_, rootSubtree, postOrderNodes);
            for (NodeId leafId : postOrderNodes) {
                if (leafId == InvalidNode || leafId == maxtree_->topology().getRoot() || !maxtree_->topology().isNode(leafId) ||
                    !maxtree_->topology().isAlive(leafId)) {
                    continue;
                }
                assert(maxtree_->topology().isLeaf(leafId));
                updateTree(mintree_, leafId);
                maxtree_->pruneNode(leafId);
            }
        }
    }

    /**
     * @brief Prunes one max-tree leaf and updates this adjuster's min-tree.
     *
     * @param leafId Identifier of the leaf node.
     */
    void pruneMaxLeafAndUpdateMinTree(NodeId leafId) {
        if (leafId == InvalidNode || leafId == maxtree_->topology().getRoot() || !maxtree_->topology().isNode(leafId) ||
            !maxtree_->topology().isAlive(leafId)) {
            return;
        }
        assert(maxtree_->topology().isLeaf(leafId));
        updateTree(mintree_, leafId);
        maxtree_->pruneNode(leafId);
    }

    /**
     * @brief Prunes one certified max-tree leaf without rediscovering its workset.
     *
     * The three certificate fields must describe the current leaf and min-tree
     * state. They are consumed synchronously before the primal leaf is pruned.
     *
     * @param leafId Leaf identifier used by the operation.
     * @param certifiedSupport Certified support used by the operation.
     * @param certifiedNodeCa Certified node ca used by the operation.
     * @param certifiedAdjacentNodes Certified adjacent nodes used by the operation.
     * @param certifiedWholeSupportOwner Certified whole support owner used by the operation.
     */
    void pruneMaxLeafAndUpdateMinTree(NodeId leafId, std::span<const NodeId> certifiedSupport, NodeId certifiedNodeCa,
                                      std::span<const NodeId> certifiedAdjacentNodes, NodeId certifiedWholeSupportOwner = InvalidNode) {
        if (leafId == InvalidNode || leafId == maxtree_->topology().getRoot() || !maxtree_->topology().isNode(leafId) ||
            !maxtree_->topology().isAlive(leafId)) {
            return;
        }
        assert(maxtree_->topology().isLeaf(leafId));
        updateTree(mintree_, leafId, certifiedSupport, certifiedNodeCa, certifiedAdjacentNodes, certifiedWholeSupportOwner);
        maxtree_->pruneNode(leafId);
    }

    /**
     * @brief Prunes min-tree subtrees and updates this adjuster's max-tree leaf by leaf.
     *
     * This is the symmetric operation of `pruneMaxTreeAndUpdateMinTree`.
     *
     * @param nodesToPrune Node identifiers selected for pruning.
     */
    void pruneMinTreeAndUpdateMaxTree(std::vector<NodeId>& nodesToPrune) {
        std::vector<NodeId> postOrderNodes;
        for (NodeId rootSubtree : nodesToPrune) {
            if (rootSubtree == InvalidNode || rootSubtree == mintree_->topology().getRoot() || !mintree_->topology().isNode(rootSubtree) ||
                !mintree_->topology().isAlive(rootSubtree)) {
                continue;
            }

            postOrderNodes.clear();
            collectPostOrderNodes(*mintree_, rootSubtree, postOrderNodes);
            for (NodeId leafId : postOrderNodes) {
                if (leafId == InvalidNode || leafId == mintree_->topology().getRoot() || !mintree_->topology().isNode(leafId) ||
                    !mintree_->topology().isAlive(leafId)) {
                    continue;
                }
                assert(mintree_->topology().isLeaf(leafId));
                updateTree(maxtree_, leafId);
                mintree_->pruneNode(leafId);
            }
        }
    }

    /**
     * @brief Prunes one min-tree leaf and updates this adjuster's max-tree.
     *
     * @param leafId Identifier of the leaf node.
     */
    void pruneMinLeafAndUpdateMaxTree(NodeId leafId) {
        if (leafId == InvalidNode || leafId == mintree_->topology().getRoot() || !mintree_->topology().isNode(leafId) ||
            !mintree_->topology().isAlive(leafId)) {
            return;
        }
        assert(mintree_->topology().isLeaf(leafId));
        updateTree(maxtree_, leafId);
        mintree_->pruneNode(leafId);
    }

    /**
     * @brief Prunes one certified min-tree leaf without rediscovering its workset.
     *
     * @param leafId Leaf identifier used by the operation.
     * @param certifiedSupport Certified support used by the operation.
     * @param certifiedNodeCa Certified node ca used by the operation.
     * @param certifiedAdjacentNodes Certified adjacent nodes used by the operation.
     * @param certifiedWholeSupportOwner Certified whole support owner used by the operation.
     */
    void pruneMinLeafAndUpdateMaxTree(NodeId leafId, std::span<const NodeId> certifiedSupport, NodeId certifiedNodeCa,
                                      std::span<const NodeId> certifiedAdjacentNodes, NodeId certifiedWholeSupportOwner = InvalidNode) {
        if (leafId == InvalidNode || leafId == mintree_->topology().getRoot() || !mintree_->topology().isNode(leafId) ||
            !mintree_->topology().isAlive(leafId)) {
            return;
        }
        assert(mintree_->topology().isLeaf(leafId));
        updateTree(maxtree_, leafId, certifiedSupport, certifiedNodeCa, certifiedAdjacentNodes, certifiedWholeSupportOwner);
        mintree_->pruneNode(leafId);
    }

    /**
     * @brief Returns the min-tree associated with this adjuster.
     *
     * @return The min-tree associated with this adjuster.
     */
    tree_t* getMinTree() const { return mintree_; }

    /**
     * @brief Returns the max-tree associated with this adjuster.
     *
     * @return The max-tree associated with this adjuster.
     */
    tree_t* getMaxTree() const { return maxtree_; }
};

} // namespace mmcfilters::adjust
