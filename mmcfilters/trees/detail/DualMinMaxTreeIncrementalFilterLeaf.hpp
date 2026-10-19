#pragma once

#include "../ValuedMorphologicalTree.hpp"
#include "TreeKernelAccess.hpp"
#include "../../utils/CommittedGridAccess.hpp"
#include "../../utils/GenerationStampSet.hpp"

#include <algorithm>
#include <array>
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

namespace mmcfilters::detail {

/**
 * @brief Internal component-tree leaf updater for ExtremumResidualTreeByCTsBuilder.
 *
 * Adapted from the MorphoTreeAdjust DualMinMaxTreeIncrementalFilterLeaf.
 * One operation levels a live non-root primal leaf, updates the dual component
 * tree to represent the same image, and contracts the primal leaf. The builder
 * supplies the nonempty leaf support and the unique dual node whose proper part
 * contains that support.
 * An optional whole-part certificate enables a bulk proper-part transfer.
 *
 * Boundary contacts seed dual ancestor paths. Nodes in the affected altitude
 * interval are grouped into buckets, merged level by level, and reattached at
 * the primal parent level. Empty donor nodes are then contracted. Dense buckets
 * serve small integral altitude types; sparse buckets support wider and
 * floating-point types without quantization.
 *
 * The builder owns both trees and the adjacency and establishes their matching
 * domains and all node/pixel preconditions. This class borrows them and does not
 * validate caller input. Tests of liveness during the sweep skip nodes consumed
 * by earlier merges. Reads follow the same path in CHECKED and UNCHECKED.
 * Established valued-tree edit sessions publish through the generic move-only
 * proof protocol; assertion-enabled builds retain the complete edit oracle.
 *
 * The reported candidate nodes refresh the builder agenda. Proper-part
 * transfers and the final dual node update its residual-support associations.
 */
template <AltitudeValue altitude_t> class DualMinMaxTreeIncrementalFilterLeaf {
  public:
    /** One complete dual proper-part transfer, in execution order. */
    struct ProperPartTransfer {
        NodeId source; ///< Dual node whose complete proper part was transferred.
        NodeId target; ///< Dual node receiving that proper part.
    };

  private:
    /** @brief Defines the `tree_t` alias used by the component. */
    using tree_t = ValuedMorphologicalTree<altitude_t>;
    /** @brief Defines the `editor_t` alias used by the component. */
    using editor_t = ValuedMorphologicalTreeEditor<altitude_t>;
    /** Direct reads after the builder established node and pixel domains. */
    using TreeAccess = mmcfilters::detail::TreeKernelAccess;

    /**
     * @brief Decides whether altitude buckets can be represented by a dense array.
     *
     * MorphoTreeAdjust assumes an 8-bit non-negative altitude domain and uses one
     * vector per possible level. In mmcfilters, the altitude type can vary, so this helper
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
     * This mirrors `MergedNodesCollection` in MorphoTreeAdjust: during
     * one leaf update, adjacent seed components are climbed toward `nodeCa`, and
     * the nodes that must be merged are placed in altitude buckets. The sweep then
     * consumes those buckets in descending order for max-tree updates and ascending
     * order for min-tree updates.
     *
     * The collection also stores frontier nodes outside the `[a, b]` interval.
     * Those frontier nodes are not merged, but must be reattached to the union node
     * created at level `b` so the dual tree remains connected.
     * The ancestor traversal deduplicates nodes before inserting them here.
     */
    class MergedNodesCollection {
      private:
        /**
         * @brief Selects dense or sparse bucket storage for the configured altitude type.
         */
        template <bool dense, typename AltitudeType> struct StorageSelector;

        /**
         * @brief Dense storage: one bucket for each representable altitude value.
         */
        template <typename AltitudeType> struct StorageSelector<true, AltitudeType> {
            /** @brief Domain size. */
            static constexpr std::size_t domain_size = static_cast<std::size_t>(static_cast<long long>(std::numeric_limits<AltitudeType>::max()) -
                                                                                static_cast<long long>(std::numeric_limits<AltitudeType>::lowest()) + 1);
            /** @brief Defines the `type` alias used by the component. */
            using type = std::array<std::vector<NodeId>, domain_size>;
        };

        /**
         * @brief Sparse storage used when the altitude domain is too large to index directly.
         */
        template <typename AltitudeType> struct StorageSelector<false, AltitudeType> {
            /** @brief Defines the `type` alias used by the component. */
            using type = std::map<AltitudeType, std::vector<NodeId>>;
        };

        /** @brief Defines the `storage_t` alias used by the component. */
        using storage_t = typename StorageSelector<use_dense_levels, altitude_t>::type;

        /** @brief Merge buckets keyed by altitude level. */
        storage_t mergeNodesByLevelStorage_;

        /** @brief Levels that currently have non-empty merge buckets in this update. */
        std::vector<altitude_t> mergeLevels_;

        /** @brief Boundary nodes outside the merge interval that must be reattached at `b`. */
        std::vector<NodeId> frontierNodesAboveB_;

        /** @brief Cursor into `mergeLevels_` while the level sweep is running. */
        int currentMergeLevelIndex_ = 0;

        /** @brief Direction flag: max-tree sweeps high-to-low, min-tree sweeps low-to-high. */
        bool isMaxtree_ = false;

      public:
        /**
         * @brief Converts an altitude value into an array index for dense storage.
         *
         * @param level Altitude level.
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
         * @param index Zero-based index.
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
            frontierNodesAboveB_.clear();
            currentMergeLevelIndex_ = 0;
        }

        /**
         * @brief Returns the merge bucket for `level`, creating it for sparse storage.
         *
         * @param level Altitude level.
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
         * @brief Exposes frontier nodes outside the merge interval.
         *
         * @return Mutable reference to the updated object.
         */
        std::vector<NodeId>& getFrontierNodesAboveB() { return frontierNodesAboveB_; }

        /**
         * @brief Adds a frontier node already deduplicated by the ancestor traversal.
         *
         * @param nodeId Dense internal node identifier.
         */
        void addFrontierNodeAboveB(NodeId nodeId) { frontierNodesAboveB_.push_back(nodeId); }

        /**
         * @brief Adds `nodeId` to the merge bucket matching its current altitude.
         *
         * The ancestor traversal has already deduplicated `nodeId`.
         *
         * @param tree Tree topology.
         * @param nodeId Dense internal node identifier.
         */
        void addMergeNode(const tree_t& tree, NodeId nodeId) {
            auto& bucket = getMergedNodes(tree.nodeAltitudes()[static_cast<std::size_t>(nodeId)]);
            bucket.push_back(nodeId);
        }

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

    /** @brief Deduplicates adjacency seeds and ancestors before collection insertion. */
    GenerationStampSet climbedNodeMarks_;

    /** @brief Dual-tree nodes that became candidate-relevant in the last update. */
    std::vector<NodeId> lastCandidateNodes_;

    /** Complete proper-part transfers needed to update support associations. */
    std::vector<ProperPartTransfer> lastProperPartTransfers_;

    /** Dual proper part containing the complete selected support after the last update. */
    NodeId lastSelectedSupportNode_ = InvalidNode;

    /**
     * @brief Returns the read-only topology view associated with a valued tree.
     *
     * @param tree Tree topology.
     * @return The read-only topology view associated with a valued tree.
     */
    static const MorphologicalTree& topologyOf(const tree_t* tree) {
        return tree->topology();
    }

    /**
     * @brief Returns the current altitude of one tree node.
     *
     * @param tree Tree topology.
     * @param nodeId Dense internal node identifier.
     * @return The current altitude of one tree node.
     */
    static altitude_t nodeAltitude(const tree_t* tree, NodeId nodeId) noexcept {
        return tree->nodeAltitudes()[static_cast<std::size_t>(nodeId)];
    }

    /**
     * @brief Returns the primal tree for an update whose dual tree has polarity `isMaxtree`.
     *
     * @param isMaxtree Flag controlling is maxtree.
     * @return The primal tree for an update whose dual tree has polarity isMaxtree.
     */
    tree_t* getPrimalTree(bool isMaxtree) { return isMaxtree ? mintree_ : maxtree_; }

    /**
     * @brief Clears all per-step marks and scratch buffers before a leaf update.
     *
     * @param isMaxtree Flag controlling is maxtree.
     */
    void initializeUpdateStepState(bool isMaxtree) {
        mergeNodesByLevel_.resetCollection(isMaxtree);
        removedMarks_.resetAll();
        removedNodesPendingAbsorption_.clear();
        climbedNodeMarks_.resetAll();
        lastCandidateNodes_.clear();
        lastProperPartTransfers_.clear();
        lastSelectedSupportNode_ = InvalidNode;
    }

    /**
     * @brief Detaches a live non-root node from its parent without releasing it.
     *
     * This mirrors the efficient editor usage in `DualMinMaxTreeIncrementalFilter`:
     * when the caller only needs a temporary detached subtree, `editor.detach`
     * avoids the parent/child relation check required by `removeChild`.
     *
     * @param topology Tree topology.
     * @param editor Active tree editor.
     * @param nodeId Dense internal node identifier.
     */
    void detachIfNonRoot(const MorphologicalTree& topology, editor_t& editor, NodeId nodeId) {
        if (topology.isRoot(nodeId)) {
            return;
        }
        editor.detach(nodeId);
    }

    /**
     * @brief Detaches a node from its current parent and optionally releases it.
     *
     * Non-release detaches use the lighter `TreeEditor::detach` path. Release
     * operations keep using `removeChild(..., true)` because they must return the
     * now-empty node slot to the free list of the tree.
     *
     * @param tree Tree topology.
     * @param editor Active tree editor.
     * @param nodeId Dense internal node identifier.
     * @param releaseNode Node identifier.
     */
    void disconnect(tree_t* tree, editor_t& editor, NodeId nodeId, bool releaseNode) {
        const MorphologicalTree& topology = topologyOf(tree);
        if (!releaseNode) {
            detachIfNonRoot(topology, editor, nodeId);
            return;
        }
        if (topology.isRoot(nodeId)) {
            return;
        }

        const NodeId parentId = TreeAccess::nodeParent(topology, nodeId);
        if (parentId == InvalidNode || parentId == nodeId) {
            return;
        }
        editor.removeChild(parentId, nodeId, releaseNode);
    }

    /**
     * @brief Moves every child and every direct proper part from `childId` to `parentId`.
     *
     * After this operation, `childId` is structurally empty but still attached
     * until the caller explicitly disconnects or releases it.
     *
     * @param tree Tree topology.
     * @param editor Active tree editor.
     * @param parentId Identifier of the parent node.
     * @param childId Identifier of the child node.
     */
    void mergedParentAndChildren(tree_t* tree, editor_t& editor, NodeId parentId, NodeId childId) {
        const MorphologicalTree& topology = topologyOf(tree);
        editor.moveChildren(parentId, childId);
        if (TreeAccess::properPartCardinality(topology, childId) > 0) {
            lastProperPartTransfers_.push_back({childId, parentId});
            if (lastSelectedSupportNode_ == childId) {
                lastSelectedSupportNode_ = parentId;
            }
        }
        editor.mergeProperParts(parentId, childId);
    }

    /**
     * @brief Moves the selected support `C` to the current union node.
     *
     * The support is represented by the direct proper parts of `leafId` in the
     * primal tree. Each pixel is looked up in the dual tree and moved from its
     * current smallest node to `unionNode`. Donor nodes that become empty are marked for
     * deferred absorption; if the emptied donor is `nodeCa`, the caller must handle
     * the special replacement/root cases after the sweep.
     *
     * @param dualTree Dual tree.
     * @param editor Active tree editor.
     * @param unionNode Node identifier.
     * @param support Direct proper-part support of the selected primal leaf.
     * @param nodeCa Unique dual node whose proper part contains the selected support.
     * @param nodeCaRemoved Whether the common-ancestor node was removed.
     * @param certifiedWholeSupportSmallestNode Optional certified smallest node of the complete
     *        support, enabling one bulk proper-part transfer.
     */
    void moveSelectedProperPartsToNode(tree_t* dualTree, editor_t& editor, NodeId unionNode, std::span<const PixelId> support, NodeId nodeCa,
                                       bool& nodeCaRemoved, NodeId certifiedWholeSupportSmallestNode) {
        const MorphologicalTree& dualTopology = topologyOf(dualTree);
        lastSelectedSupportNode_ = unionNode;
        if (certifiedWholeSupportSmallestNode != InvalidNode && certifiedWholeSupportSmallestNode != unionNode) {
            lastProperPartTransfers_.push_back({certifiedWholeSupportSmallestNode, unionNode});
            editor.mergeProperParts(unionNode, certifiedWholeSupportSmallestNode);
            if (TreeAccess::properPartCardinality(dualTopology, certifiedWholeSupportSmallestNode) == 0 &&
                !removedMarks_.isMarked(static_cast<std::size_t>(certifiedWholeSupportSmallestNode))) {
                removedMarks_.mark(static_cast<std::size_t>(certifiedWholeSupportSmallestNode));
                removedNodesPendingAbsorption_.push_back(certifiedWholeSupportSmallestNode);
                if (certifiedWholeSupportSmallestNode == nodeCa) {
                    nodeCaRemoved = true;
                }
            }
            return;
        }
        for (PixelId pixelId : support) {
            const NodeId smallestNodeId = TreeAccess::smallestNodeMap(dualTopology, pixelId);
            if (smallestNodeId == unionNode) {
                continue;
            }

            editor.movePixelToProperPart(unionNode, smallestNodeId, pixelId);
            if (TreeAccess::properPartCardinality(dualTopology, smallestNodeId) == 0 && !removedMarks_.isMarked(static_cast<std::size_t>(smallestNodeId))) {
                removedMarks_.mark(static_cast<std::size_t>(smallestNodeId));
                removedNodesPendingAbsorption_.push_back(smallestNodeId);
                if (smallestNodeId == nodeCa) {
                    nodeCaRemoved = true;
                }
            }
        }
    }

    /**
     * @brief Tests whether a marked empty node is still alive and safe to absorb.
     *
     * @param tree Tree topology.
     * @param nodeId Dense internal node identifier.
     * @return True if a marked empty node is still alive and safe to absorb; otherwise false.
     */
    bool canAbsorbRemovedNode(tree_t* tree, NodeId nodeId) const {
        return nodeId != InvalidNode && TreeAccess::isAlive(topologyOf(tree), nodeId) &&
               removedMarks_.isMarked(static_cast<std::size_t>(nodeId)) && TreeAccess::properPartCardinality(topologyOf(tree), nodeId) == 0;
    }

    /**
     * @brief Contracts an empty root by promoting its altitude-compatible child.
     *
     * If a removed node is currently the root, its children cannot simply be moved
     * to a parent. The child with the extremal altitude for the dual tree polarity
     * becomes the new root and all other children are reattached below it.
     *
     * @param tree Tree topology.
     * @param editor Active tree editor.
     * @param removedNodeId Node identifier.
     */
    void absorbRemovedRootNode(tree_t* tree, editor_t& editor, NodeId removedNodeId) {
        const bool isMaxtree = tree == maxtree_;
        const MorphologicalTree& topology = topologyOf(tree);
        const NodeId firstChild = TreeAccess::firstChild(topology, removedNodeId);
        if (firstChild == InvalidNode) {
            return;
        }

        NodeId newRoot = firstChild;
        for (NodeId childId = firstChild; childId != InvalidNode; childId = TreeAccess::nextSibling(topology, childId)) {
            if ((isMaxtree && nodeAltitude(tree, childId) < nodeAltitude(tree, newRoot)) ||
                (!isMaxtree && nodeAltitude(tree, childId) > nodeAltitude(tree, newRoot))) {
                newRoot = childId;
            }
        }

        for (NodeId childId = firstChild; childId != InvalidNode;) {
            const NodeId nextChildId = TreeAccess::nextSibling(topology, childId);
            if (childId != newRoot && TreeAccess::nodeParent(topology, childId) != newRoot) {
                detachIfNonRoot(topology, editor, childId);
                editor.attach(newRoot, childId);
            }
            childId = nextChildId;
        }
        editor.setRoot(newRoot);
        editor.releaseNode(removedNodeId);
    }

    /**
     * @brief Contracts one empty non-root node into its current parent.
     *
     * Children and direct proper parts of the node are first moved to the parent,
     * then the empty node is released. If the parent becomes empty as a result of
     * the merge, the parent is marked so the post-order absorption pass can keep
     * contracting upward.
     *
     * @param tree Tree topology.
     * @param editor Active tree editor.
     * @param removedNodeId Node identifier.
     * @return The surviving parent empty non-root node into its current parent.
     */
    NodeId absorbRemovedNonRootNode(tree_t* tree, editor_t& editor, NodeId removedNodeId) {
        const MorphologicalTree& topology = topologyOf(tree);
        const NodeId parentId = TreeAccess::nodeParent(topology, removedNodeId);
        if (parentId == removedNodeId) {
            return InvalidNode;
        }

        mergedParentAndChildren(tree, editor, parentId, removedNodeId);
        disconnect(tree, editor, removedNodeId, true);

        if (TreeAccess::properPartCardinality(topology, parentId) == 0) {
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
     * @param tree Tree topology.
     * @param editor Active tree editor.
     * @param removedNodeIds Node identifier.
     */
    void absorbRemovedNodes(tree_t* tree, editor_t& editor, const std::vector<NodeId>& removedNodeIds) {
        struct Frame {
            NodeId nodeId = InvalidNode;
            NodeId nextChildId = InvalidNode;
        };

        if (removedNodeIds.empty()) {
            return;
        }

        std::vector<Frame> stack;
        stack.reserve(std::max<std::size_t>(64, removedNodeIds.size()));
        const auto makeFrame = [tree](NodeId nodeId) {
            if (!TreeAccess::isAlive(topologyOf(tree), nodeId)) {
                return Frame{nodeId, InvalidNode};
            }
            return Frame{nodeId, TreeAccess::firstChild(topologyOf(tree), nodeId)};
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
                frame.nextChildId = TreeAccess::nextSibling(topologyOf(tree), childId);
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
     * @brief Contracts an established primal leaf through the existing edit primitives.
     * @param tree Primal tree containing an attached non-root leaf.
     * @param leafId Live leaf whose support is transferred to its parent.
     */
    static void pruneLeaf(tree_t& tree, NodeId leafId) {
        const NodeId parent = TreeAccess::nodeParent(tree.topology(), leafId);
        auto editor = mmcfilters::detail::beginEstablishedValuedEdit(tree);
        editor.mergeProperParts(parent, leafId);
        editor.removeChild(parent, leafId, true);
        auto proof = editor.proveIncremental();
        editor.commit(std::move(proof));
    }

  public:
    /**
     * @brief Builds the leaf adjuster over a fixed valuedTree min-tree / max-tree pair.
     *
     * The adjuster does not own the trees or the adjacency relation. Prune
     * methods always operate on this same pair; they do not rebind trees per call.
     * Both trees must represent the same image domain and must use node ids from
     * the internal dense node-id space of this repository.
     *
     * @param minTree Min-tree.
     * @param maxTree Max-tree.
     * @param graph Region or adjacency graph.
     */
    DualMinMaxTreeIncrementalFilterLeaf(tree_t* minTree, tree_t* maxTree, const RegularGridAdjacency2D& graph)
        : mintree_(minTree), maxtree_(maxTree), graph_(&graph),
          removedMarks_(static_cast<std::size_t>(std::max(minTree->topology().numInternalNodeSlots(), maxTree->topology().numInternalNodeSlots()))),
          climbedNodeMarks_(static_cast<std::size_t>(std::max(minTree->topology().numInternalNodeSlots(), maxTree->topology().numInternalNodeSlots()))) {}

    /**
     * @brief Returns dual-tree nodes that may affect the candidate agenda of the builder.
     *
     * SDRT construction uses this after each leaf prune to refresh only affected
     * dual candidates. The vector is cleared at the beginning of each update and
     * receives the union node created at every consumed merge level, matching the
     * MorphoTreeAdjust leaf-builder behavior.
     *
     * @return Dual-tree nodes that may affect the candidate agenda of the builder.
     */
    const std::vector<NodeId>& getLastCandidateNodes() const { return lastCandidateNodes_; }

    /**
     * @brief Exposes complete dual proper-part transfers from the last update.
     * @return Borrowed transfer records in execution order, valid until the next update.
     */
    [[nodiscard]] std::span<const ProperPartTransfer> getLastProperPartTransfers() const noexcept { return lastProperPartTransfers_; }

    /**
     * @brief Identifies the dual proper part containing the support after the last update.
     * @return Destination established by selected-pixel moves and subsequent whole-part merges.
     */
    [[nodiscard]] NodeId getLastSelectedSupportNode() const noexcept { return lastSelectedSupportNode_; }

  private:
    /**
     * @brief Builds merge buckets and frontiers for one removed leaf.
     *
     * Starting from the supplied primal leaf support `C`, the method scans its
     * graph neighbors and maps them to their smallest nodes in `dualTree`.
     * Neighbors mapped to `nodeCa` are skipped: the complete support lies in that
     * proper part, and the merge sweep excludes the altitude of `nodeCa`.
     * This avoids a pixel-mark buffer and a separate pass over `C`. Other dual
     * nodes seed ancestor paths in the altitude interval governed by `nodeCa` and `b`.
     * Each accepted node is marked before insertion. The dual tree is unchanged
     * during collection, so one mark set suffices for both seeds and ancestors.
     *
     * Nodes whose level is inside the merge side of the interval are inserted in
     * `mergeNodesByLevel_`. Boundary nodes outside the merge side are recorded as
     * frontiers and later attached to the union node created at level `b`.
     *
     * @param dualTree Dual tree.
     * @param support Direct proper-part support of the selected primal leaf.
     * @param nodeCa Unique dual node whose proper part contains the selected support.
     * @param b Primal parent altitude reached by this leveling.
     * @param isMaxtree Flag controlling is maxtree.
     */
    void buildMergedAndNestedCollections(tree_t* dualTree, std::span<const PixelId> support, NodeId nodeCa, altitude_t b, bool isMaxtree) {
        const MorphologicalTree& dualTopology = topologyOf(dualTree);
        const altitude_t altitudeCa = nodeAltitude(dualTree, nodeCa);

        const auto addAdjacentNode = [&](NodeId nodeQ) {
            NodeId n = nodeQ;
            while (!climbedNodeMarks_.isMarked(static_cast<std::size_t>(n))) {
                const altitude_t levelCurrent = nodeAltitude(dualTree, n);
                if (!((isMaxtree && levelCurrent >= altitudeCa) || (!isMaxtree && levelCurrent <= altitudeCa))) {
                    break;
                }

                climbedNodeMarks_.mark(static_cast<std::size_t>(n));

                if ((isMaxtree && levelCurrent <= b) || (!isMaxtree && levelCurrent >= b)) {
                    mergeNodesByLevel_.addMergeNode(*dualTree, n);
                } else {
                    NodeId parentId = TreeAccess::nodeParent(dualTopology, n);
                    if (parentId == n) {
                        parentId = InvalidNode;
                    }
                    if (!(parentId != InvalidNode &&
                          ((isMaxtree && nodeAltitude(dualTree, parentId) > b) || (!isMaxtree && nodeAltitude(dualTree, parentId) < b)))) {
                        mergeNodesByLevel_.addFrontierNodeAboveB(n);
                    }
                }

                const NodeId parentId = TreeAccess::nodeParent(dualTopology, n);
                if (parentId == n) {
                    break;
                }
                n = parentId;
            }
        };

        for (PixelId pixel : support) {
            for (PixelId neighbor : mmcfilters::detail::CommittedGridAccess::neighbors(*graph_, pixel)) {
                const NodeId nodeQ = TreeAccess::smallestNodeMap(dualTopology, neighbor);
                if (nodeQ != nodeCa) {
                    addAdjacentNode(nodeQ);
                }
            }
        }
    }

    /**
     * @brief Updates the dual tree for the selected primal leaf leveling.
     *
     * This is the algorithmic core of leaf mode. The method:
     *
     * - identifies the parent level `b` of the primal leaf;
     * - uses the supplied dual node `nodeCa` whose proper part contains the leaf support;
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
     * @param dualTree Dual tree.
     * @param leafId Live attached non-root leaf with a nonempty proper part.
     * @param support Nonempty proper-part support supplied by the builder.
     * @param nodeCa Unique live dual node whose proper part contains the supplied support.
     * @param certifiedWholeSupportSmallestNode Optional certified smallest node of the complete
     *        support, enabling one bulk proper-part transfer.
     */
    void updateTree(tree_t* dualTree, NodeId leafId, std::span<const PixelId> support, NodeId nodeCa,
                    NodeId certifiedWholeSupportSmallestNode) {
        const bool isMaxtree = dualTree == maxtree_;
        initializeUpdateStepState(isMaxtree);
        tree_t* primalTree = getPrimalTree(isMaxtree);
        const MorphologicalTree& dualTopology = topologyOf(dualTree);
        const MorphologicalTree& primalTopology = topologyOf(primalTree);

        const NodeId leafParent = TreeAccess::nodeParent(primalTopology, leafId);
        const altitude_t b = nodeAltitude(primalTree, leafParent);
        const altitude_t altitudeCa = nodeAltitude(dualTree, nodeCa);
        lastSelectedSupportNode_ = nodeCa;
        buildMergedAndNestedCollections(dualTree, support, nodeCa, b, isMaxtree);

        editor_t editor = mmcfilters::detail::beginEstablishedValuedEdit(*dualTree);
        altitude_t mergeLevel = mergeNodesByLevel_.firstMergeLevel();
        NodeId unionNode = InvalidNode;
        NodeId previousUnionNode = InvalidNode;
        bool nodeCaRemoved = false;

        while (mergeNodesByLevel_.hasMergeLevel() && ((isMaxtree && mergeLevel > altitudeCa) || (!isMaxtree && mergeLevel < altitudeCa))) {
            auto& mergeNodesAtLevel = mergeNodesByLevel_.getMergedNodes(mergeLevel);
            unionNode = InvalidNode;

            for (auto nodeId : mergeNodesAtLevel) {
                if (!TreeAccess::isAlive(dualTopology, nodeId)) {
                    continue;
                }
                if (unionNode == InvalidNode || nodeId < unionNode) {
                    unionNode = nodeId;
                }
            }

            for (auto nodeId : mergeNodesAtLevel) {
                if (!TreeAccess::isAlive(dualTopology, nodeId)) {
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
                moveSelectedProperPartsToNode(dualTree, editor, unionNode, support, nodeCa, nodeCaRemoved, certifiedWholeSupportSmallestNode);

                for (auto nodeId : mergeNodesByLevel_.getFrontierNodesAboveB()) {
                    if (!TreeAccess::isAlive(dualTopology, nodeId)) {
                        continue;
                    }
                    disconnect(dualTree, editor, nodeId, false);
                    editor.attach(unionNode, nodeId);
                }
            }

            if (previousUnionNode != InvalidNode && TreeAccess::isAlive(dualTopology, previousUnionNode) && TreeAccess::nodeParent(dualTopology, previousUnionNode) != unionNode) {
                if (!dualTopology.isRoot(previousUnionNode)) {
                    disconnect(dualTree, editor, previousUnionNode, false);
                }
                editor.attach(unionNode, previousUnionNode);
            }

            previousUnionNode = unionNode;
            mergeLevel = mergeNodesByLevel_.nextMergeLevel();
        }

        const NodeId finalUnionNode = previousUnionNode;
        if (finalUnionNode == InvalidNode || !TreeAccess::isAlive(dualTopology, finalUnionNode)) {
            auto proof = editor.proveIncremental();
            editor.commit(std::move(proof));
            return;
        }

        if (nodeCaRemoved) {
            if (!dualTopology.isRoot(nodeCa)) {
                const NodeId parentIdNodeCa = TreeAccess::nodeParent(dualTopology, nodeCa);
                if (!dualTopology.isRoot(finalUnionNode)) {
                    disconnect(dualTree, editor, finalUnionNode, false);
                }
                editor.attach(parentIdNodeCa, finalUnionNode);

                for (NodeId n = TreeAccess::firstChild(dualTopology, nodeCa); n != InvalidNode;) {
                    const NodeId next = TreeAccess::nextSibling(dualTopology, n);
                    if (n != finalUnionNode && TreeAccess::nodeParent(dualTopology, n) != finalUnionNode) {
                        detachIfNonRoot(dualTopology, editor, n);
                        editor.attach(finalUnionNode, n);
                    }
                    n = next;
                }

                disconnect(dualTree, editor, nodeCa, true);
            } else {
                NodeId newRoot = finalUnionNode;
                for (NodeId n = TreeAccess::firstChild(dualTopology, nodeCa); n != InvalidNode; n = TreeAccess::nextSibling(dualTopology, n)) {
                    if ((isMaxtree && nodeAltitude(dualTree, n) < nodeAltitude(dualTree, newRoot)) ||
                        (!isMaxtree && nodeAltitude(dualTree, n) > nodeAltitude(dualTree, newRoot))) {
                        newRoot = n;
                    }
                }

                if (newRoot != finalUnionNode) {
                    if (!dualTopology.isRoot(finalUnionNode)) {
                        disconnect(dualTree, editor, finalUnionNode, false);
                    }
                    editor.attach(newRoot, finalUnionNode);
                }

                for (NodeId n = TreeAccess::firstChild(dualTopology, nodeCa); n != InvalidNode;) {
                    const NodeId next = TreeAccess::nextSibling(dualTopology, n);
                    if (n != newRoot && TreeAccess::nodeParent(dualTopology, n) != newRoot) {
                        detachIfNonRoot(dualTopology, editor, n);
                        editor.attach(newRoot, n);
                    }
                    n = next;
                }
                editor.setRoot(newRoot);
                editor.releaseNode(nodeCa);
            }
        } else if (finalUnionNode != nodeCa) {
            if (!dualTopology.isRoot(finalUnionNode)) {
                disconnect(dualTree, editor, finalUnionNode, false);
            }
            editor.attach(nodeCa, finalUnionNode);
        }

        absorbRemovedNodes(dualTree, editor, removedNodesPendingAbsorption_);
        auto proof = editor.proveIncremental();
        editor.commit(std::move(proof));
    }

  public:
    /**
     * @brief Prunes one certified max-tree leaf without rediscovering its workset.
     *
     * The support and certificates describe the current leaf and min-tree.
     * They are consumed synchronously before the primal leaf is pruned.
     *
     * @param leafId Live attached non-root leaf with a nonempty proper part.
     * @param certifiedSupport Nonempty selected leaf support supplied by the builder.
     * @param certifiedNodeCa Unique live dual node whose proper part contains the support.
     * @param certifiedWholeSupportSmallestNode Dual node whose entire proper part is selected, or InvalidNode.
     */
    void pruneMaxLeafAndUpdateMinTree(NodeId leafId, std::span<const PixelId> certifiedSupport, NodeId certifiedNodeCa,
                                     NodeId certifiedWholeSupportSmallestNode) {
        updateTree(mintree_, leafId, certifiedSupport, certifiedNodeCa, certifiedWholeSupportSmallestNode);
        pruneLeaf(*maxtree_, leafId);
    }

    /**
     * @brief Prunes one certified min-tree leaf without rediscovering its workset.
     *
     * @param leafId Live attached non-root leaf with a nonempty proper part.
     * @param certifiedSupport Nonempty selected leaf support supplied by the builder.
     * @param certifiedNodeCa Unique live dual node whose proper part contains the support.
     * @param certifiedWholeSupportSmallestNode Dual node whose entire proper part is selected, or InvalidNode.
     */
    void pruneMinLeafAndUpdateMaxTree(NodeId leafId, std::span<const PixelId> certifiedSupport, NodeId certifiedNodeCa,
                                     NodeId certifiedWholeSupportSmallestNode) {
        updateTree(maxtree_, leafId, certifiedSupport, certifiedNodeCa, certifiedWholeSupportSmallestNode);
        pruneLeaf(*mintree_, leafId);
    }
};

} // namespace mmcfilters::detail
