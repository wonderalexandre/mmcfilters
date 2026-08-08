#pragma once

#include "DynamicTreeAttributeComputer.hpp"
#include "../WeightedMorphologicalTree.hpp"
#include "../../utils/GenerationStampSet.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <type_traits>
#include <utility>
#include <vector>

#ifndef MMCFILTERS_COMPONENT_TREE_ADJUSTMENT_DENSE_MAX_BITS
#define MMCFILTERS_COMPONENT_TREE_ADJUSTMENT_DENSE_MAX_BITS 8
#endif

namespace mmcfilters::adjust {

/**
 * @brief Incremental updater for a paired min-tree / max-tree state.
 *
 * This class ports the Higra dual component-tree adjustment algorithm to this
 * repository's `WeightedMorphologicalTree<T>` representation. It implements the
 * update-rather-than-rebuild strategy used for efficient connected alternating
 * sequential filters: when a rooted subtree is removed from one component tree,
 * the dual tree is updated in place so both trees remain consistent with the
 * same filtered image.
 *
 * The algorithmic structure is intentionally kept the same:
 *
 * - a rooted subtree is selected for pruning in the primal tree;
 * - its full support defines the proper-part set `C`;
 * - valid graph-adjacent components around `C` are collected in the dual tree;
 * - each relevant ancestor chain is climbed once;
 * - collected nodes are bucketed by altitude and swept in tree order;
 * - local unions, frontier reattachments, and empty-node contractions update
 *   the dual tree without rebuilding it globally.
 *
 * Public usage model:
 *
 * - construct the helper from an externally owned min-tree / max-tree pair and
 *   their shared adjacency relation;
 * - optionally register incremental attribute computers and attribute buffers
 *   to refresh after local edits;
 * - call either `pruneMaxTreeAndUpdateMinTree(...)` or
 *   `pruneMinTreeAndUpdateMaxTree(...)`.
 *
 * State managed by this class:
 *
 * - two mutable weighted component trees (`mintree_`, `maxtree_`);
 * - one adjacency relation shared by both trees;
 * - optional incremental attribute computers and two external attribute
 *   buffers;
 * - transient merge buckets, frontier roots, proper-part sets, and mark sets
 *   reused across update steps.
 *
 * Ownership and lifetime:
 *
 * - this class does not own the trees, adjacency relation, attribute computers,
 *   or attribute buffers;
 * - callers must keep these objects alive for the whole lifetime of the filter;
 * - external buffers must be indexed on this project's internal node-id space.
 *
 * Representation notes:
 *
 * - node ids are this project's dense internal node ids, not Higra global ids;
 * - altitudes are read from `WeightedMorphologicalTree<T>::getAltitude`;
 * - this is a mutable weighted-tree boundary, not a `WeightedTreeView`
 *   boundary: topology, proper-part ownership, owned altitude, and
 *   optional dynamic attributes are updated in lockstep;
 * - the `altitude_t` template parameter is the owned tree altitude type and
 *   also controls merge-bucket keys/backends;
 * - low-level edits go through `WeightedTreeEditor<T>` and publish with the
 *   generic move-only proof protocol; the update establishes its invariants
 *   during the existing algorithmic passes, while assertion-enabled builds
 *   retain the complete validation oracle;
 * - optional attribute buffers are kept in sync through
 *   `DynamicTreeAttributeComputer<std::uint8_t>`.
 *
 * Internal notation follows the original adjustment routine:
 *
 * - `nodeCa`: implementation-side representative of the extremal node in the
 *   dual tree that contains the set `C` at level `a`;
 * - `b`: altitude of the parent of the subtree removed in the primal tree;
 * - `frontierNodesAboveB`: roots that leave the merge interval and are
 *   reattached when the sweep reaches the boundary level;
 * - `properPartSetC_`: pixels directly removed or relocated in the dual tree.
 *
 * Backend selection for `mergeNodesByLevel`:
 *
 * - dense array-based buckets for small integral altitude domains, up to
 *   `MMCFILTERS_COMPONENT_TREE_ADJUSTMENT_DENSE_MAX_BITS` bits;
 * - sparse ordered-map buckets for larger integral domains and floating-point
 *   altitudes.
 *
 * The expected cost remains local to the affected region. The dominant terms
 * are the same as in the Higra implementation: collecting `C`, visiting
 * adjacent seeds, climbing relevant ancestor chains once, sweeping the active
 * altitude buckets, and contracting the nodes actually emptied by the step.
 * The class never rebuilds the dual tree globally as part of an adjustment.
 */
template <AltitudeValue altitude_t> class DualMinMaxTreeIncrementalFilter {
  private:
    /** @brief Defines the `tree_t` alias used by the component. */
    using tree_t = WeightedMorphologicalTree<altitude_t>;
    /** @brief Defines the `editor_t` alias used by the component. */
    using editor_t = WeightedTreeEditor<altitude_t>;
    /** @brief Defines the `attribute_computer_t` alias used by the component. */
    using attribute_computer_t = DynamicTreeAttributeComputer<altitude_t>;

    /**
     * @brief Returns `true` iff the altitude type should use dense buckets.
     * @details Dense buckets are enabled only for small integral altitude
     * domains. Larger integral domains and floating-point types use the sparse
     * backend to avoid allocating a bucket for every possible value.
     *
     * @return True iff the altitude type should use dense buckets.
     */
    static constexpr bool usesDenseLevels() {
        if constexpr (std::is_integral_v<altitude_t>) {
            using unsigned_altitude_t = std::make_unsigned_t<altitude_t>;
            return std::numeric_limits<unsigned_altitude_t>::digits <= MMCFILTERS_COMPONENT_TREE_ADJUSTMENT_DENSE_MAX_BITS;
        } else {
            return false;
        }
    }

    /**
     * @brief True iff `mergeNodesByLevel` uses the dense bucket backend.
     */
    static constexpr bool use_dense_levels = usesDenseLevels();

    /**
     * @brief Per-step collection of merge buckets and frontier roots.
     *
     * Nodes adjacent to `C` are grouped by altitude for the subsequent sweep.
     * Small integral altitude domains use a dense bucket array; larger integral
     * domains and floating-point instantiations use a sparse ordered map. This
     * matches the Higra backend choice and avoids imposing an 8-bit altitude
     * domain on this project.
     */
    class MergedNodesCollection {
      private:
        template <bool dense, typename altitude_type> struct StorageSelector;

        /**
         * @brief Dense bucket storage for small integral altitude domains.
         */
        template <typename altitude_type> struct StorageSelector<true, altitude_type> {
            /** @brief Stores the domain size. */
            static constexpr std::size_t domain_size = static_cast<std::size_t>(static_cast<long long>(std::numeric_limits<altitude_type>::max()) -
                                                                                static_cast<long long>(std::numeric_limits<altitude_type>::lowest()) + 1);
            /** @brief Defines the `type` alias used by the component. */
            using type = std::array<std::vector<NodeId>, domain_size>;
        };

        /**
         * @brief Sparse ordered storage for large or non-integral altitude domains.
         */
        template <typename altitude_type> struct StorageSelector<false, altitude_type> {
            /** @brief Defines the `type` alias used by the component. */
            using type = std::map<altitude_type, std::vector<NodeId>>;
        };

        /** @brief Defines the `storage_t` alias used by the component. */
        using storage_t = typename StorageSelector<use_dense_levels, altitude_t>::type;

        /** @brief Stores the merge nodes by level storage. */
        storage_t mergeNodesByLevelStorage_;
        /** @brief Stores the merge levels. */
        std::vector<altitude_t> mergeLevels_;
        /** @brief Stores the frontier nodes above b. */
        std::vector<NodeId> frontierNodesAboveB_;
        /** @brief Stores the collected node marks. */
        GenerationStampSet collectedNodeMarks_;
        /** @brief Stores the merge bucket node marks. */
        GenerationStampSet mergeBucketNodeMarks_;
        /** @brief Stores the adjacent seed marks. */
        GenerationStampSet adjacentSeedMarks_;
        /** @brief Stores the max bucket size. */
        std::size_t maxBucketSize_ = 0;
        /** @brief Stores the current merge level index. */
        int currentMergeLevelIndex_ = 0;
        /** @brief Indicates whether this side represents the max-tree. */
        bool isMaxtree_ = false;

      public:
        /**
         * @brief Maps an altitude value to an order-preserving dense bucket index.
         * @details Signed integral domains are shifted by their lowest value so
         * negative altitudes still map to non-negative array indices.
         *
         * @param level Altitude level used by the operation.
         * @return The mapped altitude value to an order-preserving dense bucket index.
         */
        static std::size_t denseBucketIndex(altitude_t level) {
            if constexpr (std::is_signed_v<altitude_t>) {
                return static_cast<std::size_t>(static_cast<long long>(level) - static_cast<long long>(std::numeric_limits<altitude_t>::lowest()));
            } else {
                return static_cast<std::size_t>(level);
            }
        }

        /**
         * @brief Converts a dense bucket index back to its altitude value.
         *
         * @param index Zero-based index used by the operation.
         * @return The converted dense bucket index back to its altitude value.
         */
        static altitude_t levelFromDenseBucketIndex(std::size_t index) {
            if constexpr (std::is_signed_v<altitude_t>) {
                return static_cast<altitude_t>(static_cast<long long>(index) + static_cast<long long>(std::numeric_limits<altitude_t>::lowest()));
            } else {
                return static_cast<altitude_t>(index);
            }
        }

        /**
         * @brief Creates a merged-node collection sized for `maxNodes`.
         * @param maxNodes Size of the internal node-id space used by mark sets.
         */
        explicit MergedNodesCollection(int maxNodes = 0)
            : collectedNodeMarks_(static_cast<size_t>(std::max(maxNodes, 0))), mergeBucketNodeMarks_(static_cast<size_t>(std::max(maxNodes, 0))),
              adjacentSeedMarks_(static_cast<size_t>(std::max(maxNodes, 0))) {}

        /**
         * @brief Resets the transient state of one adjustment step.
         * @param isMaxtree `true` for max-tree updates, `false` for min-tree updates.
         */
        void resetCollection(bool isMaxtree) {
            isMaxtree_ = isMaxtree;
            for (altitude_t level : mergeLevels_) {
                if constexpr (use_dense_levels) {
                    // Dense backend reuses one bucket per altitude value and
                    // clears only levels that were active in the previous step.
                    mergeNodesByLevelStorage_[denseBucketIndex(level)].clear();
                } else {
                    // Sparse backend stores only touched levels, so only those
                    // buckets need to be cleared here.
                    auto it = mergeNodesByLevelStorage_.find(level);
                    if (it != mergeNodesByLevelStorage_.end()) {
                        it->second.clear();
                    }
                }
            }
            mergeLevels_.clear();
            frontierNodesAboveB_.clear();
            collectedNodeMarks_.resetAll();
            mergeBucketNodeMarks_.resetAll();
            adjacentSeedMarks_.resetAll();
            maxBucketSize_ = 0;
            currentMergeLevelIndex_ = 0;
        }

        /**
         * @brief Returns the bucket associated with one altitude value.
         *
         * @param level Altitude level used by the operation.
         * @return The bucket associated with one altitude value.
         */
        std::vector<NodeId>& getMergedNodes(const altitude_t& level) {
            if constexpr (use_dense_levels) {
                return mergeNodesByLevelStorage_[denseBucketIndex(level)];
            } else {
                return mergeNodesByLevelStorage_[level];
            }
        }

        /**
         * @brief Returns the frontier roots collected above `b`.
         * @details Each stored node is the root of the first branch that leaves
         * the merge interval while climbing from a valid adjacent seed.
         *
         * @return The frontier roots collected above b.
         */
        std::vector<NodeId>& getFrontierNodesAboveB() { return frontierNodesAboveB_; }

        /**
         * @brief Returns the largest bucket observed in the current build.
         * @details This value is used to reserve the sweep worklist.
         *
         * @return The largest bucket observed in the current build.
         */
        std::size_t getMaxBucketSize() const { return maxBucketSize_; }

        /**
         * @brief Marks one adjacent seed at most once in the current step.
         *
         * @param nodeId Identifier of the node used by the operation.
         * @return `true` if the seed was accepted now; `false` if it had already been seen.
         *
         */
        bool markAdjacentSeed(NodeId nodeId) {
            if (adjacentSeedMarks_.isMarked(static_cast<size_t>(nodeId))) {
                return false;
            }
            adjacentSeedMarks_.mark(static_cast<size_t>(nodeId));
            return true;
        }

        /**
         * @brief Registers one frontier root above `b` only once.
         *
         * @param nodeId Identifier of the node used by the operation.
         */
        void addFrontierNodeAboveB(NodeId nodeId) {
            if (!collectedNodeMarks_.isMarked(static_cast<size_t>(nodeId))) {
                frontierNodesAboveB_.push_back(nodeId);
                collectedNodeMarks_.mark(static_cast<size_t>(nodeId));
            }
        }

        /**
         * @brief Inserts one node into the bucket of its own altitude.
         * @details Each visited node is registered only in the bucket of its own
         * level, without materializing the full ancestor path in advance.
         *
         * @param tree Tree topology used by the operation.
         * @param nodeId Identifier of the node used by the operation.
         */
        void addMergeNode(const tree_t& tree, NodeId nodeId) {
            if (collectedNodeMarks_.isMarked(static_cast<size_t>(nodeId))) {
                return;
            }
            auto& bucket = getMergedNodes(tree.getAltitude(nodeId));
            bucket.push_back(nodeId);
            maxBucketSize_ = std::max(maxBucketSize_, bucket.size());
            collectedNodeMarks_.mark(static_cast<size_t>(nodeId));
            mergeBucketNodeMarks_.mark(static_cast<size_t>(nodeId));
        }

        /**
         * @brief Tests whether a node was inserted in a merge bucket.
         *
         * @param nodeId Identifier of the node used by the operation.
         * @return True if a node was inserted in a merge bucket; otherwise false.
         */
        bool isMergeNode(NodeId nodeId) const { return mergeBucketNodeMarks_.isMarked(static_cast<size_t>(nodeId)); }

        /**
         * @brief Builds the ordered list of active levels and returns the first one.
         * @details The order matches the sweep direction of the current tree:
         * descending levels for max-tree updates and ascending levels for
         * min-tree updates.
         *
         * @return The resulting ordered list of active levels and returns the first one.
         */
        altitude_t firstMergeLevel() {
            mergeLevels_.clear();
            if constexpr (use_dense_levels) {
                // Dense backend scans the discrete altitude domain and keeps
                // only non-empty levels of the current step.
                for (std::size_t i = 0; i < mergeNodesByLevelStorage_.size(); ++i) {
                    if (!mergeNodesByLevelStorage_[i].empty()) {
                        mergeLevels_.push_back(levelFromDenseBucketIndex(i));
                    }
                }
            } else {
                // Sparse backend iterates only over levels that were instantiated.
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
            return mergeLevels_[static_cast<size_t>(currentMergeLevelIndex_)];
        }

        /**
         * @brief Returns `true` while the current level iterator is valid.
         *
         * @return True while the current level iterator is valid.
         */
        bool hasMergeLevel() const {
            return !mergeLevels_.empty() && currentMergeLevelIndex_ >= 0 && currentMergeLevelIndex_ < static_cast<int>(mergeLevels_.size());
        }

        /**
         * @brief Advances the ordered level iterator and returns the next level.
         *
         * @return Next level in the ordered sweep.
         */
        altitude_t nextMergeLevel() {
            currentMergeLevelIndex_ = isMaxtree_ ? currentMergeLevelIndex_ - 1 : currentMergeLevelIndex_ + 1;
            if (!hasMergeLevel()) {
                return altitude_t{};
            }
            return mergeLevels_[static_cast<size_t>(currentMergeLevelIndex_)];
        }
    };

    // Dynamic primal/dual trees and shared domain adjacency.
    /** @brief Stores the mintree. */
    tree_t* mintree_ = nullptr;
    /** @brief Stores the maxtree. */
    tree_t* maxtree_ = nullptr;
    /** @brief Stores the graph. */
    const RegularGridAdjacency2D* graph_ = nullptr;

    // Transient state of the current adjustment step.
    /** @brief Stores the merge nodes by level. */
    MergedNodesCollection mergeNodesByLevel_;
    /** @brief Stores the removed marks. */
    GenerationStampSet removedMarks_;
    /** @brief Stores the pixels in c marks. */
    GenerationStampSet pixelsInCMarks_;
    /** @brief Stores the climbed node marks. */
    GenerationStampSet climbedNodeMarks_;
    /** @brief Stores the attribute update marks. */
    GenerationStampSet attributeUpdateMarks_;
    /** @brief Stores the proper part set c. */
    std::vector<NodeId> properPartSetC_;
    /** @brief Stores the nodes pending removal. */
    std::vector<NodeId> nodesPendingRemoval_;
    /** @brief Stores the removed nodes pending absorption. */
    std::vector<NodeId> removedNodesPendingAbsorption_;
    /** @brief Stores the altitude ca. */
    altitude_t altitudeCa_ = altitude_t{};

    // Incremental attribute computation and external attribute buffers.
    /** @brief Stores the attribute computer min. */
    const attribute_computer_t* attributeComputerMin_ = nullptr;
    /** @brief Stores the attribute computer max. */
    const attribute_computer_t* attributeComputerMax_ = nullptr;
    /** @brief Stores the attribute buffer min. */
    std::vector<double>* attributeBufferMin_ = nullptr;
    /** @brief Stores the attribute buffer max. */
    std::vector<double>* attributeBufferMax_ = nullptr;

    /**
     * @brief Returns the topology view of one weighted tree.
     *
     * @param tree Tree topology used by the operation.
     * @return The topology view of one weighted tree.
     */
    static const MorphologicalTree& topologyOf(const tree_t* tree) {
        assert(tree != nullptr);
        return tree->topology();
    }

    /**
     * @brief Returns the attribute buffer associated with one tree.
     *
     * @param tree Tree topology used by the operation.
     * @return `nullptr` when no incremental attribute computer is configured.
     *
     */
    typename attribute_computer_t::buffer_type* getAttributeBuffer(tree_t* tree) const {
        const auto* computer = tree == maxtree_ ? attributeComputerMax_ : attributeComputerMin_;
        if (computer == nullptr) {
            return nullptr;
        }
        return tree == maxtree_ ? attributeBufferMax_ : attributeBufferMin_;
    }

    /**
     * @brief Returns the attribute computer associated with one tree.
     *
     * @param tree Tree topology used by the operation.
     * @return The attribute computer associated with one tree.
     */
    const attribute_computer_t* getAttributeComputer(tree_t* tree) const {
        if (tree == nullptr) {
            return nullptr;
        }
        return tree == maxtree_ ? attributeComputerMax_ : attributeComputerMin_;
    }

    /**
     * @brief Forwards a node-removal event to the incremental attribute computer of `tree`.
     * @details The notification is skipped when no computer is configured or
     * when the removed node id is invalid.
     *
     * @param tree Tree topology used by the operation.
     * @param nodeId Identifier of the node used by the operation.
     */
    void notifyNodeRemoved(tree_t* tree, NodeId nodeId) const {
        const auto* computer = getAttributeComputer(tree);
        if (computer != nullptr && tree != nullptr && nodeId != InvalidNode) {
            computer->onNodeRemoved(nodeId, *tree);
        }
    }

    /**
     * @brief Forwards a bulk proper-part transfer to the incremental attribute computer of `tree`.
     *
     * @param tree Tree topology used by the operation.
     * @param targetNodeId Node identifier represented by `targetNodeId`.
     * @param sourceNodeId Node identifier represented by `sourceNodeId`.
     */
    void notifyMoveProperParts(tree_t* tree, NodeId targetNodeId, NodeId sourceNodeId) const {
        const auto* computer = getAttributeComputer(tree);
        if (computer != nullptr && tree != nullptr) {
            computer->onMoveProperParts(targetNodeId, sourceNodeId, *tree);
        }
    }

    /**
     * @brief Forwards a single proper-part transfer to the incremental attribute computer of `tree`.
     *
     * @param tree Tree topology used by the operation.
     * @param targetNodeId Node identifier represented by `targetNodeId`.
     * @param sourceNodeId Node identifier represented by `sourceNodeId`.
     * @param pixelId Pixel identifier used by the operation.
     */
    void notifyMoveProperPart(tree_t* tree, NodeId targetNodeId, NodeId sourceNodeId, NodeId pixelId) const {
        const auto* computer = getAttributeComputer(tree);
        if (computer != nullptr && tree != nullptr) {
            computer->onMoveProperPart(targetNodeId, sourceNodeId, pixelId, *tree);
        }
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
     * @brief Refreshes one marked node attribute after local topology stabilization.
     *
     * Nodes outside the active altitude interval are unmarked without
     * recomputation. This mirrors the original implementation: attributes are
     * recomputed only on nodes whose local child/proper-part state can have
     * changed in the active sweep.
     *
     * @param tree Tree topology used by the operation.
     * @param nodeId Identifier of the node used by the operation.
     */
    void computeAttributeOnTreeNode(tree_t* tree, NodeId nodeId) {
        const auto* computer = getAttributeComputer(tree);
        if (computer == nullptr || tree == nullptr || nodeId == InvalidNode || !topologyOf(tree).isNode(nodeId) || !topologyOf(tree).isAlive(nodeId)) {
            return;
        }

        auto* buffer = getAttributeBuffer(tree);
        assert(buffer != nullptr);
        if (!attributeUpdateMarks_.isMarked(static_cast<size_t>(nodeId))) {
            return;
        }

        const altitude_t level = nodeAltitude(tree, nodeId);
        if ((tree == maxtree_ && level < altitudeCa_) || (tree == mintree_ && level > altitudeCa_)) {
            attributeUpdateMarks_.unmark(static_cast<size_t>(nodeId));
            return;
        }

        computer->computeAttributeOnNode(*tree, nodeId, *buffer);
        attributeUpdateMarks_.unmark(static_cast<size_t>(nodeId));
    }

    /**
     * @brief Clears the per-step attribute refresh marks and active reference level.
     */
    void resetAttributeUpdateMarks() {
        attributeUpdateMarks_.resetAll();
        altitudeCa_ = altitude_t{};
    }

    /**
     * @brief Marks one alive node for later attribute recomputation.
     *
     * @param tree Tree topology used by the operation.
     * @param nodeId Identifier of the node used by the operation.
     */
    void markAttributeUpdate(tree_t* tree, NodeId nodeId) {
        if (getAttributeComputer(tree) == nullptr || tree == nullptr || nodeId == InvalidNode) {
            return;
        }
        if (!topologyOf(tree).isNode(nodeId) || !topologyOf(tree).isAlive(nodeId)) {
            return;
        }
        attributeUpdateMarks_.mark(static_cast<size_t>(nodeId));
    }

    /**
     * @brief Detaches a node from its parent, optionally releasing it.
     * @details Delegates to `WeightedTreeEditor<std::uint8_t>::removeChild`, preserving the
     * root case and ignoring nodes with invalid parents.
     *
     * @param tree Tree topology used by the operation.
     * @param editor Active tree editor used by the operation.
     * @param nodeId Identifier of the node used by the operation.
     * @param releaseNode Node identifier represented by `releaseNode`.
     */
    void disconnect(tree_t* tree, editor_t& editor, NodeId nodeId, bool releaseNode) {
        assert(tree != nullptr);
        const MorphologicalTree& topology = topologyOf(tree);
        if (topology.isRoot(nodeId)) {
            return;
        }
        const NodeId parentId = topology.getNodeParent(nodeId);
        if (parentId == InvalidNode || parentId == nodeId) {
            return;
        }
        editor.removeChild(parentId, nodeId, releaseNode);
        if (releaseNode) {
            notifyNodeRemoved(tree, nodeId);
        }
    }

    /**
     * @brief Moves the proper-part set `C` onto the current union node.
     *
     * Donor nodes that become empty are not removed immediately. They are marked
     * and later contracted in post-order after the sweep has finished reconnecting
     * the local hierarchy.
     *
     * @param dualTree Dual tree used by the operation.
     * @param editor Active tree editor used by the operation.
     * @param unionNode Node identifier represented by `unionNode`.
     * @param properPartSetC Proper-part data represented by `properPartSetC`.
     */
    void moveSelectedProperPartsToNode(tree_t* dualTree, editor_t& editor, NodeId unionNode, const std::vector<NodeId>& properPartSetC) {
        for (NodeId pixelId : properPartSetC) {
            const NodeId ownerId = topologyOf(dualTree).getProperPartOwner(pixelId);
            if (ownerId == InvalidNode || ownerId == unionNode) {
                continue;
            }

            notifyMoveProperPart(dualTree, unionNode, ownerId, pixelId);
            editor.moveProperPart(unionNode, ownerId, pixelId);

            if (topologyOf(dualTree).isAlive(ownerId) && topologyOf(dualTree).getNumProperParts(ownerId) == 0 && ownerId != unionNode &&
                !removedMarks_.isMarked(static_cast<size_t>(ownerId))) {
                removedMarks_.mark(static_cast<size_t>(ownerId));
                removedNodesPendingAbsorption_.push_back(ownerId);
            }
        }
    }

    /**
     * @brief Tests whether a marked empty node is still absorbable.
     *
     * @param dualTree Dual tree used by the operation.
     * @param nodeId Identifier of the node used by the operation.
     * @return True if a marked empty node is still absorbable; otherwise false.
     */
    bool canAbsorbRemovedNode(tree_t* dualTree, NodeId nodeId) const {
        return dualTree != nullptr && nodeId != InvalidNode && topologyOf(dualTree).isNode(nodeId) && topologyOf(dualTree).isAlive(nodeId) &&
               removedMarks_.isMarked(static_cast<size_t>(nodeId)) && topologyOf(dualTree).getNumProperParts(nodeId) == 0;
    }

    /**
     * @brief Contracts one empty non-root node into its current parent.
     *
     * @param dualTree Dual tree used by the operation.
     * @param editor Active tree editor used by the operation.
     * @param removedNodeId Node identifier represented by `removedNodeId`.
     * @return The surviving parent empty non-root node into its current parent.
     */
    NodeId absorbRemovedNonRootNode(tree_t* dualTree, editor_t& editor, NodeId removedNodeId) {
        const MorphologicalTree& topology = topologyOf(dualTree);
        const NodeId parentId = topology.getNodeParent(removedNodeId);
        if (parentId == InvalidNode || parentId == removedNodeId || !topology.isAlive(parentId)) {
            return InvalidNode;
        }

        const bool parentChanged = topology.getFirstChild(removedNodeId) != InvalidNode;
        editor.moveChildren(parentId, removedNodeId);
        notifyMoveProperParts(dualTree, parentId, removedNodeId);
        editor.moveProperParts(parentId, removedNodeId);
        disconnect(dualTree, editor, removedNodeId, true);
        if (parentChanged) {
            markAttributeUpdate(dualTree, parentId);
        }
        computeAttributeOnTreeNode(dualTree, parentId);

        if (topologyOf(dualTree).isAlive(parentId) && topologyOf(dualTree).getNumProperParts(parentId) == 0) {
            removedMarks_.mark(static_cast<size_t>(parentId));
            return parentId;
        }
        return InvalidNode;
    }

    /**
     * @brief Contracts an empty root by promoting the altitude-compatible child.
     *
     * @param dualTree Dual tree used by the operation.
     * @param editor Active tree editor used by the operation.
     * @param removedNodeId Node identifier represented by `removedNodeId`.
     */
    void absorbRemovedRootNode(tree_t* dualTree, editor_t& editor, NodeId removedNodeId) {
        const bool isMaxtree = dualTree == maxtree_;
        const MorphologicalTree& topology = topologyOf(dualTree);
        const NodeId firstChild = topology.getFirstChild(removedNodeId);
        if (firstChild == InvalidNode) {
            return;
        }

        NodeId newRoot = firstChild;
        for (NodeId childId = firstChild; childId != InvalidNode; childId = topologyOf(dualTree).getNextSibling(childId)) {
            if ((isMaxtree && nodeAltitude(dualTree, childId) < nodeAltitude(dualTree, newRoot)) ||
                (!isMaxtree && nodeAltitude(dualTree, childId) > nodeAltitude(dualTree, newRoot))) {
                newRoot = childId;
            }
        }

        bool rootChanged = false;
        for (NodeId childId = firstChild; childId != InvalidNode;) {
            const NodeId next = topologyOf(dualTree).getNextSibling(childId);
            if (childId != newRoot && !topologyOf(dualTree).hasChild(newRoot, childId)) {
                editor.detach(childId);
                editor.attach(newRoot, childId);
                rootChanged = true;
            }
            childId = next;
        }

        editor.setRoot(newRoot);
        editor.releaseNode(removedNodeId);
        notifyNodeRemoved(dualTree, removedNodeId);
        if (rootChanged) {
            markAttributeUpdate(dualTree, newRoot);
        }
        computeAttributeOnTreeNode(dualTree, newRoot);
    }

    /**
     * @brief Contracts all still-empty marked nodes in post-order.
     *
     * @param dualTree Dual tree used by the operation.
     * @param editor Active tree editor used by the operation.
     * @param removedNodeIds Node identifier represented by `removedNodeIds`.
     */
    void absorbRemovedNodes(tree_t* dualTree, editor_t& editor, const std::vector<NodeId>& removedNodeIds) {
        struct Frame {
            NodeId nodeId = InvalidNode;
            NodeId nextChildId = InvalidNode;
        };

        if (dualTree == nullptr || removedNodeIds.empty()) {
            return;
        }

        std::vector<Frame> stack;
        stack.reserve(std::max<std::size_t>(64, removedNodeIds.size()));
        const auto makeFrame = [dualTree](NodeId nodeId) {
            const MorphologicalTree& topology = topologyOf(dualTree);
            if (nodeId == InvalidNode || !topology.isNode(nodeId) || !topology.isAlive(nodeId)) {
                return Frame{nodeId, InvalidNode};
            }
            return Frame{nodeId, topology.getFirstChild(nodeId)};
        };

        for (auto it = removedNodeIds.rbegin(); it != removedNodeIds.rend(); ++it) {
            if (*it != InvalidNode) {
                stack.push_back(makeFrame(*it));
            }
        }

        while (!stack.empty()) {
            Frame& frame = stack.back();
            NodeId currentNodeId = InvalidNode;

            if (!canAbsorbRemovedNode(dualTree, frame.nodeId)) {
                stack.pop_back();
            } else if (frame.nextChildId != InvalidNode) {
                const NodeId childId = frame.nextChildId;
                frame.nextChildId = topologyOf(dualTree).getNextSibling(childId);
                stack.push_back(makeFrame(childId));
            } else {
                currentNodeId = frame.nodeId;
                stack.pop_back();
            }

            if (canAbsorbRemovedNode(dualTree, currentNodeId)) {
                if (topologyOf(dualTree).isRoot(currentNodeId)) {
                    absorbRemovedRootNode(dualTree, editor, currentNodeId);
                } else {
                    const NodeId parentId = absorbRemovedNonRootNode(dualTree, editor, currentNodeId);
                    if (parentId != InvalidNode) {
                        stack.push_back(makeFrame(parentId));
                    }
                }
            }
        }
    }

    /**
     * @brief Collapses a chain of empty marked nodes immediately below the root.
     *
     * Root contraction is special because the surviving child must become the
     * root candidate. This helper repeatedly promotes the altitude-compatible
     * child while preserving all other grandchildren under the promoted node.
     *
     * @param dualTree Dual tree used by the operation.
     * @param editor Active tree editor used by the operation.
     * @param rootId Identifier of the root node.
     * @param childId Identifier of the child node.
     * @return The surviving node chain of empty marked nodes immediately below the root.
     */
    NodeId collapseRemovedRootBranch(tree_t* dualTree, editor_t& editor, NodeId rootId, NodeId childId) {
        assert(dualTree != nullptr);
        const bool isMaxtree = dualTree == maxtree_;

        NodeId current = childId;
        while (current != InvalidNode && topologyOf(dualTree).isNode(current) && topologyOf(dualTree).isAlive(current) &&
               topologyOf(dualTree).getNodeParent(current) == rootId && removedMarks_.isMarked(static_cast<size_t>(current)) &&
               topologyOf(dualTree).getNumProperParts(current) == 0) {
            const NodeId firstGrandchild = topologyOf(dualTree).getFirstChild(current);
            if (firstGrandchild == InvalidNode) {
                break;
            }

            NodeId promoted = firstGrandchild;
            for (NodeId grandchildId = firstGrandchild; grandchildId != InvalidNode; grandchildId = topologyOf(dualTree).getNextSibling(grandchildId)) {
                if ((isMaxtree && nodeAltitude(dualTree, grandchildId) < nodeAltitude(dualTree, promoted)) ||
                    (!isMaxtree && nodeAltitude(dualTree, grandchildId) > nodeAltitude(dualTree, promoted))) {
                    promoted = grandchildId;
                }
            }

            if (!topologyOf(dualTree).isRoot(promoted)) {
                disconnect(dualTree, editor, promoted, false);
            }
            editor.attach(rootId, promoted);

            for (NodeId grandchildId = topologyOf(dualTree).getFirstChild(current); grandchildId != InvalidNode;) {
                const NodeId next = topologyOf(dualTree).getNextSibling(grandchildId);
                if (grandchildId != promoted && !topologyOf(dualTree).hasChild(promoted, grandchildId)) {
                    if (!topologyOf(dualTree).isRoot(grandchildId)) {
                        disconnect(dualTree, editor, grandchildId, false);
                    }
                    editor.attach(promoted, grandchildId);
                }
                grandchildId = next;
            }

            disconnect(dualTree, editor, current, true);
            markAttributeUpdate(dualTree, promoted);
            computeAttributeOnTreeNode(dualTree, promoted);
            current = promoted;
        }
        return current;
    }

    /**
     * @brief Finishes one dual-tree update after the level sweep.
     *
     * If `nodeCa` was emptied by the move of `C`, the final union node is
     * promoted into its place; otherwise the final union node is reattached
     * below `nodeCa`. The method then contracts every empty node accumulated
     * during the step.
     *
     * @param dualTree Dual tree used by the operation.
     * @param editor Active tree editor used by the operation.
     * @param nodeCa Identifier of the common-ancestor node.
     * @param finalUnionNode Node identifier represented by `finalUnionNode`.
     */
    void finalizeUpdateTreeAndContractRemovedNodes(tree_t* dualTree, editor_t& editor, NodeId nodeCa, NodeId finalUnionNode) {
        if (dualTree == nullptr) {
            return;
        }

        const bool isMaxtree = dualTree == maxtree_;

        if (finalUnionNode != InvalidNode && topologyOf(dualTree).isAlive(finalUnionNode)) {
            if (removedMarks_.isMarked(static_cast<size_t>(nodeCa))) {
                if (!topologyOf(dualTree).isRoot(nodeCa)) {
                    const NodeId nodeCaParentId = topologyOf(dualTree).getNodeParent(nodeCa);
                    bool finalUnionNodeChanged = false;

                    if (!topologyOf(dualTree).isRoot(finalUnionNode)) {
                        disconnect(dualTree, editor, finalUnionNode, false);
                    }
                    editor.attach(nodeCaParentId, finalUnionNode);

                    for (NodeId childId = topologyOf(dualTree).getFirstChild(nodeCa); childId != InvalidNode;) {
                        const NodeId next = topologyOf(dualTree).getNextSibling(childId);
                        if (childId != finalUnionNode && !topologyOf(dualTree).hasChild(finalUnionNode, childId)) {
                            if (!topologyOf(dualTree).isRoot(childId)) {
                                disconnect(dualTree, editor, childId, false);
                            }
                            editor.attach(finalUnionNode, childId);
                            finalUnionNodeChanged = true;
                        }
                        childId = next;
                    }

                    disconnect(dualTree, editor, nodeCa, true);
                    if (finalUnionNodeChanged) {
                        markAttributeUpdate(dualTree, finalUnionNode);
                    }
                    markAttributeUpdate(dualTree, nodeCaParentId);
                    computeAttributeOnTreeNode(dualTree, finalUnionNode);
                    computeAttributeOnTreeNode(dualTree, nodeCaParentId);
                } else {
                    NodeId survivingFinalUnionNode = finalUnionNode;
                    for (NodeId childId = topologyOf(dualTree).getFirstChild(nodeCa); childId != InvalidNode;) {
                        const NodeId next = topologyOf(dualTree).getNextSibling(childId);
                        const NodeId normalizedChild = collapseRemovedRootBranch(dualTree, editor, nodeCa, childId);
                        if (childId == finalUnionNode && normalizedChild != InvalidNode) {
                            survivingFinalUnionNode = normalizedChild;
                        }
                        childId = next;
                    }

                    NodeId candidateRootId = survivingFinalUnionNode;
                    for (NodeId childId = topologyOf(dualTree).getFirstChild(nodeCa); childId != InvalidNode;
                         childId = topologyOf(dualTree).getNextSibling(childId)) {
                        if ((isMaxtree && nodeAltitude(dualTree, childId) < nodeAltitude(dualTree, candidateRootId)) ||
                            (!isMaxtree && nodeAltitude(dualTree, childId) > nodeAltitude(dualTree, candidateRootId))) {
                            candidateRootId = childId;
                        }
                    }

                    bool candidateRootChanged = false;
                    if (candidateRootId != survivingFinalUnionNode) {
                        if (!topologyOf(dualTree).isRoot(survivingFinalUnionNode)) {
                            disconnect(dualTree, editor, survivingFinalUnionNode, false);
                        }
                        editor.attach(candidateRootId, survivingFinalUnionNode);
                        candidateRootChanged = true;
                    }

                    for (NodeId childId = topologyOf(dualTree).getFirstChild(nodeCa); childId != InvalidNode;) {
                        const NodeId next = topologyOf(dualTree).getNextSibling(childId);
                        if (childId != candidateRootId && !topologyOf(dualTree).hasChild(candidateRootId, childId)) {
                            if (!topologyOf(dualTree).isRoot(childId)) {
                                disconnect(dualTree, editor, childId, false);
                            }
                            editor.attach(candidateRootId, childId);
                            candidateRootChanged = true;
                        }
                        childId = next;
                    }

                    editor.setRoot(candidateRootId);
                    editor.releaseNode(nodeCa);
                    notifyNodeRemoved(dualTree, nodeCa);
                    if (candidateRootChanged) {
                        markAttributeUpdate(dualTree, candidateRootId);
                    }
                    computeAttributeOnTreeNode(dualTree, candidateRootId);
                }
            } else {
                if (!topologyOf(dualTree).isRoot(finalUnionNode)) {
                    disconnect(dualTree, editor, finalUnionNode, false);
                }
                editor.attach(nodeCa, finalUnionNode);
                markAttributeUpdate(dualTree, nodeCa);
                computeAttributeOnTreeNode(dualTree, nodeCa);
            }
        }

        absorbRemovedNodes(dualTree, editor, removedNodesPendingAbsorption_);
    }

    /**
     * @brief Returns the primal tree corresponding to an update of the given dual tree.
     * @details Updating the max-tree means the primal operation was performed
     * on the min-tree, and conversely.
     *
     * @param isMaxtree Flag controlling is maxtree.
     * @return The primal tree corresponding to an update of the given dual tree.
     */
    tree_t* getPrimalTree(bool isMaxtree) { return isMaxtree ? mintree_ : maxtree_; }

    /**
     * @brief Moves children that are outside the active merge interval.
     *
     * Merge-bucket children are detached because they will be processed by their
     * own level in the sweep. The remaining children are moved to the current
     * union node to preserve the hierarchy around the local edit.
     *
     * @param tree Tree topology used by the operation.
     * @param editor Active tree editor used by the operation.
     * @param targetNodeId Node identifier represented by `targetNodeId`.
     * @param sourceNodeId Node identifier represented by `sourceNodeId`.
     */
    void reattachOutsideIntervalChildren(tree_t* tree, editor_t& editor, NodeId targetNodeId, NodeId sourceNodeId) {
        assert(tree != nullptr);
        if (sourceNodeId == targetNodeId) {
            for (NodeId childId = topologyOf(tree).getFirstChild(sourceNodeId); childId != InvalidNode;) {
                const NodeId next = topologyOf(tree).getNextSibling(childId);
                if (mergeNodesByLevel_.isMergeNode(childId)) {
                    editor.detach(childId);
                }
                childId = next;
            }
            return;
        }

        for (NodeId childId = topologyOf(tree).getFirstChild(sourceNodeId); childId != InvalidNode;) {
            const NodeId next = topologyOf(tree).getNextSibling(childId);
            if (mergeNodesByLevel_.isMergeNode(childId)) {
                editor.detach(childId);
            }
            childId = next;
        }
        editor.moveChildren(targetNodeId, sourceNodeId);
    }

    /**
     * @brief Builds the merge buckets and frontier roots for one update step.
     *
     * Each pixel of `C` contributes graph-neighbor seeds outside `C`. Valid
     * seeds are climbed toward the root while they remain inside the altitude
     * interval induced by `C`. A generation-stamped mark prevents revisiting
     * the same ancestor chain more than once in the step.
     *
     * Nodes whose altitude lies inside the interval are inserted into
     * `mergeNodesByLevel_`. Nodes that leave the interval through the `b`
     * boundary are stored as frontier roots and reattached when the sweep
     * reaches `b`.
     *
     * @param tree Tree topology used by the operation.
     * @param properPartSetC Proper-part data represented by `properPartSetC`.
     * @param nodeCa Identifier of the common-ancestor node.
     * @param b Branch or side selector used by the operation.
     * @param isMaxtree Flag controlling is maxtree.
     */
    void buildMergedAndNestedCollections(const tree_t& tree, const std::vector<NodeId>& properPartSetC, NodeId nodeCa, altitude_t b, bool isMaxtree) {
        mergeNodesByLevel_.resetCollection(isMaxtree);
        assert(nodeCa != InvalidNode);
        climbedNodeMarks_.resetAll();
        const altitude_t altitudeCa = nodeAltitude(&tree, nodeCa);

        for (NodeId p : properPartSetC) {
            for (NodeId q : graph_->getNeighborIndices(p)) {
                if (pixelsInCMarks_.isMarked(static_cast<size_t>(q))) {
                    continue;
                }

                const NodeId nodeQ = tree.topology().getProperPartOwner(q);
                if (nodeQ == InvalidNode) {
                    continue;
                }

                const altitude_t altitudeQ = nodeAltitude(&tree, nodeQ);
                const bool validSeed = (isMaxtree && altitudeQ >= altitudeCa) || (!isMaxtree && altitudeQ <= altitudeCa);
                if (!validSeed) {
                    continue;
                }

                if (mergeNodesByLevel_.markAdjacentSeed(nodeQ)) {
                    NodeId nodeSubtree = nodeQ;
                    NodeId n = nodeQ;
                    while (n != InvalidNode && tree.topology().isAlive(n) && !climbedNodeMarks_.isMarked(static_cast<size_t>(n))) {
                        const altitude_t levelCurrent = nodeAltitude(&tree, n);
                        if (!((isMaxtree && levelCurrent >= altitudeCa) || (!isMaxtree && levelCurrent <= altitudeCa))) {
                            break;
                        }

                        climbedNodeMarks_.mark(static_cast<size_t>(n));
                        nodeSubtree = n;

                        if ((isMaxtree && levelCurrent <= b) || (!isMaxtree && levelCurrent >= b)) {
                            mergeNodesByLevel_.addMergeNode(tree, nodeSubtree);
                        } else {
                            NodeId parentId = tree.topology().getNodeParent(nodeSubtree);
                            if (parentId == nodeSubtree) {
                                parentId = InvalidNode;
                            }
                            if (!(parentId != InvalidNode &&
                                  ((isMaxtree && nodeAltitude(&tree, parentId) > b) || (!isMaxtree && nodeAltitude(&tree, parentId) < b)))) {
                                mergeNodesByLevel_.addFrontierNodeAboveB(nodeSubtree);
                            }
                        }

                        const NodeId parentId = tree.topology().getNodeParent(n);
                        if (parentId == n) {
                            break;
                        }
                        n = parentId;
                    }
                }
            }
        }
    }

    /**
     * @brief Updates the dual tree after a rooted subtree is removed in the primal tree.
     *
     * Detailed step sequence:
     *
     * - collect all proper parts of the primal subtree into `properPartSetC_`;
     * - find `nodeCa`, the extremal dual node currently owning pixels of `C`;
     * - collect adjacent seeds and level buckets around `C`;
     * - sweep active levels from the outside toward `altitudeCa_`, merging each
     *   non-empty bucket into a single current union node;
     * - at level `b`, move the selected proper parts and attach frontier roots;
     * - reattach the previous level union node under the current one;
     * - finalize by promoting/reattaching the final union node and contracting
     *   nodes emptied by the update.
     *
     * All tree mutations are performed through `WeightedTreeEditor<T>`. The
     * generic move-only proof at the end records that this algorithm
     * established the edit invariants during its existing passes and avoids
     * duplicate global validation in the Release hot loop.
     *
     * @param dualTree Dual tree used by the operation.
     * @param subtreeRoot Root of the subtree being processed.
     */
    void updateTree(tree_t* dualTree, NodeId subtreeRoot) {
        assert(dualTree != nullptr);
        assert(subtreeRoot != InvalidNode);
        resetAttributeUpdateMarks();

        const bool isMaxtree = dualTree == maxtree_;
        tree_t* primalTree = getPrimalTree(isMaxtree);
        assert(primalTree != nullptr);
        assert(primalTree->topology().isNode(subtreeRoot));
        assert(primalTree->topology().isAlive(subtreeRoot));

        const NodeId subtreeParentId = primalTree->topology().getNodeParent(subtreeRoot);
        assert(subtreeParentId != InvalidNode && subtreeParentId != subtreeRoot);
        const altitude_t b = nodeAltitude(primalTree, subtreeParentId);

        NodeId nodeCa = InvalidNode;
        altitudeCa_ = altitude_t{};

        properPartSetC_.clear();
        properPartSetC_.reserve(64);
        pixelsInCMarks_.resetAll();
        for (NodeId subtreeNodeId : primalTree->topology().getNodeSubtree(subtreeRoot)) {
            for (NodeId p : primalTree->topology().getProperParts(subtreeNodeId)) {
                properPartSetC_.push_back(p);
                pixelsInCMarks_.mark(static_cast<size_t>(p));

                const NodeId nodeP = dualTree->topology().getProperPartOwner(p);
                if (nodeP == InvalidNode) {
                    continue;
                }
                const altitude_t altitudeP = nodeAltitude(dualTree, nodeP);
                if (nodeCa == InvalidNode || ((isMaxtree && altitudeP < altitudeCa_) || (!isMaxtree && altitudeP > altitudeCa_))) {
                    altitudeCa_ = altitudeP;
                    nodeCa = nodeP;
                }
            }
        }
        if (properPartSetC_.empty()) {
            return;
        }
        assert(nodeCa != InvalidNode);

        buildMergedAndNestedCollections(*dualTree, properPartSetC_, nodeCa, b, isMaxtree);

        editor_t editor = mmcfilters::detail::beginEstablishedWeightedEdit(*dualTree);
        altitude_t currentMergeLevel = mergeNodesByLevel_.firstMergeLevel();
        NodeId currentUnionNode = InvalidNode;
        NodeId previousLevelUnionNode = InvalidNode;
        removedMarks_.resetAll();
        removedNodesPendingAbsorption_.clear();
        nodesPendingRemoval_.reserve(mergeNodesByLevel_.getMaxBucketSize());

        while (mergeNodesByLevel_.hasMergeLevel() && ((isMaxtree && currentMergeLevel > altitudeCa_) || (!isMaxtree && currentMergeLevel < altitudeCa_))) {
            auto& nodesAtCurrentLevel = mergeNodesByLevel_.getMergedNodes(currentMergeLevel);
            currentUnionNode = InvalidNode;
            nodesPendingRemoval_.clear();

            for (NodeId nodeId : nodesAtCurrentLevel) {
                if (!dualTree->topology().isAlive(nodeId)) {
                    continue;
                }

                if (currentUnionNode == InvalidNode) {
                    if (removedMarks_.isMarked(static_cast<size_t>(nodeId))) {
                        nodesPendingRemoval_.push_back(nodeId);
                        continue;
                    }
                    currentUnionNode = nodeId;
                    disconnect(dualTree, editor, currentUnionNode, false);
                    reattachOutsideIntervalChildren(dualTree, editor, currentUnionNode, currentUnionNode);

                    for (NodeId pendingNodeId : nodesPendingRemoval_) {
                        reattachOutsideIntervalChildren(dualTree, editor, currentUnionNode, pendingNodeId);
                        notifyMoveProperParts(dualTree, currentUnionNode, pendingNodeId);
                        editor.moveProperParts(currentUnionNode, pendingNodeId);
                        disconnect(dualTree, editor, pendingNodeId, true);
                    }
                    nodesPendingRemoval_.clear();
                    continue;
                }

                reattachOutsideIntervalChildren(dualTree, editor, currentUnionNode, nodeId);
                notifyMoveProperParts(dualTree, currentUnionNode, nodeId);
                editor.moveProperParts(currentUnionNode, nodeId);
                disconnect(dualTree, editor, nodeId, true);
            }

            if (currentUnionNode == InvalidNode) {
                for (NodeId nodeId : nodesPendingRemoval_) {
                    if (!dualTree->topology().isAlive(nodeId) || dualTree->topology().isRoot(nodeId)) {
                        continue;
                    }
                    const NodeId parentId = dualTree->topology().getNodeParent(nodeId);
                    editor.moveChildren(parentId, nodeId);
                    notifyMoveProperParts(dualTree, parentId, nodeId);
                    editor.moveProperParts(parentId, nodeId);
                    disconnect(dualTree, editor, nodeId, true);
                }
                currentMergeLevel = mergeNodesByLevel_.nextMergeLevel();
                currentUnionNode = previousLevelUnionNode;
                continue;
            }

            if (currentMergeLevel == b) {
                moveSelectedProperPartsToNode(dualTree, editor, currentUnionNode, properPartSetC_);
                for (NodeId nodeId : mergeNodesByLevel_.getFrontierNodesAboveB()) {
                    disconnect(dualTree, editor, nodeId, false);
                    editor.attach(currentUnionNode, nodeId);
                }
            }

            if (previousLevelUnionNode != InvalidNode) {
                if (dualTree->topology().isAlive(previousLevelUnionNode) && !dualTree->topology().hasChild(currentUnionNode, previousLevelUnionNode)) {
                    if (!dualTree->topology().isRoot(previousLevelUnionNode)) {
                        disconnect(dualTree, editor, previousLevelUnionNode, false);
                    }
                    editor.attach(currentUnionNode, previousLevelUnionNode);
                }
            }

            markAttributeUpdate(dualTree, currentUnionNode);
            computeAttributeOnTreeNode(dualTree, currentUnionNode);

            previousLevelUnionNode = currentUnionNode;
            currentMergeLevel = mergeNodesByLevel_.nextMergeLevel();
        }

        finalizeUpdateTreeAndContractRemovedNodes(dualTree, editor, nodeCa, previousLevelUnionNode);
        auto proof = editor.proveIncremental();
        editor.commit(std::move(proof));
    }

  public:
    /**
     * @brief Reports whether this altitude type uses dense per-level buckets.
     *
     * @return True if this altitude type uses dense per-level buckets; otherwise false.
     */
    static constexpr bool usesDenseLevelBackend() { return use_dense_levels; }

    /**
     * @brief Returns the bit-width threshold for dense bucket selection.
     *
     * @return The bit-width threshold for dense bucket selection.
     */
    static constexpr int denseLevelBackendMaxBits() { return MMCFILTERS_COMPONENT_TREE_ADJUSTMENT_DENSE_MAX_BITS; }

    /**
     * @brief Creates an adjuster for externally owned min/max-tree states.
     * @param mintree Mutable min-tree state.
     * @param maxtree Mutable max-tree state.
     * @param graph Shared image-domain adjacency relation.
     *
     * The constructor sizes all generation-stamped mark sets from the current
     * node/proper-part spaces of the two trees. The tree topology is expected to
     * keep these slot spaces stable under local releases performed by the
     * adjuster.
     */
    DualMinMaxTreeIncrementalFilter(tree_t* mintree, tree_t* maxtree, const RegularGridAdjacency2D& graph)
        : mintree_(mintree), maxtree_(maxtree), graph_(&graph), mergeNodesByLevel_(std::max(mintree ? mintree->topology().getNumInternalNodeSlots() : 0,
                                                                                            maxtree ? maxtree->topology().getNumInternalNodeSlots() : 0)),
          removedMarks_(static_cast<size_t>(
              std::max(mintree ? mintree->topology().getNumInternalNodeSlots() : 0, maxtree ? maxtree->topology().getNumInternalNodeSlots() : 0))),
          pixelsInCMarks_(static_cast<size_t>(
              std::max(mintree ? mintree->topology().getNumTotalProperParts() : 0, maxtree ? maxtree->topology().getNumTotalProperParts() : 0))),
          climbedNodeMarks_(static_cast<size_t>(
              std::max(mintree ? mintree->topology().getNumInternalNodeSlots() : 0, maxtree ? maxtree->topology().getNumInternalNodeSlots() : 0))),
          attributeUpdateMarks_(static_cast<size_t>(
              std::max(mintree ? mintree->topology().getNumInternalNodeSlots() : 0, maxtree ? maxtree->topology().getNumInternalNodeSlots() : 0))) {
        assert(mintree_ != nullptr);
        assert(maxtree_ != nullptr);
        assert(graph_ != nullptr);
    }

    /**
     * @brief Registers incremental attribute computers and their external buffers.
     *
     * The buffers are not owned by the adjuster. They must remain alive and
     * indexed by internal `NodeId` while pruning/update calls are executed.
     *
     * @param computerMin Flag controlling computer min.
     * @param computerMax Flag controlling computer max.
     * @param bufferMin Min-tree attribute buffer.
     * @param bufferMax Max-tree attribute buffer.
     */
    void setAttributeComputer(const attribute_computer_t& computerMin, const attribute_computer_t& computerMax, std::vector<double>& bufferMin,
                              std::vector<double>& bufferMax) {
        attributeComputerMin_ = &computerMin;
        attributeComputerMax_ = &computerMax;
        attributeBufferMin_ = &bufferMin;
        attributeBufferMax_ = &bufferMax;
    }

    /**
     * @brief Registers one shared incremental attribute computer for both trees.
     *
     * @param computer Flag controlling computer.
     * @param bufferMin Min-tree attribute buffer.
     * @param bufferMax Max-tree attribute buffer.
     */
    void setAttributeComputer(const attribute_computer_t& computer, std::vector<double>& bufferMin, std::vector<double>& bufferMax) {
        setAttributeComputer(computer, computer, bufferMin, bufferMax);
    }

    /**
     * @brief Prunes max-tree subtrees and updates the min-tree before each prune.
     *
     * Invalid, dead, and root nodes are ignored. For each valid max-tree root,
     * the min-tree is updated first against the subtree that is about to be
     * removed; then the max-tree subtree is pruned.
     *
     * @param nodesToPrune Node identifiers selected for pruning.
     */
    void pruneMaxTreeAndUpdateMinTree(const std::vector<NodeId>& nodesToPrune) {
        assert(mintree_ != nullptr);
        assert(maxtree_ != nullptr);
        for (NodeId rootSubtree : nodesToPrune) {
            if (rootSubtree == InvalidNode || rootSubtree == maxtree_->topology().getRoot() || !maxtree_->topology().isNode(rootSubtree) ||
                !maxtree_->topology().isAlive(rootSubtree)) {
                continue;
            }
            updateTree(mintree_, rootSubtree);
            maxtree_->pruneNode(rootSubtree);
        }
    }

    /**
     * @brief Prunes min-tree subtrees and updates the max-tree before each prune.
     *
     * Invalid, dead, and root nodes are ignored. For each valid min-tree root,
     * the max-tree is updated first against the subtree that is about to be
     * removed; then the min-tree subtree is pruned.
     *
     * @param nodesToPrune Node identifiers selected for pruning.
     */
    void pruneMinTreeAndUpdateMaxTree(const std::vector<NodeId>& nodesToPrune) {
        assert(mintree_ != nullptr);
        assert(maxtree_ != nullptr);
        for (NodeId rootSubtree : nodesToPrune) {
            if (rootSubtree == InvalidNode || rootSubtree == mintree_->topology().getRoot() || !mintree_->topology().isNode(rootSubtree) ||
                !mintree_->topology().isAlive(rootSubtree)) {
                continue;
            }
            updateTree(maxtree_, rootSubtree);
            mintree_->pruneNode(rootSubtree);
        }
    }
};

} // namespace mmcfilters::adjust
