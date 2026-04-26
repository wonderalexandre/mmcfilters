#pragma once

#include "DynamicTreeAttributeComputer.hpp"
#include "../WeightedMorphologicalTree.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <type_traits>
#include <vector>

#ifndef MMCFILTERS_COMPONENT_TREE_ADJUSTMENT_DENSE_MAX_BITS
#define MMCFILTERS_COMPONENT_TREE_ADJUSTMENT_DENSE_MAX_BITS 8
#endif

namespace mmcfilters::adjust {

/**
 * @brief Incremental updater for a paired min-tree / max-tree state.
 *
 * This class ports the Higra dual component-tree adjustment algorithm to this
 * repository's `WeightedMorphologicalTree` representation. The algorithmic
 * structure is intentionally kept the same:
 *
 * - a rooted subtree is selected for pruning in the primal tree;
 * - its full support defines the proper-part set `C`;
 * - graph-adjacent components around `C` are collected in the dual tree;
 * - collected nodes are bucketed by altitude and swept in tree order;
 * - local unions, frontier reattachments, and empty-node contractions update
 *   the dual tree without rebuilding it globally.
 *
 * Representation notes:
 *
 * - node ids are this project's dense internal node ids, not Higra global ids;
 * - altitudes are read from `WeightedMorphologicalTree::getAltitude`;
 * - low-level edits go through `WeightedTreeEditor` so the topology remains
 *   validated at the end of each update step;
 * - optional attribute buffers are kept in sync through
 *   `DynamicTreeAttributeComputer`.
 *
 * The expected cost remains local to the affected region. The dominant terms
 * are the same as in the Higra implementation: collecting `C`, visiting
 * adjacent seeds, climbing relevant ancestor chains once, sweeping the active
 * altitude buckets, and contracting the nodes actually emptied by the step.
 */
template<typename altitude_t = AltitudeType>
class DualMinMaxTreeIncrementalFilter {
private:
    using tree_t = WeightedMorphologicalTree;

    static constexpr bool usesDenseLevels() {
        if constexpr (std::is_integral_v<altitude_t>) {
            using unsigned_altitude_t = std::make_unsigned_t<altitude_t>;
            return std::numeric_limits<unsigned_altitude_t>::digits <= MMCFILTERS_COMPONENT_TREE_ADJUSTMENT_DENSE_MAX_BITS;
        } else {
            return false;
        }
    }

    static constexpr bool use_dense_levels = usesDenseLevels();

    /**
     * @brief Per-step collection of merge buckets and frontier roots.
     *
     * Nodes adjacent to `C` are grouped by altitude for the subsequent sweep.
     * Small integral altitude domains use a dense bucket array; larger integral
     * domains and floating-point instantiations use a sparse ordered map. This
     * matches the Higra backend choice and avoids imposing an 8-bit altitude
     * domain on this project, where `AltitudeType` is currently `int`.
     */
    class MergedNodesCollection {
    private:
        template<bool dense, typename altitude_type>
        struct StorageSelector;

        template<typename altitude_type>
        struct StorageSelector<true, altitude_type> {
            static constexpr std::size_t domain_size =
                static_cast<std::size_t>(
                    static_cast<long long>(std::numeric_limits<altitude_type>::max()) -
                    static_cast<long long>(std::numeric_limits<altitude_type>::lowest()) + 1);
            using type = std::array<std::vector<NodeId>, domain_size>;
        };

        template<typename altitude_type>
        struct StorageSelector<false, altitude_type> {
            using type = std::map<altitude_type, std::vector<NodeId>>;
        };

        using storage_t = typename StorageSelector<use_dense_levels, altitude_t>::type;

        storage_t mergeNodesByLevelStorage_;
        std::vector<altitude_t> mergeLevels_;
        std::vector<NodeId> frontierNodesAboveB_;
        GenerationStampSet collectedNodeMarks_;
        GenerationStampSet mergeBucketNodeMarks_;
        GenerationStampSet adjacentSeedMarks_;
        std::size_t maxBucketSize_ = 0;
        int currentMergeLevelIndex_ = 0;
        bool isMaxtree_ = false;

    public:
        static std::size_t denseBucketIndex(altitude_t level) {
            if constexpr (std::is_signed_v<altitude_t>) {
                return static_cast<std::size_t>(
                    static_cast<long long>(level) -
                    static_cast<long long>(std::numeric_limits<altitude_t>::lowest()));
            } else {
                return static_cast<std::size_t>(level);
            }
        }

        static altitude_t levelFromDenseBucketIndex(std::size_t index) {
            if constexpr (std::is_signed_v<altitude_t>) {
                return static_cast<altitude_t>(
                    static_cast<long long>(index) +
                    static_cast<long long>(std::numeric_limits<altitude_t>::lowest()));
            } else {
                return static_cast<altitude_t>(index);
            }
        }

        /**
         * @brief Creates temporary marks sized for the dense node-id domain.
         */
        explicit MergedNodesCollection(int maxNodes = 0)
            : collectedNodeMarks_(static_cast<size_t>(std::max(maxNodes, 0))),
              mergeBucketNodeMarks_(static_cast<size_t>(std::max(maxNodes, 0))),
              adjacentSeedMarks_(static_cast<size_t>(std::max(maxNodes, 0))) {}

        /**
         * @brief Clears transient state before one adjustment step.
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
            collectedNodeMarks_.resetAll();
            mergeBucketNodeMarks_.resetAll();
            adjacentSeedMarks_.resetAll();
            maxBucketSize_ = 0;
            currentMergeLevelIndex_ = 0;
        }

        std::vector<NodeId>& getMergedNodes(const altitude_t& level) {
            if constexpr (use_dense_levels) {
                return mergeNodesByLevelStorage_[denseBucketIndex(level)];
            } else {
                return mergeNodesByLevelStorage_[level];
            }
        }

        std::vector<NodeId>& getFrontierNodesAboveB() { return frontierNodesAboveB_; }
        std::size_t getMaxBucketSize() const { return maxBucketSize_; }

        bool markAdjacentSeed(NodeId nodeId) {
            if (adjacentSeedMarks_.isMarked(static_cast<size_t>(nodeId))) {
                return false;
            }
            adjacentSeedMarks_.mark(static_cast<size_t>(nodeId));
            return true;
        }

        void addFrontierNodeAboveB(NodeId nodeId) {
            if (!collectedNodeMarks_.isMarked(static_cast<size_t>(nodeId))) {
                frontierNodesAboveB_.push_back(nodeId);
                collectedNodeMarks_.mark(static_cast<size_t>(nodeId));
            }
        }

        /**
         * @brief Adds a node to the bucket corresponding to its own altitude.
         */
        void addMergeNode(const tree_t& tree, NodeId nodeId) {
            if (collectedNodeMarks_.isMarked(static_cast<size_t>(nodeId))) {
                return;
            }
            auto& bucket = getMergedNodes(static_cast<altitude_t>(tree.getAltitude(nodeId)));
            bucket.push_back(nodeId);
            maxBucketSize_ = std::max(maxBucketSize_, bucket.size());
            collectedNodeMarks_.mark(static_cast<size_t>(nodeId));
            mergeBucketNodeMarks_.mark(static_cast<size_t>(nodeId));
        }

        bool isMergeNode(NodeId nodeId) const {
            return mergeBucketNodeMarks_.isMarked(static_cast<size_t>(nodeId));
        }

        /**
         * @brief Materializes the ordered list of active levels and returns the first one.
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
            return mergeLevels_[static_cast<size_t>(currentMergeLevelIndex_)];
        }

        bool hasMergeLevel() const {
            return !mergeLevels_.empty() &&
                   currentMergeLevelIndex_ >= 0 &&
                   currentMergeLevelIndex_ < static_cast<int>(mergeLevels_.size());
        }

        altitude_t nextMergeLevel() {
            currentMergeLevelIndex_ = isMaxtree_ ? currentMergeLevelIndex_ - 1 : currentMergeLevelIndex_ + 1;
            if (!hasMergeLevel()) {
                return altitude_t{};
            }
            return mergeLevels_[static_cast<size_t>(currentMergeLevelIndex_)];
        }
    };

    tree_t* mintree_ = nullptr;
    tree_t* maxtree_ = nullptr;
    AdjacencyRelation* graph_ = nullptr;

    MergedNodesCollection mergeNodesByLevel_;
    GenerationStampSet removedMarks_;
    GenerationStampSet pixelsInCMarks_;
    GenerationStampSet climbedNodeMarks_;
    GenerationStampSet attributeUpdateMarks_;
    std::vector<NodeId> properPartSetC_;
    std::vector<NodeId> nodesPendingRemoval_;
    std::vector<NodeId> removedNodesPendingAbsorption_;
    altitude_t altitudeCa_ = altitude_t{};

    const DynamicTreeAttributeComputer* attributeComputerMin_ = nullptr;
    const DynamicTreeAttributeComputer* attributeComputerMax_ = nullptr;
    std::vector<double>* attributeBufferMin_ = nullptr;
    std::vector<double>* attributeBufferMax_ = nullptr;

    static const MorphologicalTree& topologyOf(const tree_t* tree) {
        assert(tree != nullptr);
        return tree->topology();
    }

    DynamicTreeAttributeComputer::buffer_type* getAttributeBuffer(tree_t* tree) const {
        const auto* computer = tree == maxtree_ ? attributeComputerMax_ : attributeComputerMin_;
        if (computer == nullptr) {
            return nullptr;
        }
        return tree == maxtree_ ? attributeBufferMax_ : attributeBufferMin_;
    }

    const DynamicTreeAttributeComputer* getAttributeComputer(tree_t* tree) const {
        if (tree == nullptr) {
            return nullptr;
        }
        return tree == maxtree_ ? attributeComputerMax_ : attributeComputerMin_;
    }

    void notifyNodeRemoved(tree_t* tree, NodeId nodeId) const {
        const auto* computer = getAttributeComputer(tree);
        if (computer != nullptr && tree != nullptr && nodeId != InvalidNode) {
            computer->onNodeRemoved(nodeId, *tree);
        }
    }

    void notifyMoveProperParts(tree_t* tree, NodeId targetNodeId, NodeId sourceNodeId) const {
        const auto* computer = getAttributeComputer(tree);
        if (computer != nullptr && tree != nullptr) {
            computer->onMoveProperParts(targetNodeId, sourceNodeId, *tree);
        }
    }

    void notifyMoveProperPart(tree_t* tree, NodeId targetNodeId, NodeId sourceNodeId, NodeId pixelId) const {
        const auto* computer = getAttributeComputer(tree);
        if (computer != nullptr && tree != nullptr) {
            computer->onMoveProperPart(targetNodeId, sourceNodeId, pixelId, *tree);
        }
    }

    altitude_t nodeAltitude(const tree_t* tree, NodeId nodeId) const {
        assert(tree != nullptr);
        return static_cast<altitude_t>(tree->getAltitude(nodeId));
    }

    /**
     * @brief Refreshes one marked node attribute after local topology stabilization.
     */
    void computeAttributeOnTreeNode(tree_t* tree, NodeId nodeId) {
        const auto* computer = getAttributeComputer(tree);
        if (computer == nullptr || tree == nullptr || nodeId == InvalidNode ||
            !topologyOf(tree).isNode(nodeId) || !topologyOf(tree).isAlive(nodeId)) {
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

    void resetAttributeUpdateMarks() {
        attributeUpdateMarks_.resetAll();
        altitudeCa_ = altitude_t{};
    }

    void markAttributeUpdate(tree_t* tree, NodeId nodeId) {
        if (getAttributeComputer(tree) == nullptr || tree == nullptr || nodeId == InvalidNode) {
            return;
        }
        if (!topologyOf(tree).isNode(nodeId) || !topologyOf(tree).isAlive(nodeId)) {
            return;
        }
        attributeUpdateMarks_.mark(static_cast<size_t>(nodeId));
    }

    void disconnect(tree_t* tree, WeightedTreeEditor& editor, NodeId nodeId, bool releaseNode) {
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
     */
    void moveSelectedProperPartsToNode(tree_t* dualTree, WeightedTreeEditor& editor, NodeId unionNode, const std::vector<NodeId>& properPartSetC) {
        for (NodeId pixelId : properPartSetC) {
            const NodeId ownerId = topologyOf(dualTree).getSmallestComponent(pixelId);
            if (ownerId == InvalidNode || ownerId == unionNode) {
                continue;
            }

            notifyMoveProperPart(dualTree, unionNode, ownerId, pixelId);
            editor.moveProperPart(unionNode, ownerId, pixelId);

            if (topologyOf(dualTree).isAlive(ownerId) &&
                topologyOf(dualTree).getNumProperParts(ownerId) == 0 &&
                ownerId != unionNode &&
                !removedMarks_.isMarked(static_cast<size_t>(ownerId))) {
                removedMarks_.mark(static_cast<size_t>(ownerId));
                removedNodesPendingAbsorption_.push_back(ownerId);
            }
        }
    }

    bool canAbsorbRemovedNode(tree_t* dualTree, NodeId nodeId) const {
        return dualTree != nullptr &&
               nodeId != InvalidNode &&
               topologyOf(dualTree).isNode(nodeId) &&
               topologyOf(dualTree).isAlive(nodeId) &&
               removedMarks_.isMarked(static_cast<size_t>(nodeId)) &&
               topologyOf(dualTree).getNumProperParts(nodeId) == 0;
    }

    /**
     * @brief Contracts one empty non-root node into its current parent.
     */
    NodeId absorbRemovedNonRootNode(tree_t* dualTree, WeightedTreeEditor& editor, NodeId removedNodeId) {
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
     */
    void absorbRemovedRootNode(tree_t* dualTree, WeightedTreeEditor& editor, NodeId removedNodeId) {
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
     */
    void absorbRemovedNodes(tree_t* dualTree, WeightedTreeEditor& editor, const std::vector<NodeId>& removedNodeIds) {
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

    NodeId collapseRemovedRootBranch(tree_t* dualTree, WeightedTreeEditor& editor, NodeId rootId, NodeId childId) {
        assert(dualTree != nullptr);
        const bool isMaxtree = dualTree == maxtree_;

        NodeId current = childId;
        while (current != InvalidNode &&
               topologyOf(dualTree).isNode(current) &&
               topologyOf(dualTree).isAlive(current) &&
               topologyOf(dualTree).getNodeParent(current) == rootId &&
               removedMarks_.isMarked(static_cast<size_t>(current)) &&
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

    void finalizeUpdateTreeAndContractRemovedNodes(tree_t* dualTree, WeightedTreeEditor& editor, NodeId nodeCa, NodeId finalUnionNode) {
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
                    for (NodeId childId = topologyOf(dualTree).getFirstChild(nodeCa); childId != InvalidNode; childId = topologyOf(dualTree).getNextSibling(childId)) {
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

    tree_t* getPrimalTree(bool isMaxtree) {
        return isMaxtree ? mintree_ : maxtree_;
    }

    void reattachOutsideIntervalChildren(tree_t* tree, WeightedTreeEditor& editor, NodeId targetNodeId, NodeId sourceNodeId) {
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
     */
    void buildMergedAndNestedCollections(const tree_t& tree, const std::vector<NodeId>& properPartSetC, NodeId nodeCa, altitude_t b, bool isMaxtree) {
        mergeNodesByLevel_.resetCollection(isMaxtree);
        assert(nodeCa != InvalidNode);
        climbedNodeMarks_.resetAll();
        const altitude_t altitudeCa = nodeAltitude(&tree, nodeCa);

        for (NodeId p : properPartSetC) {
            for (NodeId q : graph_->getNeighborPixels(p)) {
                if (pixelsInCMarks_.isMarked(static_cast<size_t>(q))) {
                    continue;
                }

                const NodeId nodeQ = tree.topology().getSmallestComponent(q);
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
                                  ((isMaxtree && nodeAltitude(&tree, parentId) > b) ||
                                   (!isMaxtree && nodeAltitude(&tree, parentId) < b)))) {
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

                const NodeId nodeP = dualTree->topology().getSmallestComponent(p);
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

        WeightedTreeEditor editor = dualTree->edit();
        altitude_t currentMergeLevel = mergeNodesByLevel_.firstMergeLevel();
        NodeId currentUnionNode = InvalidNode;
        NodeId previousLevelUnionNode = InvalidNode;
        removedMarks_.resetAll();
        removedNodesPendingAbsorption_.clear();
        nodesPendingRemoval_.reserve(mergeNodesByLevel_.getMaxBucketSize());

        while (mergeNodesByLevel_.hasMergeLevel() &&
               ((isMaxtree && currentMergeLevel > altitudeCa_) || (!isMaxtree && currentMergeLevel < altitudeCa_))) {
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
                if (dualTree->topology().isAlive(previousLevelUnionNode) &&
                    !dualTree->topology().hasChild(currentUnionNode, previousLevelUnionNode)) {
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
        editor.commitUnchecked();
    }

public:
    /**
     * @brief Reports whether this altitude type uses dense per-level buckets.
     */
    static constexpr bool usesDenseLevelBackend() { return use_dense_levels; }

    /**
     * @brief Returns the bit-width threshold for dense bucket selection.
     */
    static constexpr int denseLevelBackendMaxBits() {
        return MMCFILTERS_COMPONENT_TREE_ADJUSTMENT_DENSE_MAX_BITS;
    }

    /**
     * @brief Creates an adjuster for externally owned min/max-tree states.
     */
    DualMinMaxTreeIncrementalFilter(tree_t* mintree, tree_t* maxtree, AdjacencyRelation& graph)
        : mintree_(mintree),
          maxtree_(maxtree),
          graph_(&graph),
          mergeNodesByLevel_(std::max(mintree ? mintree->topology().getNumInternalNodeSlots() : 0,
                                      maxtree ? maxtree->topology().getNumInternalNodeSlots() : 0)),
          removedMarks_(static_cast<size_t>(std::max(mintree ? mintree->topology().getNumInternalNodeSlots() : 0,
                                                     maxtree ? maxtree->topology().getNumInternalNodeSlots() : 0))),
          pixelsInCMarks_(static_cast<size_t>(std::max(mintree ? mintree->topology().getNumTotalProperParts() : 0,
                                                       maxtree ? maxtree->topology().getNumTotalProperParts() : 0))),
          climbedNodeMarks_(static_cast<size_t>(std::max(mintree ? mintree->topology().getNumInternalNodeSlots() : 0,
                                                         maxtree ? maxtree->topology().getNumInternalNodeSlots() : 0))),
          attributeUpdateMarks_(static_cast<size_t>(std::max(mintree ? mintree->topology().getNumInternalNodeSlots() : 0,
                                                             maxtree ? maxtree->topology().getNumInternalNodeSlots() : 0))) {
        assert(mintree_ != nullptr);
        assert(maxtree_ != nullptr);
        assert(graph_ != nullptr);
    }

    /**
     * @brief Registers incremental attribute computers and their external buffers.
     */
    void setAttributeComputer(const DynamicTreeAttributeComputer& computerMin,
                              const DynamicTreeAttributeComputer& computerMax,
                              std::vector<double>& bufferMin,
                              std::vector<double>& bufferMax) {
        attributeComputerMin_ = &computerMin;
        attributeComputerMax_ = &computerMax;
        attributeBufferMin_ = &bufferMin;
        attributeBufferMax_ = &bufferMax;
    }

    void setAttributeComputer(const DynamicTreeAttributeComputer& computer,
                              std::vector<double>& bufferMin,
                              std::vector<double>& bufferMax) {
        setAttributeComputer(computer, computer, bufferMin, bufferMax);
    }

    /**
     * @brief Prunes max-tree subtrees and updates the min-tree before each prune.
     */
    void pruneMaxTreeAndUpdateMinTree(const std::vector<NodeId>& nodesToPrune) {
        assert(mintree_ != nullptr);
        assert(maxtree_ != nullptr);
        for (NodeId rootSubtree : nodesToPrune) {
            if (rootSubtree == InvalidNode ||
                rootSubtree == maxtree_->topology().getRoot() ||
                !maxtree_->topology().isNode(rootSubtree) ||
                !maxtree_->topology().isAlive(rootSubtree)) {
                continue;
            }
            updateTree(mintree_, rootSubtree);
            maxtree_->pruneNode(rootSubtree);
        }
    }

    /**
     * @brief Prunes min-tree subtrees and updates the max-tree before each prune.
     */
    void pruneMinTreeAndUpdateMaxTree(const std::vector<NodeId>& nodesToPrune) {
        assert(mintree_ != nullptr);
        assert(maxtree_ != nullptr);
        for (NodeId rootSubtree : nodesToPrune) {
            if (rootSubtree == InvalidNode ||
                rootSubtree == mintree_->topology().getRoot() ||
                !mintree_->topology().isNode(rootSubtree) ||
                !mintree_->topology().isAlive(rootSubtree)) {
                continue;
            }
            updateTree(maxtree_, rootSubtree);
            mintree_->pruneNode(rootSubtree);
        }
    }
};

} // namespace mmcfilters::adjust
