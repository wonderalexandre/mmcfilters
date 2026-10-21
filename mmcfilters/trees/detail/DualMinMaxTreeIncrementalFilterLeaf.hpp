#pragma once

#include "../ValuedMorphologicalTree.hpp"
#include "TreeKernelAccess.hpp"
#include "../../utils/CommittedGridAccess.hpp"
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

namespace mmcfilters::detail {

/**
 * @brief Internal leaf updater for ExtremumResidualTreeByCTsBuilder, adapted from MorphoTreeAdjust.
 *
 * One update levels a non-root primal leaf, updates the dual tree to represent the same image, and contracts the primal leaf.
 * The builder supplies its nonempty support and the unique dual proper part containing it. A whole-part certificate enables a bulk transfer.
 * Both borrowed trees and the adjacency must have matching domains and established node and pixel identifiers.
 *
 * Boundary contacts seed dual ancestor paths. A sweep merges buckets in descending altitude order for a max-tree and ascending order for a min-tree.
 * The support moves at the primal parent altitude. Merge donors are released immediately; the original dual proper-part node is contracted
 * when the merged branch is reattached if its proper part became empty.
 *
 * Small integral altitude types use dense buckets; other types use sparse buckets without quantization.
 * Component minima include descendants. They are initialized once and combined in constant time per merge or attachment.
 * Partial proper-part moves to descendants preserve the minimum of the source component.
 *
 * Candidate nodes, whole proper-part transfers, and the final dual node let the builder update its agenda and pending supports.
 * Established edit sessions use the same kernel reads in CHECKED and UNCHECKED builds; assertion-enabled builds retain the complete edit oracle.
 */
template <AltitudeValue T> class DualMinMaxTreeIncrementalFilterLeaf {
  public:
    /** One complete dual proper-part transfer, in execution order. */
    struct ProperPartTransfer {
        NodeId source; ///< Dual node whose complete proper part was transferred.
        NodeId target; ///< Dual node receiving that proper part.
    };

    /**
     * @brief Borrows a fixed min-tree, max-tree, and adjacency over the same image domain.
     * @param minTree Mutable min-tree owned by the caller.
     * @param maxTree Mutable max-tree owned by the caller.
     * @param adjacency Connected adjacency matching both tree domains.
     */
    DualMinMaxTreeIncrementalFilterLeaf(ValuedMorphologicalTree<T>* minTree, ValuedMorphologicalTree<T>* maxTree, const RegularGridAdjacency2D& adjacency)
        : minTree_(minTree), maxTree_(maxTree), adjacency_(&adjacency),
          collectedNodeMarks_(static_cast<std::size_t>(std::max(minTree->topology().numInternalNodeSlots(), maxTree->topology().numInternalNodeSlots()))),
          minTreeMinPixels_(initializeComponentMinPixels(minTree->topology())), maxTreeMinPixels_(initializeComponentMinPixels(maxTree->topology())) {}

    /**
     * @brief Prunes one certified max-tree leaf without rediscovering its workset.
     *
     * The support and certificates describe the current leaf and min-tree.
     * They are consumed synchronously before the primal leaf is pruned.
     *
     * @param primalLeafNode Live attached non-root leaf with a nonempty proper part.
     * @param selectedSupport Nonempty selected leaf support supplied by the builder.
     * @param dualProperPartNode Unique live dual node whose proper part contains the support.
     * @param wholeDualProperPartNode Dual node whose entire proper part is selected, or InvalidNode.
     */
    void pruneMaxLeafAndUpdateMinTree(NodeId primalLeafNode, std::span<const PixelId> selectedSupport, NodeId dualProperPartNode,
                                      NodeId wholeDualProperPartNode) {
        updateDualTree(minTree_, primalLeafNode, selectedSupport, dualProperPartNode, wholeDualProperPartNode);
        prunePrimalLeaf(*maxTree_, primalLeafNode);
    }

    /**
     * @brief Prunes one certified min-tree leaf without rediscovering its workset.
     *
     * @param primalLeafNode Live attached non-root leaf with a nonempty proper part.
     * @param selectedSupport Nonempty selected leaf support supplied by the builder.
     * @param dualProperPartNode Unique live dual node whose proper part contains the support.
     * @param wholeDualProperPartNode Dual node whose entire proper part is selected, or InvalidNode.
     */
    void pruneMinLeafAndUpdateMaxTree(NodeId primalLeafNode, std::span<const PixelId> selectedSupport, NodeId dualProperPartNode,
                                      NodeId wholeDualProperPartNode) {
        updateDualTree(maxTree_, primalLeafNode, selectedSupport, dualProperPartNode, wholeDualProperPartNode);
        prunePrimalLeaf(*minTree_, primalLeafNode);
    }

    /**
     * @brief Reads the minimum pixel in a current component, including descendants.
     * @param tree One of the two component trees borrowed by this updater.
     * @param node Live node identifier in that tree.
     * @return Minimum row-major pixel; for a leaf, also the minimum of its proper part.
     */
    [[nodiscard]] PixelId minPixel(const ValuedMorphologicalTree<T>& tree, NodeId node) const noexcept {
        return (&tree == maxTree_ ? maxTreeMinPixels_ : minTreeMinPixels_)[node];
    }

    /**
     * @brief Returns dual nodes whose leaf status or selection key may have changed.
     * Each consumed merge level contributes its union node. The buffer is cleared at the next update.
     * @return Borrowed candidate nodes, valid until the next update.
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
     * @brief Uses dense buckets only when the integral altitude domain fits the configured bit limit.
     * @return Whether the altitude type supports dense bucket storage.
     */
    static constexpr bool usesDenseBuckets() {
        if constexpr (std::is_integral_v<T> && !std::is_same_v<std::remove_cv_t<T>, bool>) {
            return std::numeric_limits<std::make_unsigned_t<T>>::digits <= MMCFILTERS_COMPONENT_TREE_ADJUSTMENT_DENSE_MAX_BITS;
        } else {
            return false;
        }
    }

    /**
     * @brief Groups affected dual nodes by altitude and records the frontier to reattach.
     * Ancestor traversal deduplicates nodes before insertion. The sweep visits active buckets in tree order;
     * frontier nodes lie outside the merge interval and attach at the primal parent altitude.
     */
    class DualMergeRegion {
        /**
         * @brief Selects dense or sparse bucket storage for the configured altitude type.
         */
        template <bool dense, typename Altitude> struct BucketStorage;

        /**
         * @brief Dense storage: one bucket for each representable altitude value.
         */
        template <typename Altitude> struct BucketStorage<true, Altitude> {
            /** @brief Number of representable altitude values. */
            static constexpr std::size_t numAltitudes = static_cast<std::size_t>(static_cast<long long>(std::numeric_limits<Altitude>::max()) -
                                                                                 static_cast<long long>(std::numeric_limits<Altitude>::lowest()) + 1);
            /** @brief Dense merge buckets indexed by altitude. */
            using type = std::array<std::vector<NodeId>, numAltitudes>;
        };

        /**
         * @brief Sparse storage used when the altitude domain is too large to index directly.
         */
        template <typename Altitude> struct BucketStorage<false, Altitude> {
            /** @brief Sparse merge buckets keyed by altitude. */
            using type = std::map<Altitude, std::vector<NodeId>>;
        };

        /** @brief Merge buckets indexed by altitude. */
        typename BucketStorage<usesDenseBuckets(), T>::type nodesByAltitude_;

        /** @brief Altitudes with nonempty merge buckets in this update. */
        std::vector<T> activeAltitudes_;

        /** @brief Boundary nodes to reattach at the primal parent altitude. */
        std::vector<NodeId> frontierNodes_;

        /** @brief Cursor into activeAltitudes_ during the sweep. */
        int currentAltitudeIndex_ = 0;

        /** @brief Direction flag: max-tree sweeps high-to-low, min-tree sweeps low-to-high. */
        bool isMaxTree_ = false;

        /**
         * @brief Maps an altitude to a dense index, offsetting signed domains by their lowest value.
         * @param altitude Representable altitude value.
         * @return Index of the corresponding dense bucket.
         */
        static std::size_t denseBucketIndex(T altitude) {
            if constexpr (std::is_signed_v<T>) {
                return static_cast<std::size_t>(static_cast<long long>(altitude) - static_cast<long long>(std::numeric_limits<T>::lowest()));
            } else {
                return static_cast<std::size_t>(altitude);
            }
        }

        /**
         * @brief Recovers the altitude represented by a dense bucket index.
         * @param index Valid index in the dense bucket storage.
         * @return Altitude represented by the bucket.
         */
        static T altitudeFromDenseBucketIndex(std::size_t index) {
            if constexpr (std::is_signed_v<T>) {
                return static_cast<T>(static_cast<long long>(index) + static_cast<long long>(std::numeric_limits<T>::lowest()));
            } else {
                return static_cast<T>(index);
            }
        }

      public:
        /**
         * @brief Clears the previous active buckets and sets the sweep direction.
         * @param dualIsMaxTree Whether the dual tree requires a descending altitude sweep.
         */
        void reset(bool dualIsMaxTree) {
            isMaxTree_ = dualIsMaxTree;
            for (T altitude : activeAltitudes_) {
                if constexpr (usesDenseBuckets()) {
                    nodesByAltitude_[denseBucketIndex(altitude)].clear();
                } else {
                    // Active levels are keys retained by the sparse storage.
                    nodesByAltitude_.find(altitude)->second.clear();
                }
            }
            activeAltitudes_.clear();
            frontierNodes_.clear();
            currentAltitudeIndex_ = 0;
        }

        /**
         * @brief Returns the bucket at the given altitude, creating it for sparse storage.
         * @param altitude Altitude of the merge bucket.
         * @return Mutable node list for that altitude.
         */
        std::vector<NodeId>& nodesAtAltitude(const T& altitude) {
            if constexpr (usesDenseBuckets()) {
                return nodesByAltitude_[denseBucketIndex(altitude)];
            } else {
                return nodesByAltitude_[altitude];
            }
        }

        /**
         * @brief Returns the frontier nodes to reattach at the primal parent altitude.
         * @return Borrowed frontier nodes for the current update.
         */
        const std::vector<NodeId>& frontierNodes() const { return frontierNodes_; }

        /**
         * @brief Records a frontier node already deduplicated by the ancestor traversal.
         * @param node Dual node outside the merge interval.
         */
        void addFrontierNode(NodeId node) { frontierNodes_.push_back(node); }

        /**
         * @brief Adds a deduplicated dual node to the bucket at its current altitude.
         * @param tree Dual tree containing the node.
         * @param node Live dual node in the merge interval.
         */
        void addMergeNode(const ValuedMorphologicalTree<T>& tree, NodeId node) {
            auto& bucket = nodesAtAltitude(tree.nodeAltitudes()[static_cast<std::size_t>(node)]);
            bucket.push_back(node);
        }

        /**
         * @brief Collects active levels and positions the sweep cursor.
         */
        void prepareSweep() {
            activeAltitudes_.clear();
            if constexpr (usesDenseBuckets()) {
                for (std::size_t bucketIndex = 0; bucketIndex < nodesByAltitude_.size(); ++bucketIndex) {
                    if (!nodesByAltitude_[bucketIndex].empty()) {
                        activeAltitudes_.push_back(altitudeFromDenseBucketIndex(bucketIndex));
                    }
                }
            } else {
                activeAltitudes_.reserve(nodesByAltitude_.size());
                for (const auto& entry : nodesByAltitude_) {
                    if (!entry.second.empty()) {
                        activeAltitudes_.push_back(entry.first);
                    }
                }
            }

            currentAltitudeIndex_ = isMaxTree_ ? static_cast<int>(activeAltitudes_.size()) - 1 : 0;
        }

        /**
         * @brief Tests whether the sweep cursor identifies an active altitude.
         * @return Whether the current sweep position is valid.
         */
        bool hasCurrentAltitude() const { return currentAltitudeIndex_ >= 0 && static_cast<std::size_t>(currentAltitudeIndex_) < activeAltitudes_.size(); }

        /**
         * @brief Reads the altitude at the current sweep position.
         * @return Altitude of the current active bucket.
         */
        T currentAltitude() const { return activeAltitudes_[static_cast<std::size_t>(currentAltitudeIndex_)]; }

        /** @brief Advances the sweep cursor according to the tree polarity. */
        void advance() { currentAltitudeIndex_ += isMaxTree_ ? -1 : 1; }
    };

    /** @brief Mutable min-tree owned by the caller. */
    ValuedMorphologicalTree<T>* minTree_ = nullptr;

    /** @brief Mutable max-tree owned by the caller. */
    ValuedMorphologicalTree<T>* maxTree_ = nullptr;

    /** @brief Shared pixel adjacency used to find dual neighbors of the removed leaf support. */
    const RegularGridAdjacency2D* adjacency_ = nullptr;

    /** @brief Merge buckets and frontier nodes for the current dual update. */
    DualMergeRegion dualMergeRegion_;

    /** @brief Deduplicates adjacency seeds and ancestors before collection insertion. */
    GenerationStampSet collectedNodeMarks_;

    /** @brief Minimum row-major pixel in each current min-tree component, including descendants. */
    std::vector<PixelId> minTreeMinPixels_;

    /** @brief Minimum row-major pixel in each current max-tree component, including descendants. */
    std::vector<PixelId> maxTreeMinPixels_;

    /** @brief Dual nodes whose leaf status or selection key may have changed. */
    std::vector<NodeId> lastCandidateNodes_;

    /** Complete proper-part transfers needed to update the pending support indices. */
    std::vector<ProperPartTransfer> lastProperPartTransfers_;

    /** Dual proper part containing the complete selected support after the last update. */
    NodeId lastSelectedSupportNode_ = InvalidNode;

    /**
     * @brief Initializes component minima in row-major order, assigning each component once.
     * @param topology Component tree over the nonempty row-major image domain.
     * @return Minimum pixel for every live component, indexed by node slot.
     */
    static std::vector<PixelId> initializeComponentMinPixels(const MorphologicalTree& topology) {
        std::vector<PixelId> minPixels(topology.numInternalNodeSlots(), InvalidPixel);
        minPixels[topology.root()] = 0;
        if (topology.numNodes() == 1)
            return minPixels;
        for (PixelId pixel = 0; pixel < topology.numPixels(); ++pixel) {
            NodeId node = TreeKernelAccess::smallestNodeMap(topology, pixel);
            // The first pixel reaching a component is its minimum. A previously
            // assigned component also has previously assigned ancestors.
            while (minPixels[node] == InvalidPixel) {
                minPixels[node] = pixel;
                node = TreeKernelAccess::nodeParent(topology, node);
            }
        }
        return minPixels;
    }

    /**
     * @brief Reads the altitude of an established live node.
     * @param tree Component tree containing the node.
     * @param node Live node identifier in that tree.
     * @return Stored node altitude.
     */
    static T nodeAltitude(const ValuedMorphologicalTree<T>* tree, NodeId node) noexcept { return tree->nodeAltitudes()[static_cast<std::size_t>(node)]; }

    /**
     * @brief Reads the topology of a borrowed component tree.
     * @param tree Borrowed component tree.
     * @return Topology retained by the tree.
     */
    static const MorphologicalTree& topologyOf(const ValuedMorphologicalTree<T>* tree) { return tree->topology(); }

    /**
     * @brief Selects the primal tree from the polarity of the dual tree.
     * @param dualIsMaxTree Whether the dual tree is the max-tree.
     * @return Mutable component tree with the opposite polarity.
     */
    ValuedMorphologicalTree<T>* primalTreeFor(bool dualIsMaxTree) { return dualIsMaxTree ? minTree_ : maxTree_; }

    /**
     * @brief Clears the previous update results and scratch state.
     * @param dualIsMaxTree Whether the next update modifies the max-tree.
     */
    void resetUpdateState(bool dualIsMaxTree) {
        dualMergeRegion_.reset(dualIsMaxTree);
        collectedNodeMarks_.resetAll();
        lastCandidateNodes_.clear();
        lastProperPartTransfers_.clear();
        lastSelectedSupportNode_ = InvalidNode;
    }

    /**
     * @brief Releases an attached non-root node after its contents have been transferred.
     * @param topology Current topology with established node slots.
     * @param editor Active tree editor.
     * @param node Live, attached, empty non-root node.
     */
    static void releaseEmptyNode(const MorphologicalTree& topology, ValuedMorphologicalTreeEditor<T>& editor, NodeId node) {
        const NodeId parentNode = TreeKernelAccess::nodeParent(topology, node);
        assert(node != topology.root() && parentNode != InvalidNode && parentNode != node);
        editor.removeChild(parentNode, node, true);
    }

    /**
     * @brief Moves the children of sourceNode under targetNode, leaving targetNode itself in place.
     * @param topology Current topology with established node slots.
     * @param editor Active tree editor.
     * @param sourceNode Node whose child list is traversed.
     * @param targetNode Live destination, which may itself be a child of sourceNode.
     */
    static void reattachChildren(const MorphologicalTree& topology, ValuedMorphologicalTreeEditor<T>& editor, NodeId sourceNode, NodeId targetNode) {
        for (NodeId childNode = TreeKernelAccess::firstChild(topology, sourceNode); childNode != InvalidNode;) {
            const NodeId nextChildNode = TreeKernelAccess::nextSibling(topology, childNode);
            if (childNode != targetNode) {
                editor.detach(childNode);
                editor.attach(targetNode, childNode);
            }
            childNode = nextChildNode;
        }
    }

    /**
     * @brief Moves every child and the complete proper part from a donor to the union node.
     *
     * The donor remains attached and empty until the caller releases it.
     *
     * @param editor Active tree editor.
     * @param unionNode Live union node at the current merge level.
     * @param donorNode Live donor from the same bucket, with a nonempty proper part.
     */
    void mergeDonorIntoUnion(ValuedMorphologicalTreeEditor<T>& editor, NodeId unionNode, NodeId donorNode) {
        editor.moveChildren(unionNode, donorNode);
        lastProperPartTransfers_.push_back({donorNode, unionNode});
        editor.mergeProperParts(unionNode, donorNode);
    }

    /**
     * @brief Moves the selected support from its established dual proper part to the union node.
     *
     * The sweep excludes the altitude of dualProperPartNode, so its pixels still belong to
     * dualProperPartNode when the target level is reached. Only this proper part can become
     * empty. Its contraction is handled when the updated branch is reattached.
     *
     * @param dualTree Dual tree.
     * @param editor Active tree editor.
     * @param unionNode Union node at the primal parent altitude.
     * @param support Nonempty selected support contained in the proper part of dualProperPartNode.
     * @param dualProperPartNode Established dual proper-part node containing the support.
     * @param wholeDualProperPartNode dualProperPartNode when its complete proper part is selected, or InvalidNode.
     * @return Whether the dual proper part became empty.
     */
    bool moveSelectedSupport(ValuedMorphologicalTree<T>* dualTree, ValuedMorphologicalTreeEditor<T>& editor, NodeId unionNode, std::span<const PixelId> support,
                             NodeId dualProperPartNode, NodeId wholeDualProperPartNode) {
        lastSelectedSupportNode_ = unionNode;
        assert(unionNode != dualProperPartNode);
        if (wholeDualProperPartNode != InvalidNode) {
            assert(wholeDualProperPartNode == dualProperPartNode);
            lastProperPartTransfers_.push_back({dualProperPartNode, unionNode});
            editor.mergeProperParts(unionNode, dualProperPartNode);
            return true;
        }
        for (PixelId pixelId : support) {
            editor.movePixelToProperPart(unionNode, dualProperPartNode, pixelId);
        }
        return TreeKernelAccess::properPartCardinality(topologyOf(dualTree), dualProperPartNode) == 0;
    }

    /**
     * @brief Reattaches the detached updated branch and contracts its empty dual proper-part node.
     * @param dualTree Dual tree after the merge-level sweep.
     * @param editor Active tree editor.
     * @param dualProperPartNode Original dual proper-part node containing the selected support.
     * @param finalUnionNode Live union node left detached by the last active merge level.
     * @param dualProperPartEmptied Whether moving the support emptied dualProperPartNode.
     */
    void reattachMergedBranch(ValuedMorphologicalTree<T>* dualTree, ValuedMorphologicalTreeEditor<T>& editor, NodeId dualProperPartNode, NodeId finalUnionNode,
                              bool dualProperPartEmptied) {
        const MorphologicalTree& topology = topologyOf(dualTree);
        auto& componentMinPixels = dualTree == maxTree_ ? maxTreeMinPixels_ : minTreeMinPixels_;
        assert(TreeKernelAccess::nodeParent(topology, finalUnionNode) == finalUnionNode);
        if (!dualProperPartEmptied) {
            editor.attach(dualProperPartNode, finalUnionNode);
            return;
        }
        if (dualProperPartNode != topology.root()) {
            editor.attach(TreeKernelAccess::nodeParent(topology, dualProperPartNode), finalUnionNode);
            reattachChildren(topology, editor, dualProperPartNode, finalUnionNode);
            componentMinPixels[finalUnionNode] = std::min(componentMinPixels[finalUnionNode], componentMinPixels[dualProperPartNode]);
            releaseEmptyNode(topology, editor, dualProperPartNode);
            return;
        }

        const bool dualIsMaxTree = dualTree == maxTree_;
        NodeId newRoot = finalUnionNode;
        for (NodeId childNode = TreeKernelAccess::firstChild(topology, dualProperPartNode); childNode != InvalidNode;
             childNode = TreeKernelAccess::nextSibling(topology, childNode)) {
            if (dualIsMaxTree ? nodeAltitude(dualTree, childNode) < nodeAltitude(dualTree, newRoot)
                              : nodeAltitude(dualTree, childNode) > nodeAltitude(dualTree, newRoot)) {
                newRoot = childNode;
            }
        }
        if (newRoot != finalUnionNode) {
            editor.attach(newRoot, finalUnionNode);
        }
        reattachChildren(topology, editor, dualProperPartNode, newRoot);
        componentMinPixels[newRoot] = componentMinPixels[dualProperPartNode];
        editor.setRoot(newRoot);
        editor.releaseNode(dualProperPartNode);
    }

    /**
     * @brief Contracts an established primal leaf through the existing edit primitives.
     * The component of its parent is unchanged, so its stored minimum remains valid.
     * @param tree Primal tree containing an attached non-root leaf.
     * @param primalLeafNode Live leaf whose support is transferred to its parent.
     */
    static void prunePrimalLeaf(ValuedMorphologicalTree<T>& tree, NodeId primalLeafNode) {
        const NodeId parent = TreeKernelAccess::nodeParent(tree.topology(), primalLeafNode);
        auto editor = mmcfilters::detail::beginEstablishedValuedEdit(tree);
        editor.mergeProperParts(parent, primalLeafNode);
        editor.removeChild(parent, primalLeafNode, true);
        auto proof = editor.proveIncremental();
        editor.commit(std::move(proof));
    }

    /**
     * @brief Collects dual merge buckets and frontier nodes from the selected support boundary.
     * Neighbors in dualProperPartNode are skipped because its altitude is excluded from the merge sweep.
     * Other neighbors seed ancestor paths. A generation mark deduplicates seeds and ancestors before insertion.
     * The dual tree remains unchanged during collection.
     *
     * @param dualTree Dual tree before the update.
     * @param support Nonempty proper part of the selected primal leaf.
     * @param dualProperPartNode Unique dual proper part containing the support.
     * @param targetAltitude Primal parent altitude reached by this leveling.
     * @param dualIsMaxTree Whether the dual tree is a max-tree.
     */
    void collectDualMergeRegion(ValuedMorphologicalTree<T>* dualTree, std::span<const PixelId> support, NodeId dualProperPartNode, T targetAltitude,
                                bool dualIsMaxTree) {
        const MorphologicalTree& dualTopology = topologyOf(dualTree);
        const T dualProperPartAltitude = nodeAltitude(dualTree, dualProperPartNode);

        const auto collectAncestorPath = [&](NodeId seedNode) {
            NodeId ancestorNode = seedNode;
            while (!collectedNodeMarks_.isMarked(static_cast<std::size_t>(ancestorNode))) {
                const T currentAltitude = nodeAltitude(dualTree, ancestorNode);
                if (!((dualIsMaxTree && currentAltitude >= dualProperPartAltitude) || (!dualIsMaxTree && currentAltitude <= dualProperPartAltitude))) {
                    break;
                }

                collectedNodeMarks_.mark(static_cast<std::size_t>(ancestorNode));

                const NodeId parentNode = TreeKernelAccess::nodeParent(dualTopology, ancestorNode);
                if (dualIsMaxTree ? currentAltitude <= targetAltitude : currentAltitude >= targetAltitude) {
                    dualMergeRegion_.addMergeNode(*dualTree, ancestorNode);
                } else {
                    if (parentNode == ancestorNode ||
                        !(dualIsMaxTree ? nodeAltitude(dualTree, parentNode) > targetAltitude : nodeAltitude(dualTree, parentNode) < targetAltitude)) {
                        dualMergeRegion_.addFrontierNode(ancestorNode);
                    }
                }

                if (parentNode == ancestorNode) {
                    break;
                }
                ancestorNode = parentNode;
            }
        };

        for (PixelId pixel : support) {
            for (PixelId neighbor : mmcfilters::detail::CommittedGridAccess::neighbors(*adjacency_, pixel)) {
                const NodeId adjacentNode = TreeKernelAccess::smallestNodeMap(dualTopology, neighbor);
                if (adjacentNode != dualProperPartNode) {
                    collectAncestorPath(adjacentNode);
                }
            }
        }
    }

    /**
     * @brief Merges the affected dual region and moves the selected support to the primal parent altitude.
     * At each altitude, the smallest node identifier represents the union. Frontiers attach at the target altitude;
     * each later union receives the previous one. The final branch is reattached before publishing the edit.
     * The caller contracts the primal leaf after this update.
     *
     * @param dualTree Dual tree to update.
     * @param primalLeafNode Live attached non-root primal leaf.
     * @param support Nonempty selected leaf support supplied by the builder.
     * @param dualProperPartNode Unique live dual proper part containing the support.
     * @param wholeDualProperPartNode dualProperPartNode when its complete proper part is selected, or InvalidNode.
     */
    void updateDualTree(ValuedMorphologicalTree<T>* dualTree, NodeId primalLeafNode, std::span<const PixelId> support, NodeId dualProperPartNode,
                        NodeId wholeDualProperPartNode) {
        const bool dualIsMaxTree = dualTree == maxTree_;
        resetUpdateState(dualIsMaxTree);
        ValuedMorphologicalTree<T>* primalTree = primalTreeFor(dualIsMaxTree);
        const MorphologicalTree& dualTopology = topologyOf(dualTree);
        const MorphologicalTree& primalTopology = topologyOf(primalTree);
        auto& componentMinPixels = dualIsMaxTree ? maxTreeMinPixels_ : minTreeMinPixels_;

        const NodeId primalParentNode = TreeKernelAccess::nodeParent(primalTopology, primalLeafNode);
        const T targetAltitude = nodeAltitude(primalTree, primalParentNode);
        const T dualProperPartAltitude = nodeAltitude(dualTree, dualProperPartNode);
        lastSelectedSupportNode_ = dualProperPartNode;
        collectDualMergeRegion(dualTree, support, dualProperPartNode, targetAltitude, dualIsMaxTree);

        ValuedMorphologicalTreeEditor<T> editor = mmcfilters::detail::beginEstablishedValuedEdit(*dualTree);
        NodeId previousUnionNode = InvalidNode;
        bool dualProperPartEmptied = false;

        for (dualMergeRegion_.prepareSweep(); dualMergeRegion_.hasCurrentAltitude(); dualMergeRegion_.advance()) {
            const T mergeAltitude = dualMergeRegion_.currentAltitude();
            if (!(dualIsMaxTree ? mergeAltitude > dualProperPartAltitude : mergeAltitude < dualProperPartAltitude)) {
                break;
            }
            const auto& nodesAtMergeAltitude = dualMergeRegion_.nodesAtAltitude(mergeAltitude);
            // Active buckets are nonempty; each collected node occurs in exactly one bucket.
            const NodeId unionNode = *std::min_element(nodesAtMergeAltitude.begin(), nodesAtMergeAltitude.end());
            PixelId unionMinPixel = componentMinPixels[unionNode];
            for (NodeId mergeNode : nodesAtMergeAltitude) {
                assert(TreeKernelAccess::isAlive(dualTopology, mergeNode));
                if (mergeNode == unionNode) {
                    // The strict altitude interval excludes the dual root.
                    editor.detach(unionNode);
                } else {
                    assert(TreeKernelAccess::properPartCardinality(dualTopology, mergeNode) > 0 && lastSelectedSupportNode_ != mergeNode);
                    mergeDonorIntoUnion(editor, unionNode, mergeNode);
                    unionMinPixel = std::min(unionMinPixel, componentMinPixels[mergeNode]);
                    releaseEmptyNode(dualTopology, editor, mergeNode);
                }
            }
            lastCandidateNodes_.push_back(unionNode);

            if (mergeAltitude == targetAltitude) {
                dualProperPartEmptied = moveSelectedSupport(dualTree, editor, unionNode, support, dualProperPartNode, wholeDualProperPartNode);
                unionMinPixel = std::min(unionMinPixel, minPixel(*primalTree, primalLeafNode));
                for (NodeId frontierNode : dualMergeRegion_.frontierNodes()) {
                    // Frontier altitudes lie outside the merge interval.
                    assert(TreeKernelAccess::isAlive(dualTopology, frontierNode) && frontierNode != dualTopology.root());
                    editor.detach(frontierNode);
                    editor.attach(unionNode, frontierNode);
                    unionMinPixel = std::min(unionMinPixel, componentMinPixels[frontierNode]);
                }
            }

            // Each union stays detached until it is attached below the next level.
            if (previousUnionNode != InvalidNode) {
                editor.attach(unionNode, previousUnionNode);
                unionMinPixel = std::min(unionMinPixel, componentMinPixels[previousUnionNode]);
            }
            componentMinPixels[unionNode] = unionMinPixel;
            previousUnionNode = unionNode;
        }

        if (previousUnionNode != InvalidNode) {
            reattachMergedBranch(dualTree, editor, dualProperPartNode, previousUnionNode, dualProperPartEmptied);
        }
        auto proof = editor.proveIncremental();
        editor.commit(std::move(proof));
    }
};

} // namespace mmcfilters::detail
