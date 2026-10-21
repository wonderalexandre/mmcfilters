#pragma once

#include "ComponentTreeBuilder.hpp"
#include "NativeHierarchyOwned.hpp"
#include "ValuedMorphologicalTree.hpp"
#include "detail/DualMinMaxTreeIncrementalFilterLeaf.hpp"
#include "detail/TreeKernelAccess.hpp"
#include "../utils/Common.hpp"

#include <cassert>
#include <cstddef>
#include <iterator>
#include <list>
#include <optional>
#include <queue>
#include <span>
#include <utility>
#include <vector>

namespace mmcfilters {
namespace detail {
/**
 * @brief Test access to support preparation and pending support invariants.
 * @tparam T Image altitude type.
 */
template <AltitudeValue T> struct ExtremumResidualTreeTestAccess;
} // namespace detail

/**
 * Algorithm 1 of paper: "The Extremum Residual Tree: A Self-Dual Morphological Tree Based on Levelings", DGMM 2027
 *
 * ExtremumAgenda orders nonterminal leaves by support cardinality and minimum row-major pixel (Eqs. 16, 30).
 * ResidualTreeConstruction records supports, maintains Xi, and indexes pending supports by their containing proper parts in both component trees
 * (Eqs. 31-32). Selected supports come from current leaf proper parts. DualMinMaxTreeIncrementalFilterLeaf updates the dual tree after each leveling.
 * The min-tree and max-tree are temporary construction state. Each build returns an independent owning native hierarchy.
 * Callers provide a non-null, nonempty image with finite values and a connected adjacency matching its domain.
 *
 * @tparam T Image and node altitude type.
 */
template <AltitudeValue T> class ExtremumResidualTreeByCTsBuilder {
  public:
    /**
     * @brief Configures component-tree construction with a shared adjacency and row-major spatial ties.
     * @param adjacency Connected spatial adjacency matching the input image domain.
     */
    explicit ExtremumResidualTreeByCTsBuilder(RegularGridAdjacency2D adjacency) : adjacency_(std::move(adjacency)) {}

    /**
     * @brief Builds temporary min/max trees and records the extremum residual hierarchy.
     * @param image Non-null, nonempty image with finite values and the same domain as the connected adjacency.
     * @return Independent native hierarchy owning its buffers and adjacency metadata.
     */
    [[nodiscard]] NativeHierarchyOwned<T> build(const ImagePtr<T>& image) {
        auto minTree = detail::materializeNativeHierarchy<T>(ComponentTreeBuilder(adjacency_, false).build(image), std::nullopt,
                                                             detail::MaterializedVersionPolicy::PreserveLinkedConstruction);
        auto maxTree = detail::materializeNativeHierarchy<T>(ComponentTreeBuilder(adjacency_, true).build(image), std::nullopt,
                                                             detail::MaterializedVersionPolicy::PreserveLinkedConstruction);
        const auto maxTreeNodeSlots = maxTree.topology().numInternalNodeSlots();
        const auto minTreeNodeSlots = minTree.topology().numInternalNodeSlots();
        detail::DualMinMaxTreeIncrementalFilterLeaf<T> updater(&minTree, &maxTree, adjacency_);
        ExtremumAgenda agenda;
        SelectedSupport selectedSupport;
        ResidualTreeConstruction extTree(image->getSize(), maxTreeNodeSlots, minTreeNodeSlots);

        agenda.initialize(maxTree, minTree, updater);                 // Line 2.
        while (const auto selected = agenda.next(maxTree, minTree)) { // Lines 3-5.
            const auto& selectedExtremum = *selected;
            const bool isMaximum = selectedExtremum.polarity == Polarity::Maximum;
            ValuedMorphologicalTree<T>& primal = isMaximum ? maxTree : minTree;
            ValuedMorphologicalTree<T>& dual = isMaximum ? minTree : maxTree;
            selectedSupport.collectFromLeaf(selectedExtremum.nodeId, primal, dual);
            const NodeId primalParentNode = detail::TreeKernelAccess::nodeParent(primal.topology(), selectedExtremum.nodeId);
            const T nodeAltitude = primal.nodeAltitudes()[static_cast<std::size_t>(selectedExtremum.nodeId)];

            // Lines 6-11: add X_s and attach the pending supports contained in X_s as children.
            const NodeId residualNode = extTree.addSelectedSupport(selectedSupport, selectedExtremum, nodeAltitude);

            // Line 12: apply the elementary leveling.
            applyElementaryLeveling(updater, isMaximum, selectedExtremum.nodeId, selectedSupport);

            // Line 13: pending supports follow the proper-part transfers.
            extTree.updatePendingSupports(residualNode, selectedExtremum.polarity, primalParentNode, updater.getLastSelectedSupportNode(),
                                          updater.getLastProperPartTransfers());

            const Polarity dualPolarity = oppositePolarity(selectedExtremum.polarity);
            for (NodeId node : updater.getLastCandidateNodes())
                agenda.enqueue(dual, dualPolarity, node, updater);
            agenda.enqueue(primal, selectedExtremum.polarity, primalParentNode, updater);
        }

        // Lines 16-18: attach the remaining pending supports to the terminal root.
        const T terminalAltitude = minTree.nodeAltitudes()[static_cast<std::size_t>(minTree.topology().root())];
        assert(minTree.topology().numNodes() == 1 && maxTree.topology().numNodes() == 1);
        assert(maxTree.nodeAltitudes()[static_cast<std::size_t>(maxTree.topology().root())] == terminalAltitude);
        return std::move(extTree).completeTree(
            image, makeMorphologicalTreeSemantics(MorphologicalTreeKind::ExtremumResidualTree, SharedAdjacencyContext{adjacency_}), terminalAltitude);
    }

    /**
     * @brief Exposes the shared component-tree adjacency.
     * @return Borrowed adjacency retained by this builder.
     */
    [[nodiscard]] const RegularGridAdjacency2D& adjacency() const noexcept { return adjacency_; }

  private:
    friend struct detail::ExtremumResidualTreeTestAccess<T>;
    /**
     * @brief Component-tree side containing the selected regional extremum.
     */
    enum class Polarity {
        Maximum, ///< The selected max-tree leaf is lowered while the min-tree is updated.
        Minimum  ///< The selected min-tree leaf is raised while the max-tree is updated.
    };
    /**
     * @brief Returns the opposite polarity.
     * @param polarity Selected component-tree side.
     * @return Opposite component-tree side.
     */
    [[nodiscard]] static constexpr Polarity oppositePolarity(Polarity polarity) {
        return polarity == Polarity::Maximum ? Polarity::Minimum : Polarity::Maximum;
    }

    /**
     * @brief Current extremal flat zone and its deterministic selection key.
     */
    struct ExtremumCandidate {
        std::size_t supportCardinality = 0;  ///< Cardinality of the current extremum support.
        PixelId minimumPixel = InvalidPixel; ///< Smallest row-major pixel identifier in the support.
        Polarity polarity;                   ///< Component-tree side in which the flat zone is a leaf.
        NodeId nodeId;                       ///< Leaf identifier in that component tree.
    };
    /**
     * @brief Gives priority to smaller supports and then smaller spatial minima.
     */
    struct ExtremumCandidatePriority {
        /**
         * @brief Compares candidate priorities for a minimum-first queue.
         * @param lhs First candidate.
         * @param rhs Second candidate.
         * @return Whether lhs has lower priority than rhs.
         */
        bool operator()(const ExtremumCandidate& lhs, const ExtremumCandidate& rhs) const {
            if (lhs.supportCardinality != rhs.supportCardinality)
                return lhs.supportCardinality > rhs.supportCardinality;
            return lhs.minimumPixel > rhs.minimumPixel;
        }
    };
    /**
     * @brief Orders extremum candidates and discards outdated entries when selecting the next leaf.
     * Every changed leaf is enqueued after a leveling; queue entries retain immutable selection keys.
     */
    struct ExtremumAgenda {
        std::priority_queue<ExtremumCandidate, std::vector<ExtremumCandidate>, ExtremumCandidatePriority> candidates; ///< Minimum-first candidate queue.

        /**
         * @brief Queues the current selection key when a node is a nonterminal leaf.
         * @param tree Component tree containing the candidate.
         * @param polarity Component-tree side containing the candidate.
         * @param node Allocated node slot whose leaf status or support may have changed.
         * @param updater Maintained component minima for both trees.
         */
        void enqueue(const ValuedMorphologicalTree<T>& tree, Polarity polarity, NodeId node, const detail::DualMinMaxTreeIncrementalFilterLeaf<T>& updater) {
            if (!isNonterminalLeaf(tree, node))
                return;
            candidates.push({static_cast<std::size_t>(detail::TreeKernelAccess::properPartCardinality(tree.topology(), node)), updater.minPixel(tree, node),
                             polarity, node});
        }

        /**
         * @brief Queues all initial nonterminal leaves.
         * @param maxTree Initial max-tree.
         * @param minTree Initial min-tree.
         * @param updater Initialized component minima for both trees.
         */
        void initialize(const ValuedMorphologicalTree<T>& maxTree, const ValuedMorphologicalTree<T>& minTree,
                        const detail::DualMinMaxTreeIncrementalFilterLeaf<T>& updater) {
            for (NodeId node : maxTree.topology().aliveNodeIds())
                enqueue(maxTree, Polarity::Maximum, node, updater);
            for (NodeId node : minTree.topology().aliveNodeIds())
                enqueue(minTree, Polarity::Minimum, node, updater);
        }

        /**
         * @brief Pops the smallest eligible extremum, discarding outdated entries.
         * @param maxTree Current max-tree.
         * @param minTree Current min-tree.
         * @return Smallest current extremum, or no value at termination.
         */
        [[nodiscard]] std::optional<ExtremumCandidate> next(const ValuedMorphologicalTree<T>& maxTree, const ValuedMorphologicalTree<T>& minTree) {
            while (!candidates.empty()) {
                const ExtremumCandidate selected = candidates.top();
                candidates.pop();
                const ValuedMorphologicalTree<T>& primal = selected.polarity == Polarity::Maximum ? maxTree : minTree;
                if (!isNonterminalLeaf(primal, selected.nodeId))
                    continue;
                // Flat zones only merge. Equal cardinality and retention of the
                // recorded minimum certify that the selection key is unchanged.
                if (static_cast<std::size_t>(detail::TreeKernelAccess::properPartCardinality(primal.topology(), selected.nodeId)) !=
                    selected.supportCardinality)
                    continue;
                if (detail::TreeKernelAccess::smallestNodeMap(primal.topology(), selected.minimumPixel) != selected.nodeId)
                    continue;
                return selected;
            }
            return std::nullopt;
        }
    };

    /**
     * Selected extremal flat zone and dual proper parts reused by the leaf updater.
     */
    struct SelectedSupport {
        std::vector<PixelId> pixels;                  ///< Pixels in the selected extremal flat zone.
        NodeId dualProperPartNode = InvalidNode;      ///< Unique dual proper part containing the selected flat zone.
        NodeId wholeDualProperPartNode = InvalidNode; ///< Dual proper-part node when the support fills that part; otherwise InvalidNode.

        /**
         * @brief Collects a leaf proper part and establishes certificates for the dual update.
         * @param leaf Selected non-root primal leaf.
         * @param primal Component tree containing the selected extremum as a leaf.
         * @param dual Opposite component tree containing the support in one proper part.
         */
        void collectFromLeaf(NodeId leaf, const ValuedMorphologicalTree<T>& primal, const ValuedMorphologicalTree<T>& dual) {
            const auto& dualTopology = dual.topology();
            pixels.clear();
            wholeDualProperPartNode = InvalidNode;
            for (PixelId pixel : detail::TreeKernelAccess::properParts(primal.topology(), leaf)) {
                pixels.push_back(pixel);
            }
            // An extremal flat zone is nonempty and lies in one dual proper part.
            dualProperPartNode = detail::TreeKernelAccess::smallestNodeMap(dualTopology, pixels.front());
#ifndef NDEBUG
            for (PixelId pixel : pixels) {
                assert(detail::TreeKernelAccess::smallestNodeMap(dualTopology, pixel) == dualProperPartNode);
            }
#endif
            if (detail::TreeKernelAccess::properPartCardinality(dualTopology, dualProperPartNode) == static_cast<int>(pixels.size())) {
                wholeDualProperPartNode = dualProperPartNode;
            }
        }
    };
    /**
     * @brief Applies the elementary leveling to the selected extremum in both component trees.
     * @param updater Updater borrowing the synchronized min-tree and max-tree.
     * @param isMaximum Whether the selected primal leaf belongs to the max-tree.
     * @param selectedLeafNode Selected non-root primal leaf.
     * @param selectedSupport Selected pixels and dual proper-part identifiers.
     */
    [[gnu::always_inline]] static void applyElementaryLeveling(detail::DualMinMaxTreeIncrementalFilterLeaf<T>& updater, bool isMaximum, NodeId selectedLeafNode,
                                                               const SelectedSupport& selectedSupport) {
        if (isMaximum)
            updater.pruneMaxLeafAndUpdateMinTree(selectedLeafNode, selectedSupport.pixels, selectedSupport.dualProperPartNode,
                                                 selectedSupport.wholeDualProperPartNode);
        else
            updater.pruneMinLeafAndUpdateMaxTree(selectedLeafNode, selectedSupport.pixels, selectedSupport.dualProperPartNode,
                                                 selectedSupport.wholeDualProperPartNode);
    }

    /**
     * @brief Indexes pending residual supports by the component proper part containing them.
     * Pending support X is stored in residualNodesByProperPart_[Y] when the proper part of component node Y contains X (Eq. 31).
     */
    class PendingSupportsByProperPart {
        std::vector<std::list<NodeId>> residualNodesByProperPart_;        ///< Pending residual nodes indexed by component-tree proper part.
        std::vector<std::list<NodeId>::iterator> positionByResidualNode_; ///< Position in the pending support list, indexed by residual node.

      public:
        /**
         * @brief Allocates pending support lists for the fixed component-tree node slots.
         * @param componentNodeSlots Number of allocated component-tree node slots.
         */
        explicit PendingSupportsByProperPart(std::size_t componentNodeSlots) : residualNodesByProperPart_(componentNodeSlots) {}

        /**
         * @brief Returns pending residual supports contained in one component proper part.
         * @param properPartNode Component node owning the proper part.
         * @return Mutable list of pending residual nodes.
         */
        std::list<NodeId>& pendingSupportsIn(NodeId properPartNode) { return residualNodesByProperPart_[static_cast<std::size_t>(properPartNode)]; }
        /**
         * @brief Reads pending residual supports contained in one component proper part.
         * @param properPartNode Component node owning the proper part.
         * @return Borrowed list of pending residual nodes.
         */
        const std::list<NodeId>& pendingSupportsIn(NodeId properPartNode) const { return residualNodesByProperPart_[static_cast<std::size_t>(properPartNode)]; }

        /**
         * @brief Adds a pending residual support to the list of its containing proper part.
         * @param residualNode Recorded residual support awaiting an inclusion parent.
         * @param properPartNode Component-tree proper part containing the support.
         */
        void addPendingSupport(NodeId residualNode, NodeId properPartNode) {
            auto& pendingNodes = pendingSupportsIn(properPartNode);
            if (static_cast<std::size_t>(residualNode) >= positionByResidualNode_.size())
                positionByResidualNode_.resize(static_cast<std::size_t>(residualNode) + 1);
            pendingNodes.push_back(residualNode);
            positionByResidualNode_[static_cast<std::size_t>(residualNode)] = std::prev(pendingNodes.end());
        }

        /**
         * @brief Removes a pending support after its inclusion parent is recorded.
         * @param residualNode Residual node identifying the pending support to remove.
         * @param properPartNode Component node owning the proper part that contains the support.
         */
        void removePendingSupport(NodeId residualNode, NodeId properPartNode) {
            pendingSupportsIn(properPartNode).erase(positionByResidualNode_[static_cast<std::size_t>(residualNode)]);
        }

        /**
         * @brief Transfers pending supports with a whole proper-part transfer, preserving list iterators.
         * @param sourceNode Component node whose proper part was transferred.
         * @param targetNode Live component node receiving that proper part.
         */
        void transferPendingSupports(NodeId sourceNode, NodeId targetNode) {
            auto& targetPendingSupports = pendingSupportsIn(targetNode);
            targetPendingSupports.splice(targetPendingSupports.end(), pendingSupportsIn(sourceNode));
        }
    };

    /**
     * @brief Builds residual nodes, their inclusion parents, and their proper-part valuations.
     * Xi is implicit: non-root recorded nodes with nodeParents[residualNode] == InvalidNode.
     */
    struct ResidualTreeConstruction {
        PendingSupportsByProperPart maxTreePendingSupports; ///< Pending supports contained in max-tree proper parts.
        PendingSupportsByProperPart minTreePendingSupports; ///< Pending supports contained in min-tree proper parts.
        std::vector<NodeId> nodeParents{InvalidNode};       ///< Inclusion parents; slot zero is reserved for the terminal root.
        std::vector<NodeId> smallestNodeMap;                ///< First recorded support containing each pixel, or InvalidNode until assigned.
        std::vector<T> nodeAltitudes{T{}};                  ///< Proper-part valuations, with the terminal valuation reserved at slot zero.
        /**
         * @brief Allocates residual node buffers and pending support indices for both component trees.
         * @param numPixels Positive representable pixel count and upper bound on residual nodes.
         * @param maxTreeNodeSlots Initial number of max-tree node slots.
         * @param minTreeNodeSlots Initial number of min-tree node slots.
         */
        ResidualTreeConstruction(std::size_t numPixels, std::size_t maxTreeNodeSlots, std::size_t minTreeNodeSlots)
            : maxTreePendingSupports(maxTreeNodeSlots), minTreePendingSupports(minTreeNodeSlots), smallestNodeMap(numPixels, InvalidNode) {
            nodeParents.reserve(numPixels);
            nodeAltitudes.reserve(numPixels);
        }

        /**
         * @brief Adds the selected support as a residual node and attaches its pending children.
         * @param selectedSupport Selected pixels and dual proper-part identifiers.
         * @param selectedExtremum Selected leaf and component-tree polarity.
         * @param supportAltitude Valuation of the selected extremum before leveling.
         * @return Identifier of the newly recorded residual node.
         */
        NodeId addSelectedSupport(const SelectedSupport& selectedSupport, const ExtremumCandidate& selectedExtremum, T supportAltitude) {
            // Property 6 gives each recorded support a new nonempty proper part;
            // together with the terminal root this bounds the node count by pixels.
            PendingSupportsByProperPart& primalPendingSupports = pendingSupportsFor(selectedExtremum.polarity);
            PendingSupportsByProperPart& dualPendingSupports = pendingSupportsFor(oppositePolarity(selectedExtremum.polarity));
            auto& pendingChildren = primalPendingSupports.pendingSupportsIn(selectedExtremum.nodeId);
            const NodeId residualNode = static_cast<NodeId>(nodeParents.size());
            nodeParents.push_back(InvalidNode);       // Add X_s to implicit Xi, line 7.
            nodeAltitudes.push_back(supportAltitude); // eta(X_s), Eqs. 25, 33.
            for (PixelId pixel : selectedSupport.pixels) {
                if (smallestNodeMap[static_cast<std::size_t>(pixel)] == InvalidNode) {
                    smallestNodeMap[static_cast<std::size_t>(pixel)] = residualNode;
                }
            }
            while (!pendingChildren.empty()) {
                const NodeId childResidualNode = pendingChildren.front();
                nodeParents[static_cast<std::size_t>(childResidualNode)] = residualNode; // Eq. (32); remove child from implicit Xi, lines 9-10.
                dualPendingSupports.removePendingSupport(childResidualNode, selectedSupport.dualProperPartNode);
                primalPendingSupports.removePendingSupport(childResidualNode, selectedExtremum.nodeId);
            }
            return residualNode;
        }

        /**
         * @brief Transfers pending supports with their proper parts and adds the newly recorded support.
         * @tparam ProperPartTransfer Updater record containing source and target component nodes.
         * @param residualNode Newly recorded support to add to both pending support indices.
         * @param selectedPolarity Component-tree side of the selected leaf.
         * @param primalDestinationNode Live primal node receiving the selected leaf.
         * @param dualDestinationNode Live dual node containing the selected support after the update.
         * @param dualProperPartTransfers Whole dual proper-part transfers reported by the updater.
         */
        template <class ProperPartTransfer>
        void updatePendingSupports(NodeId residualNode, Polarity selectedPolarity, NodeId primalDestinationNode, NodeId dualDestinationNode,
                                   std::span<const ProperPartTransfer> dualProperPartTransfers) {
            PendingSupportsByProperPart& dualPendingSupports = pendingSupportsFor(oppositePolarity(selectedPolarity));
            // Pending supports follow whole proper-part transfers. The partial transfer
            // of X_s contains only its consumed children, already removed from Xi_s.
            for (const auto& transfer : dualProperPartTransfers)
                dualPendingSupports.transferPendingSupports(transfer.source, transfer.target);
            // Primal pruning moves the complete leaf to its parent. The dual updater
            // tracks the complete selected support through pixel moves and whole-part merges.
            pendingSupportsFor(selectedPolarity).addPendingSupport(residualNode, primalDestinationNode);
            dualPendingSupports.addPendingSupport(residualNode, dualDestinationNode);
        }

        /**
         * @brief Attaches the remaining supports and pixels to the root and returns the residual hierarchy.
         * @param image Source image defining the output pixel domain.
         * @param treeSemantics Residual-tree semantics transferred with the owned buffers.
         * @param rootAltitude Valuation of the terminal root.
         * @return Completed topology and valuation buffers.
         */
        [[nodiscard]] NativeHierarchyOwned<T> completeTree(const ImagePtr<T>& image, MorphologicalTreeSemantics treeSemantics, T rootAltitude) && {
            nodeParents[0] = 0;
            nodeAltitudes[0] = rootAltitude;
            // The remaining implicit Xi nodes become children of the terminal root.
            for (NodeId& parentNode : nodeParents)
                if (parentNode == InvalidNode)
                    parentNode = 0;
            for (NodeId& smallestNode : smallestNodeMap) {
                if (smallestNode == InvalidNode) {
                    smallestNode = 0;
                }
            }
            return makeNativeHierarchyOwned<T>(std::move(nodeParents), std::move(smallestNodeMap), std::move(nodeAltitudes), NodeId{0},
                                               GridDomain2D{image->getNumRows(), image->getNumColumns()}, std::move(treeSemantics));
        }

      private:
        /**
         * @brief Selects pending supports for one component-tree side.
         * @param polarity Component-tree side.
         * @return Mutable index of pending supports by proper part.
         */
        PendingSupportsByProperPart& pendingSupportsFor(Polarity polarity) {
            return polarity == Polarity::Maximum ? maxTreePendingSupports : minTreePendingSupports;
        }
    };

    RegularGridAdjacency2D adjacency_; ///< Shared adjacency used by both component trees.

    /**
     * @brief Tests whether a node represents an eligible nonterminal extremal flat zone.
     * @param tree Current component tree on the candidate side.
     * @param node Allocated node slot from initial construction or the leaf updater.
     * @return Whether the node is a live non-root leaf with a nonempty proper part.
     */
    [[nodiscard]] static bool isNonterminalLeaf(const ValuedMorphologicalTree<T>& tree, NodeId node) {
        const auto& topology = tree.topology();
        return detail::TreeKernelAccess::isAlive(topology, node) && node != topology.root() && detail::TreeKernelAccess::numChildren(topology, node) == 0;
    }
};

} // namespace mmcfilters
