#pragma once

#include "ComponentTreeBuilder.hpp"
#include "NativeHierarchyOwned.hpp"
#include "ValuedMorphologicalTree.hpp"
#include "detail/DualMinMaxTreeIncrementalFilterLeaf.hpp"
#include "detail/TreeKernelAccess.hpp"
#include "../utils/CommittedGridAccess.hpp"
#include "../utils/Common.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <limits>
#include <list>
#include <optional>
#include <set>
#include <span>
#include <utility>
#include <vector>

namespace mmcfilters {
namespace detail {
/** @brief Test access to support preparation and association invariants. @tparam T Image altitude type. */
template <AltitudeValue T> struct ExtremumResidualTreeTestAccess;
} // namespace detail

/**
 * Algorithm 1: construct the extremum residual tree through component trees.
 *
 * The active leaves implement Eq. (30); Xi and the two association relations
 * implement Eqs. (31-32). Only affected leaves and associations are updated.
 * Selected supports come from current component-tree proper parts. The dual
 * component-tree update is delegated to DualMinMaxTreeIncrementalFilterLeaf.
 * The min-tree and max-tree are temporary construction state. Each build
 * returns an independent owning native hierarchy. Image and graph inputs are
 * validated before evolution; construction stages consume established node,
 * support, and association domains without repeating defensive checks.
 *
 * @tparam T Image and node altitude type.
 */
template <AltitudeValue T> class ExtremumResidualTreeByCTsBuilder {
  public:
    using Altitude = T; ///< Scalar type shared by image samples and node altitudes.
    using ImagePointer = ImagePtr<T>; ///< Shared pointer to the input image.
    using Tree = ValuedMorphologicalTree<T>; ///< Temporary valued component tree modified during construction.

  private:
    friend struct detail::ExtremumResidualTreeTestAccess<T>;
    /** @brief Component-tree side containing the selected regional extremum. */
    enum class Polarity {
        Maximum, ///< The selected max-tree leaf is lowered while the min-tree is updated.
        Minimum ///< The selected min-tree leaf is raised while the max-tree is updated.
    };
    using NodeList = std::list<NodeId>; ///< Pending residual nodes associated with one component proper part.
    using Adjustment = detail::DualMinMaxTreeIncrementalFilterLeaf<T>; ///< Synchronized elementary-leveling updater.
    using TreeAccess = detail::TreeKernelAccess; ///< Direct reads of the component trees established by construction.

    /** Selected support and certificates reused by the incremental dual update. */
    struct LeafContext {
        std::vector<PixelId> supportPixels; ///< Pixels in the selected extremal flat zone.
        NodeId dualProperPartNode = InvalidNode; ///< Unique dual proper part containing the selected flat zone.
        NodeId wholeSupportSmallestNode = InvalidNode; ///< Dual node when its entire proper part is selected; otherwise InvalidNode.
    };
    /** Inverse association relation: X ◀ Y is stored in nodesByProperPart[Y]. */
    struct Associations {
        std::vector<NodeList> nodesByProperPart_; ///< Pending residual nodes indexed by component-tree proper part.
        std::vector<NodeList::iterator> positions_; ///< List position for each associated pending residual node.

      public:
        /** @brief Allocates empty association lists. @param numSlots Number of component-tree node slots. */
        explicit Associations(std::size_t numSlots) : nodesByProperPart_(numSlots) {}

        /** @brief Extends the lists for newly allocated component nodes. @param numSlots Required number of node slots. */
        void resize(std::size_t numSlots) {
            if (numSlots > nodesByProperPart_.size())
                nodesByProperPart_.resize(numSlots);
        }

        /** @brief Accesses an established association list. @param node Component-tree node slot. @return Mutable list of pending residual nodes. */
        NodeList& operator[](NodeId node) {
            return nodesByProperPart_[static_cast<std::size_t>(node)];
        }
        /** @brief Reads an established association list. @param node Component-tree node slot. @return Borrowed list of pending residual nodes. */
        const NodeList& operator[](NodeId node) const {
            return nodesByProperPart_[static_cast<std::size_t>(node)];
        }

        /**
         * @brief Associates a pending support with a component proper part.
         * @param residualNode Recorded residual support awaiting an inclusion parent.
         * @param componentNode Component-tree proper part containing the support.
         */
        void insert(NodeId residualNode, NodeId componentNode) {
            auto& nodes = (*this)[componentNode];
            if (static_cast<std::size_t>(residualNode) >= positions_.size())
                positions_.resize(static_cast<std::size_t>(residualNode) + 1);
            nodes.push_back(residualNode);
            positions_[static_cast<std::size_t>(residualNode)] = std::prev(nodes.end());
        }

        /**
         * @brief Removes a support after its inclusion parent has been recorded.
         * @param residualNode Associated residual support being consumed.
         * @param componentNode Component-tree proper part containing the association.
         */
        void erase(NodeId residualNode, NodeId componentNode) { (*this)[componentNode].erase(positions_[static_cast<std::size_t>(residualNode)]); }

        /**
         * @brief Moves all associations with a whole proper-part transfer, preserving list iterators.
         * @param source Component node whose proper part was transferred.
         * @param target Live component node receiving that proper part.
         */
        void transfer(NodeId source, NodeId target) {
            auto& destination = (*this)[target];
            destination.splice(destination.end(), (*this)[source]);
        }
    };

    /** Xi is implicit: non-root recorded nodes with parents_[node] == InvalidNode. */
    struct AssemblyState {
        Associations maxAssociations_; ///< Pending supports associated with max-tree proper parts.
        Associations minAssociations_; ///< Pending supports associated with min-tree proper parts.
        NodeId recordedSupportNode_ = InvalidNode; ///< Newly recorded support awaiting association updates.
        Polarity selectedPolarity_ = Polarity::Maximum; ///< Polarity of the newly recorded support.
        std::vector<NodeId> parents_{InvalidNode}; ///< Inclusion parents; slot zero is reserved for the terminal root.
        std::vector<NodeId> smallestNodeMap_; ///< First recorded support containing each pixel, or InvalidNode until assigned.
        std::vector<T> nodeAltitudes_{T{}}; ///< Proper-part valuations, with the terminal valuation reserved at slot zero.
        /**
         * @brief Allocates assembly buffers for a nonempty representable pixel domain.
         * @param numPixels Number of pixels and upper bound on residual nodes.
         * @param maxNodeSlots Initial number of max-tree node slots.
         * @param minNodeSlots Initial number of min-tree node slots.
         */
        AssemblyState(std::size_t numPixels, std::size_t maxNodeSlots, std::size_t minNodeSlots)
            : maxAssociations_(maxNodeSlots), minAssociations_(minNodeSlots), smallestNodeMap_(numPixels, InvalidNode) {
            if (numPixels == 0 || numPixels >= static_cast<std::size_t>(std::numeric_limits<NodeId>::max())) {
                throw std::length_error("Residual assembly requires a representable non-empty pixel domain.");
            }
            parents_.reserve(numPixels);
            nodeAltitudes_.reserve(numPixels);
        }
    };

    /** @brief Current extremal flat zone and its deterministic selection key. */
    struct Extremum {
        std::size_t supportCardinality = 0; ///< Cardinality of the current extremum support.
        PixelId spatialMinimum = InvalidPixel; ///< Smallest row-major pixel identifier in the support.
        Polarity polarity; ///< Component-tree side in which the flat zone is a leaf.
        NodeId nodeId; ///< Leaf identifier in that component tree.
    };
    /** @brief Orders current extrema by support cardinality and spatial minimum. */
    struct ExtremumLess {
        /** @brief Compares selection keys. @param lhs First candidate. @param rhs Second candidate. @return Whether lhs precedes rhs. */
        bool operator()(const Extremum& lhs, const Extremum& rhs) const {
            if (lhs.supportCardinality != rhs.supportCardinality)
                return lhs.supportCardinality < rhs.supportCardinality;
            return lhs.spatialMinimum < rhs.spatialMinimum;
        }
    };
    /** @brief Cached agenda membership for one component-tree node slot. */
    struct LeafSlot {
        bool active = false; ///< Whether the cached extremum is present in the agenda.
        Extremum extremum{0, InvalidPixel, Polarity::Maximum, InvalidNode}; ///< Last selection key and identity inserted for this slot.
    };
    /** @brief Ordered current extrema and membership caches for both component trees. */
    struct AgendaState {
        std::set<Extremum, ExtremumLess> extrema; ///< Active nonterminal extrema in deterministic selection order.
        std::vector<LeafSlot> maxSlots; ///< Max-tree agenda membership indexed by node slot.
        std::vector<LeafSlot> minSlots; ///< Min-tree agenda membership indexed by node slot.

        /**
         * @brief Creates an empty agenda with membership caches for both trees.
         * @param maxCount Initial number of max-tree node slots.
         * @param minCount Initial number of min-tree node slots.
         */
        AgendaState(std::size_t maxCount, std::size_t minCount) : maxSlots(maxCount), minSlots(minCount) {}
    };

    /** Temporary state referencing the component trees local to build(). */
    struct BuildState {
        const ImagePointer& image; ///< Original image retained for exact reconstruction checks.
        Tree& minTree; ///< Mutable min-tree of the current leveled image.
        Tree& maxTree; ///< Mutable max-tree of the current leveled image.
        Adjustment adjustment; ///< Incremental updater keeping both component trees synchronized.
        AgendaState agenda; ///< Selection agenda of current extremal flat zones.
        LeafContext context; ///< Reusable selected support and dual proper-part certificates.
        AssemblyState assembly; ///< Recorded residual supports and pending associations.

        /**
         * @brief Connects the construction state to the two local component trees.
         * @param input Original image whose valuation must be reconstructed.
         * @param min Mutable min-tree for elementary levelings.
         * @param max Mutable max-tree for elementary levelings.
         * @param adjacency Shared spatial adjacency of both trees.
         */
        BuildState(const ImagePointer& input, Tree& min, Tree& max, const RegularGridAdjacency2D& adjacency)
            : image(input), minTree(min), maxTree(max), adjustment(&min, &max, adjacency),
              agenda(max.topology().numInternalNodeSlots(), min.topology().numInternalNodeSlots()),
              assembly(input->getSize(), max.topology().numInternalNodeSlots(), min.topology().numInternalNodeSlots()) {}
    };

    /** @brief Selected extremum with tree references and values captured before leveling. */
    struct PreparedExtremum {
        Extremum selected; ///< Candidate removed by the next elementary leveling.
        Tree& primal; ///< Component tree in which the selected flat zone is a leaf.
        Tree& dual; ///< Opposite component tree updated after primal pruning.
        NodeId parent; ///< Primal parent receiving the selected leaf during pruning.
        T nodeAltitude; ///< eta(X_s), captured before elementary leveling (Eqs. 25, 33).
    };

    RegularGridAdjacency2D adjacency_; ///< Shared adjacency used by both component trees.

  public:
    /**
     * @brief Configures component-tree construction with a shared adjacency and row-major spatial ties.
     * @param adjacency Spatial adjacency matching the input image domain.
     */
    explicit ExtremumResidualTreeByCTsBuilder(RegularGridAdjacency2D adjacency) : adjacency_(std::move(adjacency)) {}

    /**
     * @brief Builds temporary min/max trees and records the extremum residual hierarchy.
     * @param image Original image reconstructed by the resulting proper-part valuations.
     * @return Independent native hierarchy owning its buffers and adjacency metadata.
     */
    [[nodiscard]] NativeHierarchyOwned<T> build(const ImagePointer& image) {
        requireResidualInputs(image, adjacency_);
        auto currentMin = detail::materializeNativeHierarchy<T>(
            ComponentTreeBuilder(adjacency_, false).build(image), std::nullopt, detail::MaterializedVersionPolicy::PreserveLinkedConstruction);
        auto currentMax = detail::materializeNativeHierarchy<T>(
            ComponentTreeBuilder(adjacency_, true).build(image), std::nullopt, detail::MaterializedVersionPolicy::PreserveLinkedConstruction);
        BuildState state(image, currentMin, currentMax, adjacency_);
        initialize(state); // Line 2.
        while (const auto selected = selectSmallestExtremum(state)) { // Lines 3-5.
            const auto extremum = prepareExtremum(state, *selected);
            // Lines 6-11: add X_s, attach its associated Xi_s children.
            static_cast<void>(recordSupport(state.assembly, state.context.supportPixels, selected->polarity, extremum.nodeAltitude,
                                            selected->nodeId, state.context.dualProperPartNode));
            levelExtremum(state, extremum); // Line 12.
            updateAssociations(state, extremum); // Line 13.
            refreshAffectedLeaves(state, extremum);
            // Line 14: the next step index follows from the number of recorded nodes.
        }
        return finalizeRoot(state); // Lines 16-18.
    }

    /** @brief Exposes the shared component-tree adjacency. @return Borrowed adjacency retained by this builder. */
    [[nodiscard]] const RegularGridAdjacency2D& adjacency() const noexcept { return adjacency_; }

  private:
    /**
     * @brief Tests whether a node represents an eligible nonterminal extremal flat zone.
     * @param tree Current component tree on the candidate side.
     * @param node Candidate node identifier.
     * @return Whether the node is a live non-root leaf with a nonempty proper part.
     */
    [[nodiscard]] static bool isResidualLeaf(const Tree& tree, NodeId node) {
        const auto& topology = tree.topology();
        return node >= 0 && topology.isNode(node) && topology.isAlive(node) && node != topology.root() && TreeAccess::numChildren(topology, node) == 0 &&
               TreeAccess::properPartCardinality(topology, node) > 0;
    }

    /** @brief Selects one side of the agenda cache. @param agenda Current selection state. @param polarity Component-tree side. @return Node-slot cache for that side. */
    static std::vector<LeafSlot>& slots(AgendaState& agenda, Polarity polarity) {
        return polarity == Polarity::Maximum ? agenda.maxSlots : agenda.minSlots;
    }

    /** @brief Removes an active cached candidate. @param agenda Current selection state. @param polarity Component-tree side. @param node Candidate node slot. */
    static void deactivate(AgendaState& agenda, Polarity polarity, NodeId node) {
        auto& side = slots(agenda, polarity);
        if (node < 0 || static_cast<std::size_t>(node) >= side.size()) return;
        auto& slot = side[node];
        if (slot.active)
            agenda.extrema.erase(slot.extremum);
        slot.active = false;
    }

    /**
     * @brief Updates agenda membership and the selection key after a local tree change.
     * @param agenda Current selection state.
     * @param tree Component tree containing the candidate.
     * @param polarity Component-tree side used to identify the candidate.
     * @param node Node slot whose leaf status or support may have changed.
     */
    static void refreshLeaf(AgendaState& agenda, const Tree& tree, Polarity polarity, NodeId node) {
        if (!isResidualLeaf(tree, node)) {
            deactivate(agenda, polarity, node);
            return;
        }
        auto& side = slots(agenda, polarity);
        if (static_cast<std::size_t>(node) >= side.size()) side.resize(static_cast<std::size_t>(node) + 1);
        auto& slot = side[node];
        const auto area = static_cast<std::size_t>(TreeAccess::properPartCardinality(tree.topology(), node));
        // Flat zones only merge: equal area and membership of the previous
        // spatial minimum certify that the cached key is still valid.
        if (slot.active && slot.extremum.supportCardinality == area &&
            TreeAccess::smallestNodeMap(tree.topology(), slot.extremum.spatialMinimum) == node) return;
        PixelId minimum = InvalidPixel;
        for (PixelId pixel : TreeAccess::properParts(tree.topology(), node))
            if (minimum == InvalidPixel || pixel < minimum) minimum = pixel;
        deactivate(agenda, polarity, node);
        const Extremum extremum{area, minimum, polarity, node};
        agenda.extrema.insert(extremum);
        slot = {true, extremum};
    }

    /** @brief Inserts all initial nonterminal extrema into the agenda. @param state Temporary state holding both initial trees. */
    static void initialize(BuildState& state) {
        for (NodeId node : state.maxTree.topology().aliveNodeIds()) refreshLeaf(state.agenda, state.maxTree, Polarity::Maximum, node);
        for (NodeId node : state.minTree.topology().aliveNodeIds()) refreshLeaf(state.agenda, state.minTree, Polarity::Minimum, node);
    }

    /** @brief Selects the first valid candidate, refreshing stale entries as needed. @param state Current component trees and agenda. @return Smallest extremum, or no value at termination. */
    [[nodiscard]] static std::optional<Extremum> selectSmallestExtremum(BuildState& state) {
        while (!state.agenda.extrema.empty()) {
            const Extremum selected = *state.agenda.extrema.begin();
            Tree& primal = selected.polarity == Polarity::Maximum ? state.maxTree : state.minTree;
            if (isResidualLeaf(primal, selected.nodeId) &&
                static_cast<std::size_t>(TreeAccess::properPartCardinality(primal.topology(), selected.nodeId)) == selected.supportCardinality &&
                TreeAccess::smallestNodeMap(primal.topology(), selected.spatialMinimum) == selected.nodeId) return selected;
            refreshLeaf(state.agenda, primal, selected.polarity, selected.nodeId);
        }
        return std::nullopt;
    }

    /**
     * @brief Collects the selected support and captures its altitude before leveling.
     * @param state Current trees and reusable support buffer.
     * @param selected Valid extremum returned by the agenda.
     * @return Primal/dual references, primal parent, and original extremum altitude.
     */
    [[nodiscard]] PreparedExtremum prepareExtremum(BuildState& state, const Extremum& selected) const {
        Tree& primal = selected.polarity == Polarity::Maximum ? state.maxTree : state.minTree;
        Tree& dual = selected.polarity == Polarity::Maximum ? state.minTree : state.maxTree;
        prepareResidualLeaf(state.context, selected.nodeId, primal, dual);
        const NodeId parent = TreeAccess::nodeParent(primal.topology(), selected.nodeId);
        return {selected, primal, dual, parent, primal.nodeAltitudes()[static_cast<std::size_t>(selected.nodeId)]};
    }

    /** @brief Prunes the selected leaf and updates the dual tree incrementally. @param state Current trees, support, and updater. @param extremum Prepared selection before leveling. */
    static void levelExtremum(BuildState& state, const PreparedExtremum& extremum) {
        const auto& context = state.context;
        if (extremum.selected.polarity == Polarity::Maximum)
            state.adjustment.pruneMaxLeafAndUpdateMinTree(extremum.selected.nodeId, context.supportPixels, context.dualProperPartNode,
                                                          context.wholeSupportSmallestNode);
        else
            state.adjustment.pruneMinLeafAndUpdateMaxTree(extremum.selected.nodeId, context.supportPixels, context.dualProperPartNode,
                                                          context.wholeSupportSmallestNode);
    }

    /** @brief Applies association transfers reported by the latest elementary leveling. @param state Current assembly and updater results. @param extremum Selection whose primal parent receives the support. */
    static void updateAssociations(BuildState& state, const PreparedExtremum& extremum) {
        updateAssociations(state.assembly, state.minTree, state.maxTree, extremum.parent, state.adjustment.getLastSelectedSupportNode(),
                           state.adjustment.getLastProperPartTransfers());
    }

    /** @brief Refreshes only primal and dual candidates touched by the latest leveling. @param state Current trees, updater, and agenda. @param extremum Selection identifying the affected primal leaf and parent. */
    static void refreshAffectedLeaves(BuildState& state, const PreparedExtremum& extremum) {
        const auto& selected = extremum.selected;
        const Polarity dualPolarity = selected.polarity == Polarity::Maximum ? Polarity::Minimum : Polarity::Maximum;
        for (NodeId node : state.adjustment.getLastCandidateNodes()) refreshLeaf(state.agenda, extremum.dual, dualPolarity, node);
        refreshLeaf(state.agenda, extremum.primal, selected.polarity, selected.nodeId);
        refreshLeaf(state.agenda, extremum.primal, selected.polarity, extremum.parent);
    }

    /** @brief Attaches remaining pending supports to the terminal full-domain root. @param state Completed evolution and residual assembly. @return Completed output buffers with the terminal valuation. */
    [[nodiscard]] NativeHierarchyOwned<T> finalizeRoot(BuildState& state) const {
        const T terminalAltitude = residualTerminalAltitude(state.minTree, state.maxTree);
        return finalize(std::move(state.assembly), state.image,
                           makeMorphologicalTreeSemantics(MorphologicalTreeKind::ExtremumResidualTree, SharedAdjacencyContext{adjacency_}),
                           terminalAltitude);
    }

    /** @brief Selects the association relation for one component-tree side. @param state Current residual assembly. @param polarity Component-tree side. @return Mutable inverse association relation. */
    static Associations& associations(AssemblyState& state, Polarity polarity) {
        return polarity == Polarity::Maximum ? state.maxAssociations_ : state.minAssociations_;
    }

    /**
     * @brief Reads a leaf proper part and establishes certificates for the dual update.
     * @param context Reusable support buffer and dual proper-part identifiers.
     * @param leaf Selected non-root primal leaf.
     * @param primal Component tree containing the selected extremum as a leaf.
     * @param dual Opposite component tree containing the support in one proper part.
     */
    static void prepareResidualLeaf(LeafContext& context, NodeId leaf, const Tree& primal, const Tree& dual) {
        const auto& dualTopology = dual.topology();
        context.supportPixels.clear();
        context.wholeSupportSmallestNode = InvalidNode;
        for (PixelId pixel : TreeAccess::properParts(primal.topology(), leaf)) {
            context.supportPixels.push_back(pixel);
        }
        // Selection establishes a nonempty extremal flat zone. Its pixels share
        // one dual proper part, so one lookup establishes the dual certificate.
        context.dualProperPartNode = TreeAccess::smallestNodeMap(dualTopology, context.supportPixels.front());
#ifndef NDEBUG
        for (PixelId pixel : context.supportPixels) {
            assert(TreeAccess::smallestNodeMap(dualTopology, pixel) == context.dualProperPartNode);
        }
#endif
        if (TreeAccess::properPartCardinality(dualTopology, context.dualProperPartNode) == static_cast<int>(context.supportPixels.size())) {
            context.wholeSupportSmallestNode = context.dualProperPartNode;
        }
    }
    /**
     * @brief Records a support, assigns its proper part, and attaches its pending children.
     * @param state Residual assembly before the selected elementary leveling.
     * @param support Pixels in the selected extremal flat zone.
     * @param polarity Side in which the selected flat zone is a leaf.
     * @param nodeAltitude Original valuation captured before leveling.
     * @param primalLeaf Selected leaf whose associations identify the pending children.
     * @param dualProperPart Dual proper part containing the complete selected support.
     * @return Identifier of the newly recorded residual node.
     */
    [[nodiscard]] static NodeId recordSupport(AssemblyState& state, std::span<const PixelId> support, Polarity polarity, T nodeAltitude,
                                            NodeId primalLeaf, NodeId dualProperPart) {
        // Property 6 gives each recorded support a new nonempty proper part;
        // together with the terminal root this bounds the node count by pixels.
        Associations& primal = associations(state, polarity);
        Associations& dual = associations(state, polarity == Polarity::Maximum ? Polarity::Minimum : Polarity::Maximum);
        auto& children = primal[primalLeaf];
        const NodeId node = static_cast<NodeId>(state.parents_.size());
        state.parents_.push_back(InvalidNode); // Add X_s to implicit Xi, line 7.
        state.nodeAltitudes_.push_back(nodeAltitude); // eta(X_s), Eqs. 25, 33.
        for (PixelId pixel : support) {
            if (state.smallestNodeMap_[static_cast<std::size_t>(pixel)] == InvalidNode) {
                state.smallestNodeMap_[static_cast<std::size_t>(pixel)] = node;
            }
        }
        while (!children.empty()) {
            const NodeId child = children.front();
            state.parents_[static_cast<std::size_t>(child)] = node; // Eq. (32); remove child from implicit Xi, lines 9-10.
            dual.erase(child, dualProperPart);
            primal.erase(child, primalLeaf);
        }
        state.recordedSupportNode_ = node;
        state.selectedPolarity_ = polarity;
        return node;
    }

    /**
     * @brief Moves pending associations with proper parts and associates the newly recorded support.
     * @tparam ProperPartTransfer Updater record containing source and target component nodes.
     * @param state Assembly containing one newly recorded support.
     * @param minTree Current min-tree after leveling.
     * @param maxTree Current max-tree after leveling.
     * @param primalProperPart Live primal node receiving the selected leaf.
     * @param dualProperPart Live dual node containing the selected support after the update.
     * @param dualTransfers Whole dual proper-part transfers reported by the updater.
     */
    template <class ProperPartTransfer>
    static void updateAssociations(AssemblyState& state, const Tree& minTree, const Tree& maxTree, NodeId primalProperPart, NodeId dualProperPart,
                            std::span<const ProperPartTransfer> dualTransfers) {
        state.maxAssociations_.resize(static_cast<std::size_t>(maxTree.topology().numInternalNodeSlots()));
        state.minAssociations_.resize(static_cast<std::size_t>(minTree.topology().numInternalNodeSlots()));
        Associations& dual = associations(state, state.selectedPolarity_ == Polarity::Maximum ? Polarity::Minimum : Polarity::Maximum);
        // Unchanged proper parts retain their associations. Every other pending
        // support follows a whole proper-part transfer. The partial transfer of
        // X_s contains only its consumed children, already removed from Xi_s.
        for (const auto& transfer : dualTransfers)
            dual.transfer(transfer.source, transfer.target);
        // Primal pruning moves the complete leaf to its parent. The dual updater
        // tracks the complete selected support through pixel moves and whole-part merges.
        associations(state, state.selectedPolarity_).insert(state.recordedSupportNode_, primalProperPart);
        dual.insert(state.recordedSupportNode_, dualProperPart);
        state.recordedSupportNode_ = InvalidNode;
    }

    /**
     * @brief Completes the terminal root using the established proper-part valuation.
     * @param state Completed assembly whose buffers are transferred to the result.
     * @param image Original image defining the pixel domain and reconstruction target.
     * @param semantics Residual-tree semantics transferred with the owned buffers.
     * @param terminalAltitude Common altitude of the terminal constant component trees.
     * @return Completed topology and valuation buffers.
     */
    [[nodiscard]] static NativeHierarchyOwned<T> finalize(AssemblyState&& state, const ImagePointer& image,
                                                                MorphologicalTreeSemantics semantics, T terminalAltitude) {
        state.parents_[0] = 0;
        state.nodeAltitudes_[0] = terminalAltitude;
        // The remaining implicit Xi nodes become children of the terminal root.
        for (NodeId& parent : state.parents_)
            if (parent == InvalidNode)
                parent = 0;
        for (NodeId& owner : state.smallestNodeMap_) {
            if (owner == InvalidNode) {
                owner = 0;
            }
        }
        return makeNativeHierarchyOwned<T>(std::move(state.parents_), std::move(state.smallestNodeMap_), std::move(state.nodeAltitudes_),
                                                  NodeId{0}, GridDomain2D{image->getNumRows(), image->getNumColumns()}, std::move(semantics));
    }

    /**
     * @brief Checks finite image samples and a matching connected domain once before evolution.
     * @param image Non-null input image with a nonempty representable domain.
     * @param adjacency Spatial adjacency whose domain must match the image.
     */
    static void requireResidualInputs(const ImagePtr<T>& image, const RegularGridAdjacency2D& adjacency) {
        if (!image)
            throw std::invalid_argument("Residual-tree construction requires a non-null image.");
        TreeAltitudeAlgorithms::validateFiniteImageAltitudes(image, "Residual-tree image");
        if (image->getNumRows() <= 0 || image->getNumColumns() <= 0 || image->getSize() <= 0 || adjacency.getNumRows() != image->getNumRows() ||
            adjacency.getNumColumns() != image->getNumColumns()) {
            throw std::invalid_argument("Residual-tree construction requires a non-empty matching image/adjacency domain.");
        }
        const auto size = static_cast<std::size_t>(image->getSize());
        std::vector<std::uint8_t> visited(size, 0);
        std::vector<PixelId> frontier{0};
        visited[0] = 1;
        std::size_t reached = 1;
        while (!frontier.empty()) {
            const PixelId pixel = frontier.back();
            frontier.pop_back();
            for (PixelId neighbor : detail::CommittedGridAccess::neighbors(adjacency, pixel)) {
                if (!visited[static_cast<std::size_t>(neighbor)]) {
                    visited[static_cast<std::size_t>(neighbor)] = 1;
                    ++reached;
                    frontier.push_back(neighbor);
                }
            }
        }
        if (reached != size)
            throw std::invalid_argument("Residual-tree adjacency must connect the image domain.");
    }

    /**
     * @brief Reads the common altitude of the terminal constant component trees.
     * @param minTree Terminal min-tree containing only its root.
     * @param maxTree Terminal max-tree containing only its root.
     * @return Common terminal altitude used to value the residual root.
     */
    [[nodiscard]] static T residualTerminalAltitude(const Tree& minTree, const Tree& maxTree) {
        const T altitude = minTree.nodeAltitudes()[static_cast<std::size_t>(minTree.topology().root())];
        assert(minTree.topology().numNodes() == 1 && maxTree.topology().numNodes() == 1);
        assert(maxTree.nodeAltitudes()[static_cast<std::size_t>(maxTree.topology().root())] == altitude);
        (void)maxTree;
        return altitude;
    }
};

} // namespace mmcfilters
