#pragma once

#include "../adjust/DualMinMaxTreeIncrementalFilterLeaf.hpp"
#include "../detail/NativeHierarchyValidationDetail.hpp"
#include "ResidualEvolution.hpp"

#include <algorithm>
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

namespace mmcfilters::sdrt {
namespace detail {
/** @brief Test access to support preparation and association invariants. @tparam T Image altitude type. */
template <AltitudeValue T> struct ExtremumResidualTreeTestAccess;
}

/**
 * Algorithm 1: construct the extremum residual tree through component trees.
 *
 * The active leaves implement Eq. (30); Xi and the two association relations
 * implement Eqs. (31-32). Only affected leaves and associations are updated.
 * Selected supports come from current component-tree proper parts. The dual
 * component-tree update is delegated to DualMinMaxTreeIncrementalFilterLeaf.
 * Result spans remain valid until the next build or result transfer.
 *
 * @tparam T Image and node altitude type.
 */
template <AltitudeValue T> class ExtremumResidualTreeByCTsBuilder {
  public:
    using Altitude = T; ///< Scalar type shared by image samples and node altitudes.
    using ImagePointer = ImagePtr<T>; ///< Shared pointer to the input image.
    using Tree = ValuedMorphologicalTree<T>; ///< Valued component tree consumed during construction.

  private:
    friend struct detail::ExtremumResidualTreeTestAccess<T>;
    using NodeList = std::list<NodeId>; ///< Pending residual nodes associated with one component proper part.
    using Adjustment = adjust::DualMinMaxTreeIncrementalFilterLeaf<T>; ///< Synchronized elementary-leveling updater.

    /** Completed output; construction buffers are owned by the local BuildState. */
    struct Result {
        int rows = 0; ///< Number of rows in the completed image domain.
        int columns = 0; ///< Number of columns in the completed image domain.
        NodeId root = InvalidNode; ///< Terminal root, or InvalidNode when no result is available.
        std::vector<NodeId> parents; ///< Inclusion parent for every recorded support.
        std::vector<NodeId> smallestNodeMap; ///< Smallest recorded support containing each pixel.
        std::vector<T> nodeAltitudes; ///< Original proper-part valuation of every recorded support.

        /** @brief Rejects access before a successful build or after result transfer. */
        void requireBuilt() const {
            if (root == InvalidNode)
                throw std::logic_error("Residual-tree builder has no completed result.");
        }

        /**
         * @brief Transfers the completed buffers and invalidates this result.
         * @param semantics Hierarchy kind, altitude order, and construction context attached by the factory.
         * @return Owning native hierarchy containing the completed topology and valuation.
         */
        [[nodiscard]] mmcfilters::detail::OwnedNativeHierarchy<T> take(MorphologicalTreeSemantics semantics) && {
            requireBuilt();
            auto hierarchy = mmcfilters::detail::makeOwnedNativeHierarchy<T>(std::move(parents), std::move(smallestNodeMap), std::move(nodeAltitudes), root,
                                                                             GridDomain2D{rows, columns}, std::move(semantics));
            root = InvalidNode;
            rows = columns = 0;
            return hierarchy;
        }
    };
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

        /** @brief Accesses a checked association list. @param node Component-tree node slot. @return Mutable list of pending residual nodes. */
        NodeList& at(NodeId node) {
            if (node < 0 || static_cast<std::size_t>(node) >= nodesByProperPart_.size())
                throw std::out_of_range("Residual association has an invalid component node.");
            return nodesByProperPart_[static_cast<std::size_t>(node)];
        }
        /** @brief Reads a checked association list. @param node Component-tree node slot. @return Borrowed list of pending residual nodes. */
        const NodeList& at(NodeId node) const {
            if (node < 0 || static_cast<std::size_t>(node) >= nodesByProperPart_.size())
                throw std::out_of_range("Residual association has an invalid component node.");
            return nodesByProperPart_[static_cast<std::size_t>(node)];
        }

        /**
         * @brief Associates a pending support with a component proper part.
         * @param residualNode Recorded residual support awaiting an inclusion parent.
         * @param componentNode Component-tree proper part containing the support.
         */
        void insert(NodeId residualNode, NodeId componentNode) {
            auto& nodes = at(componentNode);
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
        void erase(NodeId residualNode, NodeId componentNode) { at(componentNode).erase(positions_.at(static_cast<std::size_t>(residualNode))); }

        /**
         * @brief Moves all associations with a whole proper-part transfer, preserving list iterators.
         * @param source Component node whose proper part was transferred.
         * @param target Live component node receiving that proper part.
         */
        void transfer(NodeId source, NodeId target) {
            if (source == target)
                throw std::logic_error("Residual association transfer requires distinct nodes.");
            auto& destination = at(target);
            destination.splice(destination.end(), at(source));
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
        SelfDualResidualKey residualKey; ///< Support cardinality and spatial minimum used for selection.
        Polarity polarity; ///< Component-tree side in which the flat zone is a leaf.
        NodeId nodeId; ///< Leaf identifier in that component tree.
    };
    /** @brief Orders current extrema by support cardinality and spatial minimum. */
    struct ExtremumLess {
        SelfDualResidualOrder order; ///< Shared polarity-independent ordering rule.
        /** @brief Compares selection keys. @param lhs First candidate. @param rhs Second candidate. @return Whether lhs precedes rhs. */
        bool operator()(const Extremum& lhs, const Extremum& rhs) const { return order(lhs.residualKey, rhs.residualKey); }
    };
    /** @brief Cached agenda membership for one component-tree node slot. */
    struct LeafSlot {
        bool active = false; ///< Whether the cached extremum is present in the agenda.
        Extremum extremum{{}, Polarity::Maximum, InvalidNode}; ///< Last selection key and identity inserted for this slot.
    };
    /** @brief Ordered current extrema and membership caches for both component trees. */
    struct AgendaState {
        SpatialOrder spatialOrder; ///< Total pixel order used to break cardinality ties.
        std::set<Extremum, ExtremumLess> extrema; ///< Active nonterminal extrema in deterministic selection order.
        std::vector<LeafSlot> maxSlots; ///< Max-tree agenda membership indexed by node slot.
        std::vector<LeafSlot> minSlots; ///< Min-tree agenda membership indexed by node slot.

        /**
         * @brief Creates an empty agenda with membership caches for both trees.
         * @param order Total pixel order used by candidate keys.
         * @param maxCount Initial number of max-tree node slots.
         * @param minCount Initial number of min-tree node slots.
         */
        AgendaState(SpatialOrder order, std::size_t maxCount, std::size_t minCount)
            : spatialOrder(std::move(order)), extrema(ExtremumLess{SelfDualResidualOrder(spatialOrder)}), maxSlots(maxCount), minSlots(minCount) {}
    };

    /** One build owns all temporary state; only Result survives build(). */
    struct BuildState {
        const ImagePointer& image; ///< Original image retained for exact reconstruction checks.
        Tree& minTree; ///< Mutable min-tree of the current leveled image.
        Tree& maxTree; ///< Mutable max-tree of the current leveled image.
        Adjustment adjustment; ///< Incremental updater keeping both component trees synchronized.
        AgendaState agenda; ///< Selection agenda of current extremal flat zones.
        LeafContext context; ///< Reusable selected support and dual proper-part certificates.
        AssemblyState assembly; ///< Recorded residual supports and pending associations.

        /**
         * @brief Connects the temporary construction state to the two consumed trees.
         * @param input Original image whose valuation must be reconstructed.
         * @param min Mutable min-tree for elementary levelings.
         * @param max Mutable max-tree for elementary levelings.
         * @param adjacency Shared spatial adjacency of both trees.
         * @param order Total pixel order used for deterministic selection.
         */
        BuildState(const ImagePointer& input, Tree& min, Tree& max, const RegularGridAdjacency2D& adjacency, const SpatialOrder& order)
            : image(input), minTree(min), maxTree(max), adjustment(&min, &max, adjacency),
              agenda(order, max.topology().numInternalNodeSlots(), min.topology().numInternalNodeSlots()),
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
    SpatialOrder spatialOrder_; ///< Total pixel order defining deterministic selection ties.
    Result result_; ///< Completed output retained after temporary construction state is released.

  public:
    /**
     * @brief Configures component-tree construction with a shared adjacency and pixel order.
     * @param adjacency Spatial adjacency matching the input image domain.
     * @param options Spatial ordering policy for equal-cardinality extremal flat zones.
     */
    explicit ExtremumResidualTreeByCTsBuilder(RegularGridAdjacency2D adjacency, ExtremumResidualTreeOptions options = {})
        : adjacency_(std::move(adjacency)), spatialOrder_(std::move(options.spatialOrder)) {}

    /**
     * @brief Consumes the initial min/max trees and records the extremum residual hierarchy.
     * @param image Original image reconstructed by the resulting proper-part valuations.
     * @param minTree Min-tree of image, consumed and modified during construction.
     * @param maxTree Max-tree of image, consumed and modified during construction.
     */
    void build(const ImagePointer& image, Tree&& minTree, Tree&& maxTree) {
        result_ = {};
        requireResidualInputs(image, adjacency_, spatialOrder_);
        requireResidualComponentTree(minTree, image, MorphologicalTreeKind::MinTree);
        requireResidualComponentTree(maxTree, image, MorphologicalTreeKind::MaxTree);
        Tree currentMin(std::move(minTree));
        Tree currentMax(std::move(maxTree));
        BuildState state(image, currentMin, currentMax, adjacency_, spatialOrder_);
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
        result_ = finalizeRoot(state); // Lines 16-18.
    }

  public:
    /** @brief Reports the completed image height. @return Number of rows after a successful build. */
    [[nodiscard]] int rows() const {
        result_.requireBuilt();
        return result_.rows;
    }
    /** @brief Reports the completed image width. @return Number of columns after a successful build. */
    [[nodiscard]] int columns() const {
        result_.requireBuilt();
        return result_.columns;
    }
    /** @brief Identifies the terminal full-domain support. @return Root identifier of the completed residual tree. */
    [[nodiscard]] NodeId root() const {
        result_.requireBuilt();
        return result_.root;
    }
    /** @brief Exposes the shared component-tree adjacency. @return Borrowed adjacency retained by this builder. */
    [[nodiscard]] const RegularGridAdjacency2D& adjacency() const noexcept { return adjacency_; }
    /** @brief Exposes the deterministic pixel order. @return Borrowed spatial ordering policy. */
    [[nodiscard]] const SpatialOrder& spatialOrder() const noexcept { return spatialOrder_; }
    /** @brief Exposes the completed inclusion hierarchy. @return Borrowed parent identifiers indexed by residual node. */
    [[nodiscard]] std::span<const NodeId> parents() const {
        result_.requireBuilt();
        return result_.parents;
    }
    /** @brief Exposes the completed proper-part assignment. @return Borrowed smallest residual node for each pixel. */
    [[nodiscard]] std::span<const NodeId> smallestNodeMap() const {
        result_.requireBuilt();
        return result_.smallestNodeMap;
    }
    /** @brief Exposes the original proper-part valuations. @return Borrowed altitude buffer indexed by residual node. */
    [[nodiscard]] std::span<const T> nodeAltitudes() const {
        result_.requireBuilt();
        return result_.nodeAltitudes;
    }
    /**
     * @brief Transfers the completed result to the tree factory and invalidates result access.
     * @param semantics Hierarchy descriptor attached by the consuming factory.
     * @return Owning native hierarchy with buffers validated during construction.
     */
    [[nodiscard]] mmcfilters::detail::OwnedNativeHierarchy<T> takeValidatedHierarchy(MorphologicalTreeSemantics semantics) && {
        return std::move(result_).take(std::move(semantics));
    }

  private:
    /**
     * @brief Tests whether a node represents an eligible nonterminal extremal flat zone.
     * @param tree Current component tree on the candidate side.
     * @param node Candidate node identifier.
     * @return Whether the node is a live non-root leaf with a nonempty proper part.
     */
    [[nodiscard]] static bool isResidualLeaf(const Tree& tree, NodeId node) {
        const auto& topology = tree.topology();
        return node >= 0 && topology.isNode(node) && topology.isAlive(node) && node != topology.root() && topology.isLeaf(node) &&
               topology.properPartCardinality(node) > 0;
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
        if (slot.active && agenda.extrema.erase(slot.extremum) != 1)
            throw std::runtime_error("Residual leaf agenda lost an active entry.");
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
        const auto area = static_cast<std::size_t>(tree.topology().properPartCardinality(node));
        // Flat zones only merge: equal area and membership of the previous
        // spatial minimum certify that the cached key is still valid.
        if (slot.active && slot.extremum.residualKey.supportCardinality == area &&
            tree.topology().smallestNode(slot.extremum.residualKey.spatialMinimum) == node) return;
        PixelId minimum = InvalidPixel;
        for (PixelId pixel : tree.topology().properPart(node))
            if (minimum == InvalidPixel || agenda.spatialOrder.precedes(pixel, minimum)) minimum = pixel;
        deactivate(agenda, polarity, node);
        const Extremum extremum{{area, minimum}, polarity, node};
        if (!agenda.extrema.insert(extremum).second)
            throw std::runtime_error("Residual leaf agenda received duplicate support keys.");
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
                static_cast<std::size_t>(primal.topology().properPartCardinality(selected.nodeId)) == selected.residualKey.supportCardinality &&
                primal.topology().smallestNode(selected.residualKey.spatialMinimum) == selected.nodeId) return selected;
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
        if (selected.residualKey.supportCardinality != state.context.supportPixels.size())
            throw std::runtime_error("Residual leaf schedule has an inconsistent support cardinality.");
        const NodeId parent = primal.topology().parent(selected.nodeId);
        return {selected, primal, dual, parent, primal.nodeAltitude(selected.nodeId)};
    }

    /** @brief Prunes the selected leaf and updates the dual tree incrementally. @param state Current trees, support, and updater. @param extremum Prepared selection before leveling. */
    static void levelExtremum(BuildState& state, const PreparedExtremum& extremum) {
        const auto& context = state.context;
        // An empty adjacent-node span delegates boundary discovery to the updater.
        if (extremum.selected.polarity == Polarity::Maximum)
            state.adjustment.pruneMaxLeafAndUpdateMinTree(extremum.selected.nodeId, context.supportPixels, context.dualProperPartNode,
                                                          {}, context.wholeSupportSmallestNode);
        else
            state.adjustment.pruneMinLeafAndUpdateMaxTree(extremum.selected.nodeId, context.supportPixels, context.dualProperPartNode,
                                                          {}, context.wholeSupportSmallestNode);
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

    /** @brief Attaches remaining pending supports to the terminal full-domain root. @param state Completed evolution and residual assembly. @return Validated output buffers with the terminal valuation. */
    [[nodiscard]] Result finalizeRoot(BuildState& state) const {
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
        context.dualProperPartNode = InvalidNode;
        context.wholeSupportSmallestNode = InvalidNode;
        if (!isResidualLeaf(primal, leaf))
            throw std::runtime_error("Residual support preparation requires a non-root leaf with a non-empty proper part.");
        for (PixelId pixel : primal.topology().properPart(leaf)) {
            context.supportPixels.push_back(pixel);
            const NodeId smallestNode = dualTopology.smallestNode(pixel);
            if (smallestNode == InvalidNode || !dualTopology.isAlive(smallestNode)) {
                throw std::runtime_error("Residual leaf support has an invalid dual smallest node.");
            }
            // A selected extremal flat zone lies in one dual proper part.
            if (context.dualProperPartNode != InvalidNode && context.dualProperPartNode != smallestNode) {
                throw std::runtime_error("A selected flat zone must lie in one dual proper part.");
            }
            context.dualProperPartNode = smallestNode;
        }
        if (dualTopology.properPartCardinality(context.dualProperPartNode) == static_cast<int>(context.supportPixels.size())) {
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
        if (state.recordedSupportNode_ != InvalidNode)
            throw std::logic_error("Residual associations must be updated after each leveling.");
        if (support.empty()) {
            throw std::invalid_argument("Residual supports must be non-empty.");
        }
        if (state.parents_.size() >= state.smallestNodeMap_.size()) {
            throw std::length_error("Residual supports exhausted the pixel-domain node bound.");
        }
        Associations& primal = associations(state, polarity);
        Associations& dual = associations(state, polarity == Polarity::Maximum ? Polarity::Minimum : Polarity::Maximum);
        auto& children = primal.at(primalLeaf);
        static_cast<void>(dual.at(dualProperPart));
        const NodeId node = static_cast<NodeId>(state.parents_.size());
        state.parents_.push_back(InvalidNode); // Add X_s to implicit Xi, line 7.
        state.nodeAltitudes_.push_back(nodeAltitude); // eta(X_s), Eqs. 25, 33.
        bool hasProperPart = false;
        for (PixelId pixel : support) {
            if (pixel < 0 || static_cast<std::size_t>(pixel) >= state.smallestNodeMap_.size()) {
                throw std::out_of_range("Residual support contains an invalid pixel.");
            }
            if (state.smallestNodeMap_[static_cast<std::size_t>(pixel)] == InvalidNode) {
                state.smallestNodeMap_[static_cast<std::size_t>(pixel)] = node;
                hasProperPart = true;
            }
        }
        if (!hasProperPart)
            throw std::runtime_error("Residual proper parts must be non-empty (paper, Property 6).");
        while (!children.empty()) {
            const NodeId child = children.front();
            if (state.parents_[static_cast<std::size_t>(child)] != InvalidNode)
                throw std::runtime_error("A recorded residual support received a second parent.");
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
        if (state.recordedSupportNode_ == InvalidNode)
            throw std::logic_error("Residual association update requires a newly recorded support.");
        const Tree& primalTree = state.selectedPolarity_ == Polarity::Maximum ? maxTree : minTree;
        const Tree& dualTree = state.selectedPolarity_ == Polarity::Maximum ? minTree : maxTree;
        if (!primalTree.topology().isAlive(primalProperPart) || !dualTree.topology().isAlive(dualProperPart)) {
            throw std::invalid_argument("Residual association destinations must be live component-tree nodes.");
        }
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
     * @brief Completes the root and checks the descriptor and exact proper-part reconstruction.
     * @param state Completed assembly whose buffers are transferred to the result.
     * @param image Original image defining the pixel domain and reconstruction target.
     * @param semantics Residual-tree descriptor checked before buffer transfer.
     * @param terminalAltitude Common altitude of the terminal constant component trees.
     * @return Completed topology and valuation buffers.
     */
    [[nodiscard]] static Result finalize(AssemblyState&& state, const ImagePointer& image, MorphologicalTreeSemantics semantics, T terminalAltitude) {
        if (state.recordedSupportNode_ != InvalidNode)
            throw std::logic_error("Residual finalization requires completed association updates.");
        state.parents_[0] = 0;
        state.nodeAltitudes_[0] = terminalAltitude;
        // The remaining implicit Xi nodes become children of the terminal root.
        for (NodeId& parent : state.parents_)
            if (parent == InvalidNode)
                parent = 0;
        bool rootHasProperPart = false;
        for (NodeId& owner : state.smallestNodeMap_) {
            if (owner == InvalidNode) {
                owner = 0;
                rootHasProperPart = true;
            }
        }
        if (!rootHasProperPart)
            throw std::runtime_error("Residual root requires a non-empty proper part.");
        const NativeHierarchyView<T> hierarchy{
            state.parents_, state.smallestNodeMap_, state.nodeAltitudes_, 0, GridDomain2D{image->getNumRows(), image->getNumColumns()}, std::move(semantics)};
        mmcfilters::detail::NativeHierarchyValidation::validateDescriptorOnly(hierarchy);
        for (PixelId pixel = 0; pixel < image->getSize(); ++pixel) {
            const NodeId node = state.smallestNodeMap_[static_cast<std::size_t>(pixel)];
            if (state.nodeAltitudes_[static_cast<std::size_t>(node)] != (*image)[pixel]) {
                throw std::runtime_error("Residual proper-part valuation does not reconstruct the input.");
            }
        }
        return {image->getNumRows(),
                image->getNumColumns(),
                0,
                std::move(state.parents_),
                std::move(state.smallestNodeMap_),
                std::move(state.nodeAltitudes_)};
    }

    /**
     * @brief Checks finite image samples, a matching connected domain, and the spatial order once before evolution.
     * @param image Non-null input image with a nonempty representable domain.
     * @param adjacency Spatial adjacency whose domain must match the image.
     * @param order Total pixel order valid on the image domain.
     */
    static void requireResidualInputs(const ImagePtr<T>& image, const RegularGridAdjacency2D& adjacency, const SpatialOrder& order) {
        if (!image)
            throw std::invalid_argument("Residual-tree construction requires a non-null image.");
        TreeAltitudeAlgorithms::validateFiniteImageAltitudes(image, "Residual-tree image");
        if (image->getNumRows() <= 0 || image->getNumColumns() <= 0 || image->getSize() <= 0 || adjacency.getNumRows() != image->getNumRows() ||
            adjacency.getNumColumns() != image->getNumColumns()) {
            throw std::invalid_argument("Residual-tree construction requires a non-empty matching image/adjacency domain.");
        }
        const auto size = static_cast<std::size_t>(image->getSize());
        order.validateDomainSize(size);
        std::vector<std::uint8_t> visited(size, 0);
        std::vector<PixelId> frontier{0};
        visited[0] = 1;
        std::size_t reached = 1;
        while (!frontier.empty()) {
            const PixelId pixel = frontier.back();
            frontier.pop_back();
            for (PixelId neighbor : adjacency.getNeighborIndices(pixel)) {
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
     * @brief Checks the initial component-tree kind, domain, and altitude buffer shape.
     * @param tree Component tree supplied for consumption by the builder.
     * @param image Input image defining the expected domain.
     * @param kind Required min-tree or max-tree kind.
     */
    static void requireResidualComponentTree(const Tree& tree, const ImagePtr<T>& image, MorphologicalTreeKind kind) {
        const auto& topology = tree.topology();
        if (topology.kind() != kind || topology.numRows() != image->getNumRows() || topology.numColumns() != image->getNumColumns() ||
            topology.numPixels() != image->getSize()) {
            throw std::invalid_argument("Residual-tree component tree has an incompatible kind or domain.");
        }
        tree.validateNodeAltitudeBufferShape();
    }

    /**
     * @brief Checks that evolution ended at two matching constant component trees.
     * @param minTree Terminal min-tree containing only its root.
     * @param maxTree Terminal max-tree containing only its root.
     * @return Common terminal altitude used to value the residual root.
     */
    [[nodiscard]] static T residualTerminalAltitude(const Tree& minTree, const Tree& maxTree) {
        const T altitude = minTree.nodeAltitude(minTree.topology().root());
        if (minTree.topology().numNodes() != 1 || maxTree.topology().numNodes() != 1 || maxTree.nodeAltitude(maxTree.topology().root()) != altitude) {
            throw std::runtime_error("Residual evolution did not terminate at two matching constant trees.");
        }
        return altitude;
    }
};

} // namespace mmcfilters::sdrt
