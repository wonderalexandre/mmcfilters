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
// Defined only by the tests that check Eq. (31) after each elementary leveling.
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
 */
template <AltitudeValue T> class ExtremumResidualTreeByCTsBuilder {
  public:
    using Altitude = T;
    using ImagePointer = ImagePtr<T>;
    using Tree = ValuedMorphologicalTree<T>;

  private:
    friend struct detail::ExtremumResidualTreeTestAccess<T>;
    using NodeList = std::list<NodeId>;
    using Adjustment = adjust::DualMinMaxTreeIncrementalFilterLeaf<T>;

    /** Completed output; construction buffers are owned by the local BuildState. */
    struct Result {
        int rows = 0;
        int columns = 0;
        NodeId root = InvalidNode;
        std::vector<NodeId> parents;
        std::vector<NodeId> smallestNodeMap;
        std::vector<T> nodeAltitudes;

        void requireBuilt() const {
            if (root == InvalidNode)
                throw std::logic_error("Residual-tree builder has no completed result.");
        }

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
        std::vector<PixelId> supportPixels;
        NodeId dualProperPartNode = InvalidNode; ///< Unique dual proper part containing the selected flat zone.
        NodeId wholeSupportSmallestNode = InvalidNode;
    };
    /** Inverse association relation: X ◀ Y is stored in nodesByProperPart[Y]. */
    struct Associations {
        std::vector<NodeList> nodesByProperPart_;
        std::vector<NodeList::iterator> positions_;

      public:
        explicit Associations(std::size_t numSlots) : nodesByProperPart_(numSlots) {}

        void resize(std::size_t numSlots) {
            if (numSlots > nodesByProperPart_.size())
                nodesByProperPart_.resize(numSlots);
        }

        NodeList& at(NodeId node) {
            if (node < 0 || static_cast<std::size_t>(node) >= nodesByProperPart_.size())
                throw std::out_of_range("Residual association has an invalid component node.");
            return nodesByProperPart_[static_cast<std::size_t>(node)];
        }
        const NodeList& at(NodeId node) const {
            if (node < 0 || static_cast<std::size_t>(node) >= nodesByProperPart_.size())
                throw std::out_of_range("Residual association has an invalid component node.");
            return nodesByProperPart_[static_cast<std::size_t>(node)];
        }

        void insert(NodeId residualNode, NodeId componentNode) {
            auto& nodes = at(componentNode);
            if (static_cast<std::size_t>(residualNode) >= positions_.size())
                positions_.resize(static_cast<std::size_t>(residualNode) + 1);
            nodes.push_back(residualNode);
            positions_[static_cast<std::size_t>(residualNode)] = std::prev(nodes.end());
        }

        // Every child of the selected leaf also belongs to its unique dual
        // proper part. Whole-list splices preserve all stored iterators.
        void erase(NodeId residualNode, NodeId componentNode) { at(componentNode).erase(positions_.at(static_cast<std::size_t>(residualNode))); }

        void transfer(NodeId source, NodeId target) {
            if (source == target)
                throw std::logic_error("Residual association transfer requires distinct nodes.");
            auto& destination = at(target);
            destination.splice(destination.end(), at(source));
        }
    };

    /** Xi is implicit: non-root recorded nodes with parents_[node] == InvalidNode. */
    struct AssemblyState {
        Associations maxAssociations_;
        Associations minAssociations_;
        NodeId recordedSupportNode_ = InvalidNode;
        Polarity selectedPolarity_ = Polarity::Maximum;
        std::vector<NodeId> parents_{InvalidNode};
        std::vector<NodeId> smallestNodeMap_;
        std::vector<T> nodeAltitudes_{T{}};
        AssemblyState(std::size_t numPixels, std::size_t maxNodeSlots, std::size_t minNodeSlots)
            : maxAssociations_(maxNodeSlots), minAssociations_(minNodeSlots), smallestNodeMap_(numPixels, InvalidNode) {
            if (numPixels == 0 || numPixels >= static_cast<std::size_t>(std::numeric_limits<NodeId>::max())) {
                throw std::length_error("Residual assembly requires a representable non-empty pixel domain.");
            }
            parents_.reserve(numPixels);
            nodeAltitudes_.reserve(numPixels);
        }
    };

    struct Extremum {
        SelfDualResidualKey residualKey;
        Polarity polarity;
        NodeId nodeId;
    };
    struct ExtremumLess {
        SelfDualResidualOrder order;
        bool operator()(const Extremum& lhs, const Extremum& rhs) const { return order(lhs.residualKey, rhs.residualKey); }
    };
    struct LeafSlot {
        bool active = false;
        Extremum extremum{{}, Polarity::Maximum, InvalidNode};
    };
    struct AgendaState {
        SpatialOrder spatialOrder;
        std::set<Extremum, ExtremumLess> extrema;
        std::vector<LeafSlot> maxSlots;
        std::vector<LeafSlot> minSlots;

        AgendaState(SpatialOrder order, std::size_t maxCount, std::size_t minCount)
            : spatialOrder(std::move(order)), extrema(ExtremumLess{SelfDualResidualOrder(spatialOrder)}), maxSlots(maxCount), minSlots(minCount) {}
    };

    /** One build owns all temporary state; only Result survives build(). */
    struct BuildState {
        const ImagePointer& image;
        Tree& minTree;
        Tree& maxTree;
        Adjustment adjustment;
        AgendaState agenda;
        LeafContext context;
        AssemblyState assembly;

        BuildState(const ImagePointer& input, Tree& min, Tree& max, const RegularGridAdjacency2D& adjacency, const SpatialOrder& order)
            : image(input), minTree(min), maxTree(max), adjustment(&min, &max, adjacency),
              agenda(order, max.topology().numInternalNodeSlots(), min.topology().numInternalNodeSlots()),
              assembly(input->getSize(), max.topology().numInternalNodeSlots(), min.topology().numInternalNodeSlots()) {}
    };

    struct PreparedExtremum {
        Extremum selected;
        Tree& primal;
        Tree& dual;
        NodeId parent;
        T nodeAltitude; ///< eta(X_s), captured before elementary leveling (Eqs. 25, 33).
    };

    RegularGridAdjacency2D adjacency_;
    SpatialOrder spatialOrder_;
    Result result_;

  public:
    explicit ExtremumResidualTreeByCTsBuilder(RegularGridAdjacency2D adjacency, ExtremumResidualTreeOptions options = {})
        : adjacency_(std::move(adjacency)), spatialOrder_(std::move(options.spatialOrder)) {}

    /** Consume the initial min/max trees (Algorithm 1, line 1). */
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
    [[nodiscard]] int rows() const {
        result_.requireBuilt();
        return result_.rows;
    }
    [[nodiscard]] int columns() const {
        result_.requireBuilt();
        return result_.columns;
    }
    [[nodiscard]] NodeId root() const {
        result_.requireBuilt();
        return result_.root;
    }
    [[nodiscard]] const RegularGridAdjacency2D& adjacency() const noexcept { return adjacency_; }
    [[nodiscard]] const SpatialOrder& spatialOrder() const noexcept { return spatialOrder_; }
    [[nodiscard]] std::span<const NodeId> parents() const {
        result_.requireBuilt();
        return result_.parents;
    }
    [[nodiscard]] std::span<const NodeId> smallestNodeMap() const {
        result_.requireBuilt();
        return result_.smallestNodeMap;
    }
    [[nodiscard]] std::span<const T> nodeAltitudes() const {
        result_.requireBuilt();
        return result_.nodeAltitudes;
    }
    [[nodiscard]] mmcfilters::detail::OwnedNativeHierarchy<T> takeValidatedHierarchy(MorphologicalTreeSemantics semantics) && {
        return std::move(result_).take(std::move(semantics));
    }

  private:
    [[nodiscard]] static bool isResidualLeaf(const Tree& tree, NodeId node) {
        const auto& topology = tree.topology();
        return node >= 0 && topology.isNode(node) && topology.isAlive(node) && node != topology.root() && topology.isLeaf(node) &&
               topology.properPartCardinality(node) > 0;
    }

    static std::vector<LeafSlot>& slots(AgendaState& agenda, Polarity polarity) {
        return polarity == Polarity::Maximum ? agenda.maxSlots : agenda.minSlots;
    }

    static void deactivate(AgendaState& agenda, Polarity polarity, NodeId node) {
        auto& side = slots(agenda, polarity);
        if (node < 0 || static_cast<std::size_t>(node) >= side.size()) return;
        auto& slot = side[node];
        if (slot.active && agenda.extrema.erase(slot.extremum) != 1)
            throw std::runtime_error("Residual leaf agenda lost an active entry.");
        slot.active = false;
    }

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

    static void initialize(BuildState& state) {
        for (NodeId node : state.maxTree.topology().aliveNodeIds()) refreshLeaf(state.agenda, state.maxTree, Polarity::Maximum, node);
        for (NodeId node : state.minTree.topology().aliveNodeIds()) refreshLeaf(state.agenda, state.minTree, Polarity::Minimum, node);
    }

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

    [[nodiscard]] PreparedExtremum prepareExtremum(BuildState& state, const Extremum& selected) const {
        Tree& primal = selected.polarity == Polarity::Maximum ? state.maxTree : state.minTree;
        Tree& dual = selected.polarity == Polarity::Maximum ? state.minTree : state.maxTree;
        prepareResidualLeaf(state.context, selected.nodeId, primal, dual);
        if (selected.residualKey.supportCardinality != state.context.supportPixels.size())
            throw std::runtime_error("Residual leaf schedule has an inconsistent support cardinality.");
        const NodeId parent = primal.topology().parent(selected.nodeId);
        return {selected, primal, dual, parent, primal.nodeAltitude(selected.nodeId)};
    }

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

    static void updateAssociations(BuildState& state, const PreparedExtremum& extremum) {
        updateAssociations(state.assembly, state.minTree, state.maxTree, extremum.parent, state.adjustment.getLastSelectedSupportNode(),
                           state.adjustment.getLastProperPartTransfers());
    }

    static void refreshAffectedLeaves(BuildState& state, const PreparedExtremum& extremum) {
        const auto& selected = extremum.selected;
        const Polarity dualPolarity = selected.polarity == Polarity::Maximum ? Polarity::Minimum : Polarity::Maximum;
        for (NodeId node : state.adjustment.getLastCandidateNodes()) refreshLeaf(state.agenda, extremum.dual, dualPolarity, node);
        refreshLeaf(state.agenda, extremum.primal, selected.polarity, selected.nodeId);
        refreshLeaf(state.agenda, extremum.primal, selected.polarity, extremum.parent);
    }

    [[nodiscard]] Result finalizeRoot(BuildState& state) const {
        const T terminalAltitude = residualTerminalAltitude(state.minTree, state.maxTree);
        return finalize(std::move(state.assembly), state.image,
                           makeMorphologicalTreeSemantics(MorphologicalTreeKind::ExtremumResidualTree, SharedAdjacencyContext{adjacency_}),
                           terminalAltitude);
    }

    static Associations& associations(AssemblyState& state, Polarity polarity) {
        return polarity == Polarity::Maximum ? state.maxAssociations_ : state.minAssociations_;
    }

    /** Read the proper part of the leaf and prepare the incremental dual update. */
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

    /** Validate the domain once; no full image/tree validation occurs inside the leveling loop. */
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

    static void requireResidualComponentTree(const Tree& tree, const ImagePtr<T>& image, MorphologicalTreeKind kind) {
        const auto& topology = tree.topology();
        if (topology.kind() != kind || topology.numRows() != image->getNumRows() || topology.numColumns() != image->getNumColumns() ||
            topology.numPixels() != image->getSize()) {
            throw std::invalid_argument("Residual-tree component tree has an incompatible kind or domain.");
        }
        tree.validateNodeAltitudeBufferShape();
    }

    [[nodiscard]] static T residualTerminalAltitude(const Tree& minTree, const Tree& maxTree) {
        const T altitude = minTree.nodeAltitude(minTree.topology().root());
        if (minTree.topology().numNodes() != 1 || maxTree.topology().numNodes() != 1 || maxTree.nodeAltitude(maxTree.topology().root()) != altitude) {
            throw std::runtime_error("Residual evolution did not terminate at two matching constant trees.");
        }
        return altitude;
    }
};

} // namespace mmcfilters::sdrt
