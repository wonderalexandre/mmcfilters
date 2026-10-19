#include "support/TestSupport.hpp"

#include "mmcfilters/trees/saliency/HierarchySaliencyMap.hpp"
#include "mmcfilters/trees/ExtremumResidualTreeByCTsBuilder.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

// Test-only access keeps the stepwise Eq. (31) oracle outside the public API.
namespace mmcfilters::detail {
/** @brief Exposes private construction stages to the stepwise association oracle. @tparam T Image altitude type. */
template <AltitudeValue T> struct ExtremumResidualTreeTestAccess {
    using Builder = ExtremumResidualTreeByCTsBuilder<T>; ///< Production builder checked by the oracle.
    using Polarity = typename Builder::Polarity; ///< Private component-tree side exposed only to the test oracle.
    using Extremum = typename Builder::Extremum; ///< Private candidate fields exposed only to ordering tests.
    using ExtremumLess = typename Builder::ExtremumLess; ///< Production candidate comparator exposed only to ordering tests.
    using Context = typename Builder::LeafContext; ///< Selected support and dual proper-part certificates.

    /**
     * @brief Collects a selected leaf support through the production helper.
     * @param context Reusable selected-support buffer.
     * @param leaf Selected primal leaf identifier.
     * @param primal Component tree containing that leaf.
     * @param dual Opposite component tree used to establish proper-part certificates.
     */
    static void prepareLeaf(Context& context, NodeId leaf, const ValuedMorphologicalTree<T>& primal,
                            const ValuedMorphologicalTree<T>& dual) {
        Builder::prepareResidualLeaf(context, leaf, primal, dual);
    }
    /** @brief Test facade for recording supports and inspecting pending associations. */
    struct Assembler {
        typename Builder::AssemblyState state; ///< Production assembly state inspected after each test step.
        /** @brief Allocates an isolated test assembly. @param pixels Number of domain pixels. @param maxSlots Initial max-tree node slots. @param minSlots Initial min-tree node slots. */
        Assembler(std::size_t pixels, std::size_t maxSlots, std::size_t minSlots) : state(pixels, maxSlots, minSlots) {}
        /**
         * @brief Records a test support through the production assembly routine.
         * @param support Pixels in the selected flat zone.
         * @param polarity Component-tree side of the selected leaf.
         * @param nodeAltitude Original valuation of the support.
         * @param primalLeaf Selected primal leaf identifier.
         * @param dualProperPart Dual proper part containing the support.
         * @return Newly recorded residual node identifier.
         */
        NodeId recordSupport(std::span<const PixelId> support, Polarity polarity, T nodeAltitude, NodeId primalLeaf, NodeId dualProperPart) {
            return Builder::recordSupport(state, support, polarity, nodeAltitude, primalLeaf, dualProperPart);
        }
        /** @brief Counts recorded nonterminal supports. @return Number of residual nodes excluding the reserved root. */
        std::size_t numRecordedSupports() const { return state.parents_.size() - 1; }
        /** @brief Enumerates the implicit pending set. @return Recorded non-root nodes whose parent remains undefined. */
        std::vector<NodeId> pendingNodes() const {
            std::vector<NodeId> pending;
            for (NodeId node = 1; node < static_cast<NodeId>(state.parents_.size()); ++node)
                if (state.parents_[static_cast<std::size_t>(node)] == InvalidNode)
                    pending.push_back(node);
            return pending;
        }
        /** @brief Reads one inverse association list. @param polarity Component-tree side. @param componentNode Proper-part node identifier. @return Borrowed list of associated pending supports. */
        const auto& associatedNodes(Polarity polarity, NodeId componentNode) const {
            return (polarity == Polarity::Maximum ? state.maxAssociations_ : state.minAssociations_)[componentNode];
        }
        /**
         * @brief Applies production association updates to the test state.
         * @tparam Transfer Whole proper-part transfer record type.
         * @param minTree Current min-tree after leveling.
         * @param maxTree Current max-tree after leveling.
         * @param primal Live primal destination of the selected support.
         * @param dual Live dual destination of the selected support.
         * @param transfers Whole dual proper-part transfers reported by the updater.
         */
        template <class Transfer>
        void updateAssociations(const ValuedMorphologicalTree<T>& minTree, const ValuedMorphologicalTree<T>& maxTree,
                                NodeId primal, NodeId dual, std::span<const Transfer> transfers) {
            Builder::updateAssociations(state, minTree, maxTree, primal, dual, transfers);
        }
        /**
         * @brief Completes the test assembly through the production reconstruction checks.
         * @param image Original image defining the expected valuation.
         * @param adjacency Shared construction adjacency attached to the descriptor.
         * @param terminalAltitude Common terminal component-tree altitude.
         * @return Completed residual topology and valuation buffers.
         */
        auto finalize(const ImagePtr<T>& image, const RegularGridAdjacency2D& adjacency, T terminalAltitude) && {
            return Builder::finalize(std::move(state), image,
                                    makeMorphologicalTreeSemantics(MorphologicalTreeKind::ExtremumResidualTree, SharedAdjacencyContext{adjacency}),
                                    terminalAltitude);
        }
    };
};
} // namespace mmcfilters::detail

template <class T> using BuilderAccess = mmcfilters::detail::ExtremumResidualTreeTestAccess<T>;
template <class T> using TestAssembler = typename BuilderAccess<T>::Assembler;

namespace {

template <class Engine>
concept HasInfinityPixelAccessor = requires(const Engine& engine) { engine.infinityPixel(); };

template <class Assembler>
concept LvalueFinalizableAssembler = requires(Assembler& assembler, const ImagePtr<int>& image, const RegularGridAdjacency2D& adjacency) {
    assembler.finalize(image, adjacency, 0);
};

template <class Assembler>
concept RvalueFinalizableAssembler = requires(Assembler&& assembler, const ImagePtr<int>& image, const RegularGridAdjacency2D& adjacency) {
    std::move(assembler).finalize(image, adjacency, 0);
};

static_assert(!HasInfinityPixelAccessor<ExtremumResidualTreeByCTsBuilder<std::uint8_t>>, "extremum construction must not expose infinity-pixel state");
static_assert(!LvalueFinalizableAssembler<TestAssembler<int>>, "residual support assembly must not be finalizable through an lvalue");
static_assert(RvalueFinalizableAssembler<TestAssembler<int>>, "residual support assembly must be finalizable through an rvalue");
static_assert(std::is_constructible_v<ExtremumResidualTreeByCTsBuilder<std::uint8_t>, RegularGridAdjacency2D>,
              "extremum builders require only the shared adjacency");

template <class T> std::vector<NodeId> parentBuffer(const ValuedMorphologicalTree<T>& tree) {
    std::vector<NodeId> parent(static_cast<std::size_t>(tree.topology().numInternalNodeSlots()));
    for (NodeId node : tree.topology().aliveNodeIds()) {
        parent[static_cast<std::size_t>(node)] = tree.topology().parent(node);
    }
    return parent;
}

template <class T> std::vector<NodeId> smallestNodeMapBuffer(const ValuedMorphologicalTree<T>& tree) {
    std::vector<NodeId> smallestNodeMap(static_cast<std::size_t>(tree.topology().numPixels()));
    for (PixelId pixel = 0; pixel < tree.topology().numPixels(); ++pixel) {
        smallestNodeMap[static_cast<std::size_t>(pixel)] = tree.topology().smallestNode(pixel);
    }
    return smallestNodeMap;
}

template <class T> std::vector<T> altitudeBuffer(const ValuedMorphologicalTree<T>& tree) {
    std::vector<T> altitude(static_cast<std::size_t>(tree.topology().numInternalNodeSlots()));
    for (NodeId node : tree.topology().aliveNodeIds()) {
        altitude[static_cast<std::size_t>(node)] = tree.nodeAltitude(node);
    }
    return altitude;
}

template <class T> ImagePtr<T> makeTypedImage(int rows, int columns, std::initializer_list<T> values) {
    requireEqual(static_cast<int>(values.size()), rows * columns, "typed image buffer size");
    auto image = Image<T>::create(rows, columns);
    PixelId pixel = 0;
    for (T value : values) {
        (*image)[pixel++] = value;
    }
    return image;
}

template <class T> void requireExactReconstruction(const ValuedMorphologicalTree<T>& tree, const ImagePtr<T>& image, const std::string& label) {
    for (PixelId pixel = 0; pixel < image->getSize(); ++pixel) {
        const NodeId smallestNode = tree.topology().smallestNode(pixel);
        requireEqual(tree.nodeAltitude(smallestNode), (*image)[pixel], label);
    }
}

template <AltitudeValue T>
void requireRebuiltComponentTrees(const ValuedMorphologicalTree<T>& minTree, const ValuedMorphologicalTree<T>& maxTree,
                                 const ImagePtr<T>& leveledImage, const RegularGridAdjacency2D& adjacency, const std::string& label) {
    const auto components = [](const ValuedMorphologicalTree<T>& tree) {
        std::vector<std::pair<T, std::vector<PixelId>>> result;
        for (NodeId node : tree.topology().aliveNodeIds()) {
            const auto support = tree.topology().nodeSupport(node);
            std::vector<PixelId> pixels(support.begin(), support.end());
            std::ranges::sort(pixels);
            result.emplace_back(tree.nodeAltitude(node), std::move(pixels));
        }
        std::ranges::sort(result);
        return result;
    };
    const auto rebuiltMin = MorphologicalTreeFactory::createMinTree(leveledImage, adjacency);
    const auto rebuiltMax = MorphologicalTreeFactory::createMaxTree(leveledImage, adjacency);
    require(components(minTree) == components(rebuiltMin), label + ": min-tree components match a fresh construction");
    require(components(maxTree) == components(rebuiltMax), label + ": max-tree components match a fresh construction");
    for (const auto* tree : {&minTree, &maxTree}) {
        require(tree->topology().validateConnectedRootedTreeResult().ok, label + ": connected rooted topology");
        tree->validateMonotoneNodeAltitudes();
        requireExactReconstruction(*tree, leveledImage, label + ": exact leveled image");
    }
}

ImageUInt8Ptr makeRadixImage(int code) {
    auto image = ImageUInt8::create(2, 3);
    for (PixelId pixel = 0; pixel < 6; ++pixel) {
        (*image)[pixel] = static_cast<std::uint8_t>(code % 3);
        code /= 3;
    }
    return image;
}

void requireResidualSemantics(const MorphologicalTree& tree, const RegularGridAdjacency2D& adjacency, MorphologicalTreeKind expectedKind,
                              const std::string& label) {
    require(tree.kind() == expectedKind, label + ": descriptive kind");
    require(tree.nodeAltitudeOrder() == NodeAltitudeOrder::Unconstrained, label + ": unconstrained altitude order");
    const RegularGridAdjacency2D* stored = nullptr;
    if (const auto* context = tree.sharedAdjacencyContext()) {
        stored = &context->adjacency;
    }
    require(stored != nullptr, label + ": retained construction adjacency");
    requireEqual(stored->getNumRows(), adjacency.getNumRows(), label + ": adjacency rows");
    requireEqual(stored->getNumColumns(), adjacency.getNumColumns(), label + ": adjacency columns");
}

void testFactoryAndSaliencyIntegration() {
    const auto image = makeImage(3, 4,
                                 {
                                     0,
                                     1,
                                     2,
                                     1,
                                     1,
                                     4,
                                     1,
                                     0,
                                     2,
                                     1,
                                     3,
                                     1,
                                 });
    const RegularGridAdjacency2D adjacency(3, 4, 1.0);

    const auto extremum = MorphologicalTreeFactory::createExtremumResidualTreeByCTs(image, adjacency);

    requireExactReconstruction(extremum, image, "extremum residual reconstruction");
    requireResidualSemantics(extremum.topology(), adjacency, MorphologicalTreeKind::ExtremumResidualTree, "extremum residual semantics");

    const auto extremumSaliency = HierarchySaliencyMap::computeTopologicalLevelEdgeMap(extremum.topology());
    require(!extremumSaliency.empty(), "extremum residual topological saliency");
    requireThrows<std::invalid_argument>([&] { static_cast<void>(HierarchySaliencyMap::computeNormalizedAltitudeEdgeMap(extremum)); },
                                         "residual gray-level altitude is not a hierarchy saliency valuation");
}

void testContrastInversion() {
    const auto image = makeImage(1, 6, {0, 2, 0, 2, 1, 3});
    const auto inverted = makeImage(1, 6, {3, 1, 3, 1, 2, 0});
    const RegularGridAdjacency2D adjacency(1, 6, 1.0);

    const auto requireConjugateEvolution = [&](const auto& primary, const auto& conjugate, const std::string& label) {
        const MorphologicalTree& lhs = primary.topology();
        const MorphologicalTree& rhs = conjugate.topology();
        requireEqual(lhs.numNodes(), rhs.numNodes(), label + " node count");
        std::size_t previousSupportCardinality = 0;
        for (NodeId node : lhs.aliveNodeIds()) {
            requireEqual(lhs.parent(node), rhs.parent(node), label + " parent");
            requireEqual(static_cast<int>(primary.nodeAltitude(node)) + static_cast<int>(conjugate.nodeAltitude(node)), 3, label + " complemented altitude");
            const auto lhsSupportRange = lhs.nodeSupport(node);
            const auto rhsSupportRange = rhs.nodeSupport(node);
            const std::vector<NodeId> lhsSupport(lhsSupportRange.begin(), lhsSupportRange.end());
            const std::vector<NodeId> rhsSupport(rhsSupportRange.begin(), rhsSupportRange.end());
            requireVectorEqual(lhsSupport, rhsSupport, label + " selected support");
            if (node != lhs.root()) {
                requireEqual(primary.nodeResidue(node), -conjugate.nodeResidue(node), label + " signed residue");
                const auto supportCardinality = lhsSupport.size();
                require(supportCardinality >= previousSupportCardinality, label + " nondecreasing support cardinality");
                previousSupportCardinality = supportCardinality;
            }
        }
        for (PixelId pixel = 0; pixel < image->getSize(); ++pixel) {
            requireEqual(lhs.smallestNode(pixel), rhs.smallestNode(pixel), label + " smallest node");
        }
    };

    requireConjugateEvolution(MorphologicalTreeFactory::createExtremumResidualTreeByCTs(image, adjacency),
                              MorphologicalTreeFactory::createExtremumResidualTreeByCTs(inverted, adjacency), "extremum contrast");


}

void testExtremumOrder() {
    using Access = BuilderAccess<int>;
    using Polarity = Access::Polarity;
    const Access::ExtremumLess rowMajorOrder;
    const Access::Extremum first{2, 0, Polarity::Maximum, 1};
    const Access::Extremum second{2, 3, Polarity::Minimum, 2};
    const Access::Extremum sameKey{2, 0, Polarity::Minimum, 3};
    require(rowMajorOrder(first, second), "row-major residual order uses the spatial minimum");
    require(!rowMajorOrder(first, first), "equal residual keys are equivalent");
    require(!rowMajorOrder(first, sameKey) && !rowMajorOrder(sameKey, first), "polarity and node identity do not affect selection priority");

    require(rowMajorOrder({1, 3, Polarity::Maximum, 1}, {2, 0, Polarity::Minimum, 2}),
            "support cardinality takes precedence over spatial order");
    require(!rowMajorOrder(second, first), "row-major residual order is asymmetric");
}

void testResidualTreeAssemblerInIsolation() {
    using Polarity = BuilderAccess<int>::Polarity;
    const auto image = makeTypedImage<int>(1, 4, {5, 5, 3, 1});
    const RegularGridAdjacency2D adjacency(1, 4, 1.0);
    auto minTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto maxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
    detail::DualMinMaxTreeIncrementalFilterLeaf<int> adjustment(&minTree, &maxTree, adjacency);
    TestAssembler<int> assembler(4, maxTree.topology().numInternalNodeSlots(), minTree.topology().numInternalNodeSlots());
    BuilderAccess<int>::Context context;

    const std::array<PixelId, 2> firstSupport{0, 1};
    NodeId leaf = maxTree.topology().smallestNode(0);
    NodeId parent = maxTree.topology().parent(leaf);
    BuilderAccess<int>::prepareLeaf(context, leaf, maxTree, minTree);
    requireEqual(assembler.recordSupport(firstSupport, Polarity::Maximum, 5, leaf, context.dualProperPartNode), NodeId{1}, "assembler first recorded support id");
    requireVectorEqual(assembler.pendingNodes(), std::vector<NodeId>{1}, "Xi after recording");
    adjustment.pruneMaxLeafAndUpdateMinTree(leaf, context.supportPixels, context.dualProperPartNode, context.wholeSupportSmallestNode);
    assembler.updateAssociations(minTree, maxTree, parent, adjustment.getLastSelectedSupportNode(), adjustment.getLastProperPartTransfers());
    for (Polarity polarity : {Polarity::Maximum, Polarity::Minimum}) {
        const auto& tree = polarity == Polarity::Maximum ? maxTree : minTree;
        const auto& associated = assembler.associatedNodes(polarity, tree.topology().smallestNode(0));
        requireVectorEqual(std::vector<NodeId>(associated.begin(), associated.end()), std::vector<NodeId>{1}, "explicit association after leveling");
    }
    const std::array<PixelId, 3> secondSupport{0, 1, 2};
    leaf = maxTree.topology().smallestNode(0);
    parent = maxTree.topology().parent(leaf);
    BuilderAccess<int>::prepareLeaf(context, leaf, maxTree, minTree);
    requireEqual(assembler.recordSupport(secondSupport, Polarity::Maximum, 3, leaf, context.dualProperPartNode), NodeId{2}, "assembler second recorded support id");
    requireVectorEqual(assembler.pendingNodes(), std::vector<NodeId>{2},
                       "Xi removes attached child");
    adjustment.pruneMaxLeafAndUpdateMinTree(leaf, context.supportPixels, context.dualProperPartNode, context.wholeSupportSmallestNode);
    assembler.updateAssociations(minTree, maxTree, parent, adjustment.getLastSelectedSupportNode(), adjustment.getLastProperPartTransfers());
    requireEqual(assembler.numRecordedSupports(), std::size_t{2}, "assembler recorded support count");

    auto output = std::move(assembler).finalize(image, adjacency, 1);
    requireVectorEqual(output.parent, std::vector<NodeId>{0, 2, 0}, "assembler parent materialization");
    requireVectorEqual(output.smallestNodeMap, std::vector<NodeId>{1, 1, 2, 0}, "assembler smallest-node-map materialization");
    requireVectorEqual(output.nodeAltitudes, std::vector<int>{1, 5, 3}, "assembler altitude materialization");


}

// Independently check Xi and Eq. (31) against complete recorded pixel sets
// after every leveling, including pending supports outside the selected leaf.
template <AltitudeValue T> void checkResidualAssociations(const ImagePtr<T>& image, const RegularGridAdjacency2D& adjacency) {
    using Polarity = typename BuilderAccess<T>::Polarity;
    auto minTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto maxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
    detail::DualMinMaxTreeIncrementalFilterLeaf<T> adjustment(&minTree, &maxTree, adjacency);
    TestAssembler<T> assembler(image->getSize(), maxTree.topology().numInternalNodeSlots(), minTree.topology().numInternalNodeSlots());
    typename BuilderAccess<T>::Context context;
    std::vector<std::vector<PixelId>> supports(1);
    std::vector<bool> pending{false};
    std::vector<NodeId> expectedParents{0};
    std::size_t previousArea = 0;
    const auto verify = [&] {
        std::vector<NodeId> expectedXi;
        for (NodeId node = 1; node < static_cast<NodeId>(pending.size()); ++node)
            if (pending[node])
                expectedXi.push_back(node);
        auto actualXi = assembler.pendingNodes();
        std::sort(actualXi.begin(), actualXi.end());
        requireVectorEqual(actualXi, expectedXi, "Xi equals maximal recorded supports");
        for (Polarity polarity : {Polarity::Maximum, Polarity::Minimum}) {
            const auto& topology = (polarity == Polarity::Maximum ? maxTree : minTree).topology();
            std::vector<std::vector<NodeId>> expected(topology.numInternalNodeSlots());
            for (NodeId node : expectedXi) {
                const NodeId owner = topology.smallestNode(supports[node].front());
                require(topology.isAlive(owner), "pending support has a live component owner");
                for (PixelId pixel : supports[node])
                    requireEqual(topology.smallestNode(pixel), owner, "complete support satisfies the association relation");
                expected[owner].push_back(node);
            }
            for (NodeId owner = 0; owner < topology.numInternalNodeSlots(); ++owner) {
                const auto& associated = assembler.associatedNodes(polarity, owner);
                auto actual = std::vector<NodeId>(associated.begin(), associated.end());
                std::sort(actual.begin(), actual.end());
                requireVectorEqual(actual, expected[owner], "explicit relation equals full-support containment");
            }
        }
    };
    verify();
    for (;;) {
        NodeId selected = InvalidNode;
        Polarity polarity = Polarity::Maximum;
        std::pair<std::size_t, PixelId> best;
        // Literal minimum over the current leaves, independent of the agenda.
        for (Polarity side : {Polarity::Maximum, Polarity::Minimum}) {
            const auto& topology = (side == Polarity::Maximum ? maxTree : minTree).topology();
            for (NodeId node : topology.aliveNodeIds()) {
                if (node == topology.root() || !topology.isLeaf(node))
                    continue;
                const std::pair<std::size_t, PixelId> key{static_cast<std::size_t>(topology.properPartCardinality(node)),
                                                       *std::min_element(topology.properPart(node).begin(), topology.properPart(node).end())};
                if (selected == InvalidNode || key < best) {
                    selected = node;
                    polarity = side;
                    best = key;
                }
            }
        }
        if (selected == InvalidNode)
            break;
        auto& primal = polarity == Polarity::Maximum ? maxTree : minTree;
        auto& dual = polarity == Polarity::Maximum ? minTree : maxTree;
        BuilderAccess<T>::prepareLeaf(context, selected, primal, dual);
        for (PixelId pixel : context.supportPixels)
            requireEqual(dual.topology().smallestNode(pixel), context.dualProperPartNode, "complete flat zone lies in one dual proper part");
        require(best.first >= previousArea, "selected support cardinalities are non-decreasing");
        previousArea = best.first;
        const T nodeAltitude = primal.nodeAltitude(selected);
        const NodeId primalParent = primal.topology().parent(selected);
        const auto signedResidualValue = static_cast<AltitudeDifference<T>>(nodeAltitude) -
                                         static_cast<AltitudeDifference<T>>(primal.nodeAltitude(primalParent));
        require(polarity == Polarity::Maximum ? signedResidualValue > 0 : signedResidualValue < 0,
                "residual sign agrees with the selected polarity");
        for (NodeId node = 1; node < static_cast<NodeId>(pending.size()); ++node) {
            if (pending[node] && std::all_of(supports[node].begin(), supports[node].end(), [&](PixelId pixel) {
                    return std::find(context.supportPixels.begin(), context.supportPixels.end(), pixel) != context.supportPixels.end();
                })) {
                pending[node] = false;
                expectedParents[node] = static_cast<NodeId>(supports.size());
            }
        }
        supports.emplace_back(context.supportPixels);
        pending.push_back(true);
        expectedParents.push_back(0);
        static_cast<void>(assembler.recordSupport(context.supportPixels, polarity, nodeAltitude, selected, context.dualProperPartNode));
        if (polarity == Polarity::Maximum)
            adjustment.pruneMaxLeafAndUpdateMinTree(selected, context.supportPixels, context.dualProperPartNode, context.wholeSupportSmallestNode);
        else
            adjustment.pruneMinLeafAndUpdateMaxTree(selected, context.supportPixels, context.dualProperPartNode, context.wholeSupportSmallestNode);
        for (PixelId pixel : context.supportPixels) {
            requireEqual(primal.topology().smallestNode(pixel), primalParent, "primal pruning establishes containment over the complete support");
            requireEqual(dual.topology().smallestNode(pixel), adjustment.getLastSelectedSupportNode(),
                         "dual transfer destination contains the complete selected support");
        }
        assembler.updateAssociations(minTree, maxTree, primalParent, adjustment.getLastSelectedSupportNode(), adjustment.getLastProperPartTransfers());
        verify();
    }
    const auto output = std::move(assembler).finalize(image, adjacency, minTree.nodeAltitude(minTree.topology().root()));
    requireVectorEqual(output.parent, expectedParents, "finalization attaches implicit Xi to the root and preserves assigned parents");
    require(output.parent.size() <= static_cast<std::size_t>(image->getSize()), "nonempty proper parts bound the residual node count by pixels");
    std::vector<int> properPartCounts(output.parent.size(), 0);
    for (PixelId pixel = 0; pixel < image->getSize(); ++pixel) {
        const NodeId owner = output.smallestNodeMap[static_cast<std::size_t>(pixel)];
        require(owner >= 0 && static_cast<std::size_t>(owner) < output.parent.size(), "every residual pixel has a valid owner");
        ++properPartCounts[static_cast<std::size_t>(owner)];
        requireEqual(output.nodeAltitudes[static_cast<std::size_t>(owner)], (*image)[pixel], "residual proper-part valuations reconstruct the input");
    }
    for (int count : properPartCounts) {
        require(count > 0, "every residual proper part is nonempty, including the terminal root");
    }
}

void testResidualAssociationsAgainstFullSupports() {
    checkResidualAssociations(makeTypedImage<std::uint16_t>(1, 2, {1, 7}), RegularGridAdjacency2D(1, 2, 1.0));
    checkResidualAssociations(makeTypedImage<int>(1, 5, {3, 1, 2, 1, 3}), RegularGridAdjacency2D(1, 5, 1.0));
    for (int code = 0; code < 729; ++code)
        checkResidualAssociations(makeRadixImage(code), RegularGridAdjacency2D(2, 3, 1.0));
    std::mt19937 random(20261004);
    for (int trial = 0; trial < 64; ++trial) {
        auto image = Image<std::int32_t>::create(4, 4);
        for (PixelId pixel = 0; pixel < 16; ++pixel)
            (*image)[pixel] = static_cast<std::int32_t>(random() % 17) - 8;
        checkResidualAssociations(image, RegularGridAdjacency2D(4, 4, trial % 3 == 0 ? 3.2 : trial % 3 == 1 ? 1.5 : 1.0));
    }
}

void testResidualLeafPreparation() {
    const auto image = makeImage(1, 4, {2, 1, 1, 0});
    const RegularGridAdjacency2D adjacency(1, 4, 1.0);
    auto minTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto maxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
    detail::DualMinMaxTreeIncrementalFilterLeaf<std::uint8_t> adjustment(&minTree, &maxTree, adjacency);
    BuilderAccess<std::uint8_t>::Context context;
    const NodeId leaf = maxTree.topology().smallestNode(0);
    const NodeId parent = maxTree.topology().parent(leaf);
    BuilderAccess<std::uint8_t>::prepareLeaf(context, leaf, maxTree, minTree);
    requireVectorEqual(context.supportPixels, std::vector<PixelId>{0}, "initial leaf support");
    adjustment.pruneMaxLeafAndUpdateMinTree(leaf, context.supportPixels, context.dualProperPartNode, context.wholeSupportSmallestNode);
    requireRebuiltComponentTrees(minTree, maxTree, makeImage(1, 4, {1, 1, 1, 0}), adjacency,
                                "updater discovers the initial leaf boundary");
    require(maxTree.topology().isLeaf(parent), "primal contraction creates the next extremal leaf");
    BuilderAccess<std::uint8_t>::prepareLeaf(context, parent, maxTree, minTree);
    std::ranges::sort(context.supportPixels);
    requireVectorEqual(context.supportPixels, std::vector<PixelId>{0, 1, 2}, "merged leaf support from proper parts");
    adjustment.pruneMaxLeafAndUpdateMinTree(parent, context.supportPixels, context.dualProperPartNode, context.wholeSupportSmallestNode);
    requireRebuiltComponentTrees(minTree, maxTree, makeImage(1, 4, {0, 0, 0, 0}), adjacency,
                                "updater discovers the merged leaf boundary");

    // A constant domain has no nonterminal extremum; the builder never prepares its root.
    const auto constantImage = makeImage(1, 4, {7, 7, 7, 7});
    const auto constant = MorphologicalTreeFactory::createExtremumResidualTreeByCTs(constantImage, adjacency);
    requireEqual(constant.topology().numNodes(), 1, "constant evolution records only the terminal root");
    requireExactReconstruction(constant, constantImage, "constant evolution preserves its input valuation");
}

// The builder supplies the complete workset; verify both polarities and both bucket backends.
template <AltitudeValue T> void testCertifiedLeafUpdate() {
    for (bool pruneMaximum : {false, true}) {
        const auto image = makeTypedImage<T>(2, 3, {T{0}, T{2}, T{1}, T{2}, T{0}, T{0}});
        const RegularGridAdjacency2D adjacency(2, 3, 1.0);
        auto minTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
        auto maxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
        auto& primal = pruneMaximum ? maxTree : minTree;
        auto& dual = pruneMaximum ? minTree : maxTree;
        auto reference = pruneMaximum ? MorphologicalTreeFactory::createMaxTree(image, adjacency)
                                      : MorphologicalTreeFactory::createMinTree(image, adjacency);
        const NodeId leaf = primal.topology().leaves().front();
        const T targetAltitude = primal.nodeAltitude(primal.topology().parent(leaf));
        typename BuilderAccess<T>::Context context;
        BuilderAccess<T>::prepareLeaf(context, leaf, primal, dual);
        const auto expectedImage = makeTypedImage<T>(2, 3, {T{0}, T{2}, T{1}, T{2}, T{0}, T{0}});
        for (PixelId pixel : context.supportPixels) (*expectedImage)[pixel] = targetAltitude;
        detail::DualMinMaxTreeIncrementalFilterLeaf<T> adjustment(&minTree, &maxTree, adjacency);
        if (pruneMaximum)
            adjustment.pruneMaxLeafAndUpdateMinTree(leaf, context.supportPixels, context.dualProperPartNode,
                                                   context.wholeSupportSmallestNode);
        else
            adjustment.pruneMinLeafAndUpdateMaxTree(leaf, context.supportPixels, context.dualProperPartNode,
                                                   context.wholeSupportSmallestNode);
        reference.pruneNode(leaf);
        requireVectorEqual(parentBuffer(primal), parentBuffer(reference), "internal leaf contraction preserves public prune topology");
        requireEqual(primal.topology().getMutationVersion(), reference.topology().getMutationVersion(),
                     "internal leaf contraction preserves mutation versions");
        requireRebuiltComponentTrees(minTree, maxTree, expectedImage, adjacency, "certified leaf update");
    }
}

void testLeafBoundaryExcludesNonAdjacentDualBranches() {
    const RegularGridAdjacency2D adjacency(1, 4, 1.0);
    for (bool invertContrast : {false, true}) {
        const auto image = invertContrast ? makeImage(1, 4, {0, 2, 1, 2}) : makeImage(1, 4, {3, 1, 2, 1});
        auto minTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
        auto maxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
        auto& primal = invertContrast ? minTree : maxTree;
        auto& dual = invertContrast ? maxTree : minTree;
        detail::DualMinMaxTreeIncrementalFilterLeaf<std::uint8_t> adjustment(&minTree, &maxTree, adjacency);
        BuilderAccess<std::uint8_t>::Context context;
        const NodeId leaf = primal.topology().smallestNode(0);
        const NodeId adjacentNode = dual.topology().smallestNode(1);
        const NodeId nonAdjacentNode = dual.topology().smallestNode(3);
        require(adjacentNode != nonAdjacentNode, "equal-level dual branches are distinct before leveling");
        BuilderAccess<std::uint8_t>::prepareLeaf(context, leaf, primal, dual);
        requireEqual(context.wholeSupportSmallestNode, context.dualProperPartNode, "whole-part certificate remains available");
        if (invertContrast)
            adjustment.pruneMinLeafAndUpdateMaxTree(leaf, context.supportPixels, context.dualProperPartNode, context.wholeSupportSmallestNode);
        else
            adjustment.pruneMaxLeafAndUpdateMinTree(leaf, context.supportPixels, context.dualProperPartNode, context.wholeSupportSmallestNode);
        require(dual.topology().smallestNode(0) != dual.topology().smallestNode(3), "non-adjacent equal-level regions remain separate");
        const auto leveledImage = invertContrast ? makeImage(1, 4, {2, 2, 1, 2}) : makeImage(1, 4, {1, 1, 2, 1});
        requireRebuiltComponentTrees(minTree, maxTree, leveledImage, adjacency, "spatial boundary excludes the non-adjacent branch");
    }
}

void testPartialDualProperPartPreservesReattachments() {
    const RegularGridAdjacency2D adjacency(1, 4, 1.0);
    for (bool invertContrast : {false, true}) {
        const auto image = invertContrast ? makeImage(1, 4, {1, 0, 2, 0}) : makeImage(1, 4, {1, 2, 0, 2});
        auto minTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
        auto maxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
        auto& primal = invertContrast ? minTree : maxTree;
        auto& dual = invertContrast ? maxTree : minTree;
        const NodeId leaf = primal.topology().smallestNode(1);
        BuilderAccess<std::uint8_t>::Context context;
        BuilderAccess<std::uint8_t>::prepareLeaf(context, leaf, primal, dual);
        requireVectorEqual(context.supportPixels, std::vector<PixelId>{1}, "selected leaf support is one pixel");
        requireEqual(context.dualProperPartNode, dual.topology().smallestNode(3), "dual proper part also contains an unselected extremum");
        requireEqual(context.wholeSupportSmallestNode, InvalidNode, "selected support is a strict subset of the dual proper part");

        detail::DualMinMaxTreeIncrementalFilterLeaf<std::uint8_t> adjustment(&minTree, &maxTree, adjacency);
        if (invertContrast) {
            adjustment.pruneMinLeafAndUpdateMaxTree(leaf, context.supportPixels, context.dualProperPartNode, context.wholeSupportSmallestNode);
        } else {
            adjustment.pruneMaxLeafAndUpdateMinTree(leaf, context.supportPixels, context.dualProperPartNode, context.wholeSupportSmallestNode);
        }
        const auto leveledImage = invertContrast ? makeImage(1, 4, {1, 1, 2, 0}) : makeImage(1, 4, {1, 1, 0, 2});
        requireRebuiltComponentTrees(minTree, maxTree, leveledImage, adjacency, "partial dual proper part retains the required reattachments");
    }
}

template <class T> void testTypedResidualConstruction(const std::string& label) {
    const auto image = makeTypedImage<T>(2, 3, {T{-3}, T{2}, T{1}, T{4}, T{-1}, T{0}});
    const RegularGridAdjacency2D adjacency(2, 3, 1.0);
    const auto extremum = MorphologicalTreeFactory::createExtremumResidualTreeByCTs(image, adjacency);
    requireExactReconstruction(extremum, image, label + " extremum reconstruction");
}

void testAltitudeTypesAndDegenerateDomain() {
    testTypedResidualConstruction<std::int16_t>("int16");
    testTypedResidualConstruction<std::int32_t>("int32");
    testTypedResidualConstruction<float>("float");
    testTypedResidualConstruction<double>("double");

    const auto unsignedImage = makeTypedImage<std::uint16_t>(2, 3, {0, 5, 2, 7, 1, 3});
    const RegularGridAdjacency2D unsignedAdjacency(2, 3, 1.0);
    requireExactReconstruction(MorphologicalTreeFactory::createExtremumResidualTreeByCTs(unsignedImage, unsignedAdjacency), unsignedImage,
                               "uint16 extremum reconstruction");

    const auto singleton = makeTypedImage<double>(1, 1, {2.5});
    const auto extremumSingleton = MorphologicalTreeFactory::createExtremumResidualTreeByCTs(singleton, 1.0);
    requireEqual(extremumSingleton.topology().numNodes(), 1, "extremum singleton node count");
    requireExactReconstruction(extremumSingleton, singleton, "extremum singleton reconstruction");
}

void testBuilderAndInputContracts() {
    requireThrows<std::invalid_argument>(
        [] { static_cast<void>(MorphologicalTreeFactory::createExtremumResidualTreeByCTs(ImageUInt8Ptr{}, 1.0)); },
        "radius factory rejects a null image before reading its dimensions");
    const auto image = makeImage(2, 3, {0, 2, 1, 2, 0, 0});
    const RegularGridAdjacency2D adjacency(2, 3, 1.0);

    const auto factoryTree = MorphologicalTreeFactory::createExtremumResidualTreeByCTs(image, adjacency);
    ExtremumResidualTreeByCTsBuilder<std::uint8_t> extremumBuilder(adjacency);
    const auto hierarchy = extremumBuilder.build(image);
    requireVectorEqual(parentBuffer(factoryTree), hierarchy.parent, "factory and builder parents");
    requireVectorEqual(smallestNodeMapBuffer(factoryTree), hierarchy.smallestNodeMap,
                       "factory and builder smallest-node map");
    requireVectorEqual(altitudeBuffer(factoryTree), hierarchy.nodeAltitudes, "factory and builder altitudes");

    requireThrows<std::invalid_argument>([&] { static_cast<void>(extremumBuilder.build(ImageUInt8Ptr{})); },
                                         "direct residual builder rejects a null image");
    requireVectorEqual(parentBuffer(factoryTree), hierarchy.parent, "failed build preserves the previously returned hierarchy");
    const auto rebuilt = extremumBuilder.build(image);
    requireVectorEqual(hierarchy.parent, rebuilt.parent, "builder reuse after failure preserves parents");
    requireVectorEqual(hierarchy.smallestNodeMap, rebuilt.smallestNodeMap, "builder reuse after failure preserves smallest-node map");
    requireVectorEqual(hierarchy.nodeAltitudes, rebuilt.nodeAltitudes, "builder reuse after failure preserves altitudes");

    const RegularGridAdjacency2D mismatched(1, 6, 1.0);
    ExtremumResidualTreeByCTsBuilder<std::uint8_t> mismatchedBuilder(mismatched);
    requireThrows<std::invalid_argument>([&] { static_cast<void>(mismatchedBuilder.build(image)); }, "residual builder rejects a mismatched domain");

    const std::array<GridOffset2D, 1> originOnly{{{0, 0}}};
    const auto disconnected = RegularGridAdjacency2D::fromStructuringElement(2, 3, originOnly);
    ExtremumResidualTreeByCTsBuilder<std::uint8_t> disconnectedBuilder(disconnected);
    requireThrows<std::invalid_argument>([&] { static_cast<void>(disconnectedBuilder.build(image)); },
                                         "residual builder rejects a disconnected adjacency");

    if constexpr (contract::validationsEnabled) {
        auto nonFinite = makeTypedImage<float>(1, 2, {0.0F, std::numeric_limits<float>::quiet_NaN()});
        requireThrows<std::invalid_argument>([&] { static_cast<void>(MorphologicalTreeFactory::createExtremumResidualTreeByCTs(nonFinite, 1.0)); },
                                             "extremum factory rejects NaN");
        (*nonFinite)[1] = std::numeric_limits<float>::infinity();
        requireThrows<std::invalid_argument>([&] { static_cast<void>(MorphologicalTreeFactory::createExtremumResidualTreeByCTs(nonFinite, 1.0)); },
                                             "extremum factory rejects infinity");
    }
}

void testSharedArbitrarySymmetricAdjacency() {
    const std::array<GridOffset2D, 7> offsets{{
        {0, 0},
        {-1, 0},
        {1, 0},
        {0, -1},
        {0, 1},
        {0, -2},
        {0, 2},
    }};
    const RegularGridAdjacency2D adjacency = RegularGridAdjacency2D::fromStructuringElement(4, 5, offsets);
    const auto image = makeImage(4, 5,
                                 {
                                     0, 3, 1, 4, 2, 2, 1, 4, 0, 3, 4, 0, 2, 3, 1, 1, 4, 0, 2, 3,
                                 });
    const auto extremum = MorphologicalTreeFactory::createExtremumResidualTreeByCTs(image, adjacency);

    requireExactReconstruction(extremum, image, "custom-adjacency extremum reconstruction");
    requireResidualSemantics(extremum.topology(), adjacency, MorphologicalTreeKind::ExtremumResidualTree, "custom-adjacency extremum semantics");
}

void testOwnedBuildResultLifetime() {
    const NodeId* parentData = nullptr;
    const NodeId* smallestNodeData = nullptr;
    const std::uint8_t* altitudeData = nullptr;
    auto hierarchy = [&] {
        auto image = makeImage(2, 3, {0, 2, 1, 2, 0, 0});
        RegularGridAdjacency2D adjacency(2, 3, 1.0);
        ExtremumResidualTreeByCTsBuilder<std::uint8_t> builder(std::move(adjacency));
        auto first = builder.build(image);
        parentData = first.parent.data();
        smallestNodeData = first.smallestNodeMap.data();
        altitudeData = first.nodeAltitudes.data();
        std::fill_n(image->rawData(), image->getSize(), std::uint8_t{7});
        const auto second = builder.build(image);
        requireEqual(second.parent.size(), std::size_t{1}, "reused builder processes the new constant image");
        requireEqual(second.nodeAltitudes.front(), std::uint8_t{7}, "reused builder records the new terminal altitude");
        return first;
    }();
    require(hierarchy.parent.data() == parentData && hierarchy.smallestNodeMap.data() == smallestNodeData &&
                hierarchy.nodeAltitudes.data() == altitudeData,
            "returned native buffers are moved without copies");
    const NativeHierarchyView<std::uint8_t> view{hierarchy.parent, hierarchy.smallestNodeMap, hierarchy.nodeAltitudes,
                                               hierarchy.root, hierarchy.gridDomain2D, hierarchy.semantics};
    detail::NativeHierarchyValidation::validateComplete(view);
    const auto tree = detail::materializeNativeHierarchy(std::move(hierarchy));
    const auto original = makeImage(2, 3, {0, 2, 1, 2, 0, 0});
    requireExactReconstruction(tree, original, "result survives image modification, builder reuse, and destruction of construction state");
    requireResidualSemantics(tree.topology(), RegularGridAdjacency2D(2, 3, 1.0), MorphologicalTreeKind::ExtremumResidualTree,
                             "returned hierarchy owns its adjacency");
}

} // namespace

int main() {
    testFactoryAndSaliencyIntegration();
    testContrastInversion();
    testExtremumOrder();
    testResidualTreeAssemblerInIsolation();
    testResidualAssociationsAgainstFullSupports();
    testResidualLeafPreparation();
    testLeafBoundaryExcludesNonAdjacentDualBranches();
    testPartialDualProperPartPreservesReattachments();
    testCertifiedLeafUpdate<std::uint8_t>();
    testCertifiedLeafUpdate<std::int32_t>();
    testCertifiedLeafUpdate<double>();
    testAltitudeTypesAndDegenerateDomain();
    testBuilderAndInputContracts();
    testSharedArbitrarySymmetricAdjacency();
    testOwnedBuildResultLifetime();
    return 0;
}
