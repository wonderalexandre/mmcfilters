#include "support/TestSupport.hpp"

#include "mmcfilters/trees/saliency/HierarchySaliencyMap.hpp"
#include "mmcfilters/trees/sdrt/ExtremumResidualTreeByCTsBuilder.hpp"

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
namespace mmcfilters::sdrt::detail {
template <AltitudeValue T> struct ExtremumResidualTreeTestAccess {
    using Builder = ExtremumResidualTreeByCTsBuilder<T>;
    using Context = typename Builder::LeafContext;

    static void prepareLeaf(Context& context, NodeId leaf, const ValuedMorphologicalTree<T>& primal,
                            const ValuedMorphologicalTree<T>& dual) {
        Builder::prepareResidualLeaf(context, leaf, primal, dual);
    }
    struct Assembler {
        typename Builder::AssemblyState state;
        Assembler(std::size_t pixels, std::size_t maxSlots, std::size_t minSlots) : state(pixels, maxSlots, minSlots) {}
        NodeId recordSupport(std::span<const PixelId> support, Polarity polarity, T nodeAltitude, NodeId primalLeaf, NodeId dualProperPart) {
            return Builder::recordSupport(state, support, polarity, nodeAltitude, primalLeaf, dualProperPart);
        }
        std::size_t numRecordedSupports() const { return state.parents_.size() - 1; }
        std::vector<NodeId> pendingNodes() const {
            std::vector<NodeId> pending;
            for (NodeId node = 1; node < static_cast<NodeId>(state.parents_.size()); ++node)
                if (state.parents_[static_cast<std::size_t>(node)] == InvalidNode)
                    pending.push_back(node);
            return pending;
        }
        const auto& associatedNodes(Polarity polarity, NodeId componentNode) const {
            return (polarity == Polarity::Maximum ? state.maxAssociations_ : state.minAssociations_).at(componentNode);
        }
        template <class Transfer>
        void updateAssociations(const ValuedMorphologicalTree<T>& minTree, const ValuedMorphologicalTree<T>& maxTree,
                                NodeId primal, NodeId dual, std::span<const Transfer> transfers) {
            Builder::updateAssociations(state, minTree, maxTree, primal, dual, transfers);
        }
        auto finalize(const ImagePtr<T>& image, const RegularGridAdjacency2D& adjacency, T terminalAltitude) && {
            return Builder::finalize(std::move(state), image,
                                    makeMorphologicalTreeSemantics(MorphologicalTreeKind::ExtremumResidualTree, SharedAdjacencyContext{adjacency}),
                                    terminalAltitude);
        }
    };
};
} // namespace mmcfilters::sdrt::detail

template <class T> using BuilderAccess = mmcfilters::sdrt::detail::ExtremumResidualTreeTestAccess<T>;
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

static_assert(!HasInfinityPixelAccessor<sdrt::ExtremumResidualTreeByCTsBuilder<std::uint8_t>>, "extremum construction must not expose infinity-pixel state");
static_assert(!LvalueFinalizableAssembler<TestAssembler<int>>, "residual support assembly must not be finalizable through an lvalue");
static_assert(RvalueFinalizableAssembler<TestAssembler<int>>, "residual support assembly must be finalizable through an rvalue");
static_assert(!std::is_constructible_v<sdrt::ExtremumResidualTreeByCTsBuilder<std::uint8_t>, RegularGridAdjacency2D, sdrt::SpatialOrder>,
              "extremum builders must require their mode-specific option object");

template <class T> std::vector<T> toVector(std::span<const T> values) { return {values.begin(), values.end()}; }

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

    sdrt::ExtremumResidualTreeOptions customOrderOptions;
    customOrderOptions.spatialOrder = sdrt::SpatialOrder(std::vector<PixelId>{5, 4, 3, 2, 1, 0});
    requireConjugateEvolution(MorphologicalTreeFactory::createExtremumResidualTreeByCTs(image, adjacency, customOrderOptions),
                              MorphologicalTreeFactory::createExtremumResidualTreeByCTs(inverted, adjacency, customOrderOptions),
                              "custom-order extremum contrast");
}

void testSelfDualResidualSchedule() {
    const sdrt::SelfDualResidualOrder rowMajorOrder;
    const sdrt::SelfDualResidualKey first{2, 0};
    const sdrt::SelfDualResidualKey second{2, 3};
    require(rowMajorOrder.compareExtremumKeys(first, second), "row-major residual order uses the spatial minimum");
    require(!rowMajorOrder.compareExtremumKeys(first, first), "equal residual keys are equivalent");

    const sdrt::SpatialOrder reverseOrder(std::vector<PixelId>{3, 2, 1, 0});
    const sdrt::SelfDualResidualSchedule schedule(reverseOrder);
    const std::array<sdrt::SelfDualResidualKey, 3> keys{{{2, 1}, {1, 2}, {1, 3}}};
    requireEqual(schedule.selectSmallestExtremumIndex(keys), std::size_t{2}, "custom spatial order selects the least equal-cardinality support");
    const std::array<sdrt::SelfDualResidualKey, 2> duplicateKeys{{{1, 2}, {1, 2}}};
    requireThrows<std::invalid_argument>([&] { static_cast<void>(schedule.selectSmallestExtremumIndex(duplicateKeys)); },
                                         "self-dual residual schedule rejects duplicate keys");
    requireEqual(reverseOrder.spatialMinimum(std::array<PixelId, 3>{0, 2, 3}), PixelId{3}, "custom spatial minimum");
    requireThrows<std::invalid_argument>([] { static_cast<void>(sdrt::SpatialOrder(std::vector<PixelId>{0, 0})); }, "spatial order rejects duplicate pixels");
}

void testResidualTreeAssemblerInIsolation() {
    const auto image = makeTypedImage<int>(1, 4, {5, 5, 3, 1});
    const RegularGridAdjacency2D adjacency(1, 4, 1.0);
    auto minTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto maxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
    adjust::DualMinMaxTreeIncrementalFilterLeaf<int> adjustment(&minTree, &maxTree, adjacency);
    TestAssembler<int> assembler(4, maxTree.topology().numInternalNodeSlots(), minTree.topology().numInternalNodeSlots());
    BuilderAccess<int>::Context context;

    const std::array<PixelId, 2> firstSupport{0, 1};
    NodeId leaf = maxTree.topology().smallestNode(0);
    NodeId parent = maxTree.topology().parent(leaf);
    BuilderAccess<int>::prepareLeaf(context, leaf, maxTree, minTree);
    requireEqual(assembler.recordSupport(firstSupport, sdrt::Polarity::Maximum, 5, leaf, context.dualProperPartNode), NodeId{1}, "assembler first recorded support id");
    requireVectorEqual(assembler.pendingNodes(), std::vector<NodeId>{1}, "Xi after recording");
    adjustment.pruneMaxLeafAndUpdateMinTree(leaf, context.supportPixels, context.dualProperPartNode, {},
                                            context.wholeSupportSmallestNode);
    assembler.updateAssociations(minTree, maxTree, parent, adjustment.getLastSelectedSupportNode(), adjustment.getLastProperPartTransfers());
    for (sdrt::Polarity polarity : {sdrt::Polarity::Maximum, sdrt::Polarity::Minimum}) {
        const auto& tree = polarity == sdrt::Polarity::Maximum ? maxTree : minTree;
        const auto& associated = assembler.associatedNodes(polarity, tree.topology().smallestNode(0));
        requireVectorEqual(std::vector<NodeId>(associated.begin(), associated.end()), std::vector<NodeId>{1}, "explicit association after leveling");
    }
    const std::array<PixelId, 3> secondSupport{0, 1, 2};
    leaf = maxTree.topology().smallestNode(0);
    parent = maxTree.topology().parent(leaf);
    BuilderAccess<int>::prepareLeaf(context, leaf, maxTree, minTree);
    requireEqual(assembler.recordSupport(secondSupport, sdrt::Polarity::Maximum, 3, leaf, context.dualProperPartNode), NodeId{2}, "assembler second recorded support id");
    requireVectorEqual(assembler.pendingNodes(), std::vector<NodeId>{2},
                       "Xi removes attached child");
    adjustment.pruneMaxLeafAndUpdateMinTree(leaf, context.supportPixels, context.dualProperPartNode, {},
                                            context.wholeSupportSmallestNode);
    assembler.updateAssociations(minTree, maxTree, parent, adjustment.getLastSelectedSupportNode(), adjustment.getLastProperPartTransfers());
    requireEqual(assembler.numRecordedSupports(), std::size_t{2}, "assembler recorded support count");

    auto output = std::move(assembler).finalize(image, adjacency, 1);
    requireVectorEqual(output.parents, std::vector<NodeId>{0, 2, 0}, "assembler parent materialization");
    requireVectorEqual(output.smallestNodeMap, std::vector<NodeId>{1, 1, 2, 0}, "assembler smallest-node-map materialization");
    requireVectorEqual(output.nodeAltitudes, std::vector<int>{1, 5, 3}, "assembler altitude materialization");

    TestAssembler<int> invalidAssembler(4, 1, 1);
    static_cast<void>(invalidAssembler.recordSupport(firstSupport, sdrt::Polarity::Maximum, 5, 0, 0));
    using Transfer = adjust::DualMinMaxTreeIncrementalFilterLeaf<int>::ProperPartTransfer;
    invalidAssembler.updateAssociations(minTree, maxTree, maxTree.topology().root(), minTree.topology().root(), std::span<const Transfer>{});
    requireThrows<std::runtime_error>(
        [&] { static_cast<void>(invalidAssembler.recordSupport(firstSupport, sdrt::Polarity::Maximum, 5, maxTree.topology().root(), minTree.topology().root())); },
        "assembler rejects an empty residual proper part");


}

// Independently check Xi and Eq. (31) against complete recorded pixel sets
// after every leveling, including pending supports outside the selected leaf.
template <AltitudeValue T> void checkResidualAssociations(const ImagePtr<T>& image, const RegularGridAdjacency2D& adjacency) {
    auto minTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto maxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
    adjust::DualMinMaxTreeIncrementalFilterLeaf<T> adjustment(&minTree, &maxTree, adjacency);
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
        for (sdrt::Polarity polarity : {sdrt::Polarity::Maximum, sdrt::Polarity::Minimum}) {
            const auto& topology = (polarity == sdrt::Polarity::Maximum ? maxTree : minTree).topology();
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
        sdrt::Polarity polarity = sdrt::Polarity::Maximum;
        sdrt::SelfDualResidualKey best;
        const sdrt::SelfDualResidualOrder order;
        // Literal minimum over the current leaves, independent of the agenda.
        for (sdrt::Polarity side : {sdrt::Polarity::Maximum, sdrt::Polarity::Minimum}) {
            const auto& topology = (side == sdrt::Polarity::Maximum ? maxTree : minTree).topology();
            for (NodeId node : topology.aliveNodeIds()) {
                if (node == topology.root() || !topology.isLeaf(node))
                    continue;
                const sdrt::SelfDualResidualKey key{static_cast<std::size_t>(topology.properPartCardinality(node)),
                                                    *std::min_element(topology.properPart(node).begin(), topology.properPart(node).end())};
                if (selected == InvalidNode || order(key, best)) {
                    selected = node;
                    polarity = side;
                    best = key;
                }
            }
        }
        if (selected == InvalidNode)
            break;
        auto& primal = polarity == sdrt::Polarity::Maximum ? maxTree : minTree;
        auto& dual = polarity == sdrt::Polarity::Maximum ? minTree : maxTree;
        BuilderAccess<T>::prepareLeaf(context, selected, primal, dual);
        for (PixelId pixel : context.supportPixels)
            requireEqual(dual.topology().smallestNode(pixel), context.dualProperPartNode, "complete flat zone lies in one dual proper part");
        require(best.supportCardinality >= previousArea, "selected support cardinalities are non-decreasing");
        previousArea = best.supportCardinality;
        const T nodeAltitude = primal.nodeAltitude(selected);
        const NodeId primalParent = primal.topology().parent(selected);
        const auto signedResidualValue = static_cast<AltitudeDifference<T>>(nodeAltitude) -
                                         static_cast<AltitudeDifference<T>>(primal.nodeAltitude(primalParent));
        require(polarity == sdrt::Polarity::Maximum ? signedResidualValue > 0 : signedResidualValue < 0,
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
        if (polarity == sdrt::Polarity::Maximum)
            adjustment.pruneMaxLeafAndUpdateMinTree(selected, context.supportPixels, context.dualProperPartNode, {},
                                                    context.wholeSupportSmallestNode);
        else
            adjustment.pruneMinLeafAndUpdateMaxTree(selected, context.supportPixels, context.dualProperPartNode, {},
                                                    context.wholeSupportSmallestNode);
        for (PixelId pixel : context.supportPixels) {
            requireEqual(primal.topology().smallestNode(pixel), primalParent, "primal pruning establishes containment over the complete support");
            requireEqual(dual.topology().smallestNode(pixel), adjustment.getLastSelectedSupportNode(),
                         "dual transfer destination contains the complete selected support");
        }
        assembler.updateAssociations(minTree, maxTree, primalParent, adjustment.getLastSelectedSupportNode(), adjustment.getLastProperPartTransfers());
        verify();
    }
    const auto output = std::move(assembler).finalize(image, adjacency, minTree.nodeAltitude(minTree.topology().root()));
    requireVectorEqual(output.parents, expectedParents, "finalization attaches implicit Xi to the root and preserves assigned parents");
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
    adjust::DualMinMaxTreeIncrementalFilterLeaf<std::uint8_t> adjustment(&minTree, &maxTree, adjacency);
    BuilderAccess<std::uint8_t>::Context context;
    const NodeId leaf = maxTree.topology().smallestNode(0);
    const NodeId parent = maxTree.topology().parent(leaf);
    BuilderAccess<std::uint8_t>::prepareLeaf(context, leaf, maxTree, minTree);
    requireVectorEqual(context.supportPixels, std::vector<PixelId>{0}, "initial leaf support");
    const NodeId firstBoundaryNode = minTree.topology().smallestNode(1);
    adjustment.pruneMaxLeafAndUpdateMinTree(leaf, context.supportPixels, context.dualProperPartNode, {},
                                            context.wholeSupportSmallestNode);
    requireVectorEqual(adjustment.getAdjacentNodes(), std::vector<NodeId>{firstBoundaryNode}, "updater discovers initial leaf boundary");
    require(maxTree.topology().isLeaf(parent), "primal contraction creates the next extremal leaf");
    BuilderAccess<std::uint8_t>::prepareLeaf(context, parent, maxTree, minTree);
    std::ranges::sort(context.supportPixels);
    requireVectorEqual(context.supportPixels, std::vector<PixelId>{0, 1, 2}, "merged leaf support from proper parts");
    const NodeId secondBoundaryNode = minTree.topology().smallestNode(3);
    adjustment.pruneMaxLeafAndUpdateMinTree(parent, context.supportPixels, context.dualProperPartNode, {}, context.wholeSupportSmallestNode);
    requireVectorEqual(adjustment.getAdjacentNodes(), std::vector<NodeId>{secondBoundaryNode}, "updater discovers merged leaf external boundary");
    const auto leveledImage = makeImage(1, 4, {0, 0, 0, 0});
    requireExactReconstruction(minTree, leveledImage, "min-tree after consecutive delegated boundary scans");
    requireExactReconstruction(maxTree, leveledImage, "max-tree after consecutive delegated boundary scans");

    const auto constantImage = makeImage(1, 4, {7, 7, 7, 7});
    auto constantMin = MorphologicalTreeFactory::createMinTree(constantImage, adjacency);
    auto constantMax = MorphologicalTreeFactory::createMaxTree(constantImage, adjacency);
    requireThrows<std::runtime_error>(
        [&] { BuilderAccess<std::uint8_t>::prepareLeaf(context, constantMax.topology().root(), constantMax, constantMin); },
        "the terminal root is rejected as a residual leaf after context reuse");
}

void testDelegatedBoundaryExcludesNonAdjacentDualBranches() {
    const RegularGridAdjacency2D adjacency(1, 4, 1.0);
    for (bool invertContrast : {false, true}) {
        const auto image = invertContrast ? makeImage(1, 4, {0, 2, 1, 2}) : makeImage(1, 4, {3, 1, 2, 1});
        auto minTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
        auto maxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
        auto& primal = invertContrast ? minTree : maxTree;
        auto& dual = invertContrast ? maxTree : minTree;
        adjust::DualMinMaxTreeIncrementalFilterLeaf<std::uint8_t> adjustment(&minTree, &maxTree, adjacency);
        BuilderAccess<std::uint8_t>::Context context;
        const NodeId leaf = primal.topology().smallestNode(0);
        const NodeId adjacentNode = dual.topology().smallestNode(1);
        const NodeId nonAdjacentNode = dual.topology().smallestNode(3);
        require(adjacentNode != nonAdjacentNode, "equal-level dual branches are distinct before leveling");
        BuilderAccess<std::uint8_t>::prepareLeaf(context, leaf, primal, dual);
        requireEqual(context.wholeSupportSmallestNode, context.dualProperPartNode, "whole-part certificate remains available");
        if (invertContrast)
            adjustment.pruneMinLeafAndUpdateMaxTree(leaf, context.supportPixels, context.dualProperPartNode, {}, context.wholeSupportSmallestNode);
        else
            adjustment.pruneMaxLeafAndUpdateMinTree(leaf, context.supportPixels, context.dualProperPartNode, {}, context.wholeSupportSmallestNode);
        requireVectorEqual(adjustment.getAdjacentNodes(), std::vector<NodeId>{adjacentNode}, "only the spatially adjacent dual branch is discovered");
        require(dual.topology().smallestNode(0) != dual.topology().smallestNode(3), "non-adjacent equal-level regions remain separate");
        const auto leveledImage = invertContrast ? makeImage(1, 4, {2, 2, 1, 2}) : makeImage(1, 4, {1, 1, 2, 1});
        requireExactReconstruction(minTree, leveledImage, "min-tree after delegated boundary discovery");
        requireExactReconstruction(maxTree, leveledImage, "max-tree after delegated boundary discovery");
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

void testOptionsAndInputContracts() {
    const auto image = makeImage(2, 3, {0, 2, 1, 2, 0, 0});
    const RegularGridAdjacency2D adjacency(2, 3, 1.0);

    sdrt::ExtremumResidualTreeOptions extremumOptions;
    extremumOptions.spatialOrder = sdrt::SpatialOrder(std::vector<PixelId>{5, 4, 3, 2, 1, 0});
    const auto extremumWithOptions = MorphologicalTreeFactory::createExtremumResidualTreeByCTs(image, adjacency, extremumOptions);
    auto optionMinTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto optionMaxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
    sdrt::ExtremumResidualTreeByCTsBuilder<std::uint8_t> extremumBuilder(adjacency, extremumOptions);
    extremumBuilder.build(image, std::move(optionMinTree), std::move(optionMaxTree));
    require(!extremumBuilder.spatialOrder().isRowMajor(), "extremum builder retains its explicit spatial order");
    requireVectorEqual(parentBuffer(extremumWithOptions), toVector(extremumBuilder.parents()), "extremum option parents");
    requireVectorEqual(smallestNodeMapBuffer(extremumWithOptions), toVector(extremumBuilder.smallestNodeMap()),
                       "extremum option smallest-node map");
    requireVectorEqual(altitudeBuffer(extremumWithOptions), toVector(extremumBuilder.nodeAltitudes()), "extremum option altitudes");

    auto nullMinTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto nullMaxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
    requireThrows<std::invalid_argument>([&] { extremumBuilder.build(ImageUInt8Ptr{}, std::move(nullMinTree), std::move(nullMaxTree)); },
                                         "direct residual builder rejects a null image before consuming its component trees");
    requireThrows<std::logic_error>([&] { static_cast<void>(extremumBuilder.parents()); }, "a failed build clears the previous result");

    sdrt::ExtremumResidualTreeByCTsBuilder<std::uint8_t> reusable(adjacency);
    auto firstMinTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto secondMinTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    requireThrows<std::invalid_argument>([&] { reusable.build(image, std::move(firstMinTree), std::move(secondMinTree)); },
                                         "extremum builder rejects a min-tree in the max-tree slot");
    requireThrows<std::logic_error>([&] { static_cast<void>(reusable.parents()); }, "failed construction leaves no observable partial result");
    auto validMinTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto validMaxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
    reusable.build(image, std::move(validMinTree), std::move(validMaxTree));
    requireEqual(reusable.smallestNodeMap().size(), static_cast<std::size_t>(image->getSize()), "builder reuse after failure");

    const std::array<GridOffset2D, 1> originOnly{{{0, 0}}};
    const auto disconnected = RegularGridAdjacency2D::fromStructuringElement(2, 3, originOnly);
    sdrt::ExtremumResidualTreeByCTsBuilder<std::uint8_t> disconnectedBuilder(disconnected);
    auto connectedMinTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto connectedMaxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
    requireThrows<std::invalid_argument>([&] { disconnectedBuilder.build(image, std::move(connectedMinTree), std::move(connectedMaxTree)); },
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

} // namespace

int main() {
    testFactoryAndSaliencyIntegration();
    testContrastInversion();
    testSelfDualResidualSchedule();
    testResidualTreeAssemblerInIsolation();
    testResidualAssociationsAgainstFullSupports();
    testResidualLeafPreparation();
    testDelegatedBoundaryExcludesNonAdjacentDualBranches();
    testAltitudeTypesAndDegenerateDomain();
    testOptionsAndInputContracts();
    testSharedArbitrarySymmetricAdjacency();
    return 0;
}
