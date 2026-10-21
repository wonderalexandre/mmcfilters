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
/** @brief Exposes private construction stages to the stepwise pending support oracle. @tparam T Image altitude type. */
template <AltitudeValue T> struct ExtremumResidualTreeTestAccess {
    using Builder = ExtremumResidualTreeByCTsBuilder<T>;                           ///< Production builder checked by the oracle.
    using Polarity = typename Builder::Polarity;                                   ///< Private component-tree side exposed only to the test oracle.
    using ExtremumCandidate = typename Builder::ExtremumCandidate;                 ///< Private candidate fields exposed only to ordering tests.
    using ExtremumCandidatePriority = typename Builder::ExtremumCandidatePriority; ///< Production queue comparator exposed only to priority tests.
    using Context = typename Builder::SelectedSupport;                             ///< Selected support and dual proper-part certificates.

    /**
     * @brief Collects a selected leaf support through the production helper.
     * @param context Reusable selected-support buffer.
     * @param leaf Selected primal leaf identifier.
     * @param primal Component tree containing that leaf.
     * @param dual Opposite component tree used to establish proper-part certificates.
     */
    static void prepareLeaf(Context& context, NodeId leaf, const ValuedMorphologicalTree<T>& primal, const ValuedMorphologicalTree<T>& dual) {
        context.collectFromLeaf(leaf, primal, dual);
    }
    /** @brief Test facade for residual construction and pending support inspection. */
    struct Construction {
        typename Builder::ResidualTreeConstruction extTree; ///< Production residual tree construction inspected after each test step.
        /**
         * @brief Allocates an isolated residual tree construction.
         * @param numPixels Number of domain pixels.
         * @param maxTreeNodeSlots Initial max-tree node slots.
         * @param minTreeNodeSlots Initial min-tree node slots.
         */
        Construction(std::size_t numPixels, std::size_t maxTreeNodeSlots, std::size_t minTreeNodeSlots)
            : extTree(numPixels, maxTreeNodeSlots, minTreeNodeSlots) {}
        /**
         * @brief Adds a test support through the production construction routine.
         * @param selectedSupport Pixels in the selected flat zone.
         * @param selectedPolarity Component-tree side of the selected leaf.
         * @param supportAltitude Valuation of the support before leveling.
         * @param primalLeafNode Selected primal leaf identifier.
         * @return Newly recorded residual node identifier.
         */
        NodeId addSelectedSupport(const Context& selectedSupport, Polarity selectedPolarity, T supportAltitude, NodeId primalLeafNode) {
            return extTree.addSelectedSupport(selectedSupport, ExtremumCandidate{selectedSupport.pixels.size(), InvalidPixel, selectedPolarity, primalLeafNode},
                                              supportAltitude);
        }
        /** @brief Counts recorded nonterminal supports. @return Number of residual nodes excluding the reserved root. */
        std::size_t numRecordedNodes() const { return extTree.nodeParents.size() - 1; }
        /** @brief Enumerates the implicit pending set. @return Recorded non-root nodes whose parent remains undefined. */
        std::vector<NodeId> pendingResidualNodes() const {
            std::vector<NodeId> pendingNodes;
            for (NodeId residualNode = 1; residualNode < static_cast<NodeId>(extTree.nodeParents.size()); ++residualNode)
                if (extTree.nodeParents[static_cast<std::size_t>(residualNode)] == InvalidNode)
                    pendingNodes.push_back(residualNode);
            return pendingNodes;
        }
        /**
         * @brief Reads the pending supports contained in one component proper part.
         * @param polarity Component-tree side.
         * @param properPartNode Component node owning the proper part.
         * @return Borrowed list of pending supports contained in that proper part.
         */
        const auto& pendingSupportsIn(Polarity polarity, NodeId properPartNode) const {
            return (polarity == Polarity::Maximum ? extTree.maxTreePendingSupports : extTree.minTreePendingSupports).pendingSupportsIn(properPartNode);
        }
        /**
         * @brief Updates pending supports through the production routine.
         * @tparam ProperPartTransfer Whole proper-part transfer record type.
         * @param residualNode Newly recorded residual support.
         * @param selectedPolarity Component-tree side of the selected leaf.
         * @param primalDestinationNode Live primal destination of the selected support.
         * @param dualDestinationNode Live dual destination of the selected support.
         * @param dualProperPartTransfers Whole dual proper-part transfers reported by the updater.
         */
        template <class ProperPartTransfer>
        void updatePendingSupports(NodeId residualNode, Polarity selectedPolarity, NodeId primalDestinationNode, NodeId dualDestinationNode,
                                   std::span<const ProperPartTransfer> dualProperPartTransfers) {
            extTree.updatePendingSupports(residualNode, selectedPolarity, primalDestinationNode, dualDestinationNode, dualProperPartTransfers);
        }
        /**
         * @brief Completes the test hierarchy through the production construction routine.
         * @param image Original image defining the pixel domain.
         * @param adjacency Shared construction adjacency attached to the descriptor.
         * @param rootAltitude Valuation of the terminal root.
         * @return Completed residual topology and valuation buffers.
         */
        auto completeTree(const ImagePtr<T>& image, const RegularGridAdjacency2D& adjacency, T rootAltitude) && {
            return std::move(extTree).completeTree(
                image, makeMorphologicalTreeSemantics(MorphologicalTreeKind::ExtremumResidualTree, SharedAdjacencyContext{adjacency}), rootAltitude);
        }
    };
};
} // namespace mmcfilters::detail

template <class T> using BuilderAccess = mmcfilters::detail::ExtremumResidualTreeTestAccess<T>;
template <class T> using TestResidualTreeConstruction = typename BuilderAccess<T>::Construction;

namespace {

template <class Engine>
concept HasInfinityPixelAccessor = requires(const Engine& engine) { engine.infinityPixel(); };

template <class Construction>
concept LvalueCompletableConstruction =
    requires(Construction& extTree, const ImagePtr<int>& image, const RegularGridAdjacency2D& adjacency) { extTree.completeTree(image, adjacency, 0); };

template <class Construction>
concept RvalueCompletableConstruction = requires(Construction&& extTree, const ImagePtr<int>& image, const RegularGridAdjacency2D& adjacency) {
    std::move(extTree).completeTree(image, adjacency, 0);
};

static_assert(!HasInfinityPixelAccessor<ExtremumResidualTreeByCTsBuilder<std::uint8_t>>, "extremum construction must not expose infinity-pixel state");
static_assert(!LvalueCompletableConstruction<TestResidualTreeConstruction<int>>, "residual tree construction must not be completable through an lvalue");
static_assert(RvalueCompletableConstruction<TestResidualTreeConstruction<int>>, "residual tree construction must be completable through an rvalue");
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

void testExtremumCandidatePriority() {
    using Access = BuilderAccess<int>;
    using Polarity = Access::Polarity;
    const Access::ExtremumCandidatePriority lowerPriority;
    const Access::ExtremumCandidate first{2, 0, Polarity::Maximum, 1};
    const Access::ExtremumCandidate second{2, 3, Polarity::Minimum, 2};
    const Access::ExtremumCandidate sameKey{2, 0, Polarity::Minimum, 3};
    require(lowerPriority(second, first), "smaller spatial minima have higher residual priority");
    require(!lowerPriority(first, first), "equal residual keys are equivalent");
    require(!lowerPriority(first, sameKey) && !lowerPriority(sameKey, first), "polarity and node identity do not affect selection priority");

    require(lowerPriority({2, 0, Polarity::Minimum, 2}, {1, 3, Polarity::Maximum, 1}), "support cardinality takes precedence over spatial order");
    require(!lowerPriority(first, second), "residual priority is asymmetric");
}

void testResidualTreeConstructionInIsolation() {
    using Polarity = BuilderAccess<int>::Polarity;
    const auto image = makeTypedImage<int>(1, 4, {5, 5, 3, 1});
    const RegularGridAdjacency2D adjacency(1, 4, 1.0);
    auto minTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto maxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
    detail::DualMinMaxTreeIncrementalFilterLeaf<int> adjustment(&minTree, &maxTree, adjacency);
    TestResidualTreeConstruction<int> extTree(4, maxTree.topology().numInternalNodeSlots(), minTree.topology().numInternalNodeSlots());
    BuilderAccess<int>::Context context;

    NodeId leaf = maxTree.topology().smallestNode(0);
    NodeId parent = maxTree.topology().parent(leaf);
    BuilderAccess<int>::prepareLeaf(context, leaf, maxTree, minTree);
    const NodeId firstResidualNode = extTree.addSelectedSupport(context, Polarity::Maximum, 5, leaf);
    requireEqual(firstResidualNode, NodeId{1}, "residual construction first recorded support id");
    requireVectorEqual(extTree.pendingResidualNodes(), std::vector<NodeId>{1}, "Xi after recording");
    adjustment.pruneMaxLeafAndUpdateMinTree(leaf, context.pixels, context.dualProperPartNode, context.wholeDualProperPartNode);
    extTree.updatePendingSupports(firstResidualNode, Polarity::Maximum, parent, adjustment.getLastSelectedSupportNode(),
                                  adjustment.getLastProperPartTransfers());
    for (Polarity polarity : {Polarity::Maximum, Polarity::Minimum}) {
        const auto& tree = polarity == Polarity::Maximum ? maxTree : minTree;
        const auto& associated = extTree.pendingSupportsIn(polarity, tree.topology().smallestNode(0));
        requireVectorEqual(std::vector<NodeId>(associated.begin(), associated.end()), std::vector<NodeId>{1}, "explicit association after leveling");
    }
    leaf = maxTree.topology().smallestNode(0);
    parent = maxTree.topology().parent(leaf);
    BuilderAccess<int>::prepareLeaf(context, leaf, maxTree, minTree);
    const NodeId secondResidualNode = extTree.addSelectedSupport(context, Polarity::Maximum, 3, leaf);
    requireEqual(secondResidualNode, NodeId{2}, "residual construction second recorded support id");
    requireVectorEqual(extTree.pendingResidualNodes(), std::vector<NodeId>{2}, "Xi removes attached child");
    adjustment.pruneMaxLeafAndUpdateMinTree(leaf, context.pixels, context.dualProperPartNode, context.wholeDualProperPartNode);
    extTree.updatePendingSupports(secondResidualNode, Polarity::Maximum, parent, adjustment.getLastSelectedSupportNode(),
                                  adjustment.getLastProperPartTransfers());
    requireEqual(extTree.numRecordedNodes(), std::size_t{2}, "residual construction recorded support count");

    auto output = std::move(extTree).completeTree(image, adjacency, 1);
    requireVectorEqual(output.parent, std::vector<NodeId>{0, 2, 0}, "residual construction parent materialization");
    requireVectorEqual(output.smallestNodeMap, std::vector<NodeId>{1, 1, 2, 0}, "residual construction smallest-node-map materialization");
    requireVectorEqual(output.nodeAltitudes, std::vector<int>{1, 5, 3}, "residual construction altitude materialization");
}

// Independently check Xi and Eq. (31) against complete recorded pixel sets
// after every leveling, including pending supports outside the selected leaf.
template <AltitudeValue T> void checkResidualAssociations(const ImagePtr<T>& image, const RegularGridAdjacency2D& adjacency) {
    using Polarity = typename BuilderAccess<T>::Polarity;
    auto minTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto maxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
    detail::DualMinMaxTreeIncrementalFilterLeaf<T> adjustment(&minTree, &maxTree, adjacency);
    auto leveledImage = Image<T>::create(image->getNumRows(), image->getNumColumns());
    for (PixelId pixel = 0; pixel < image->getSize(); ++pixel)
        (*leveledImage)[pixel] = (*image)[pixel];
    const auto initialMaxSlots = maxTree.topology().numInternalNodeSlots();
    const auto initialMinSlots = minTree.topology().numInternalNodeSlots();
    TestResidualTreeConstruction<T> extTree(image->getSize(), maxTree.topology().numInternalNodeSlots(), minTree.topology().numInternalNodeSlots());
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
        auto actualXi = extTree.pendingResidualNodes();
        std::sort(actualXi.begin(), actualXi.end());
        requireVectorEqual(actualXi, expectedXi, "Xi equals maximal recorded supports");
        for (Polarity polarity : {Polarity::Maximum, Polarity::Minimum}) {
            const auto& tree = polarity == Polarity::Maximum ? maxTree : minTree;
            const auto& topology = tree.topology();
            requireEqual(topology.numInternalNodeSlots(), polarity == Polarity::Maximum ? initialMaxSlots : initialMinSlots,
                         "leaf updates preserve the allocated component-tree node slots");
            for (NodeId node : topology.aliveNodeIds()) {
                require(topology.properPartCardinality(node) > 0, "live component-tree nodes retain nonempty proper parts");
                PixelId expectedMinimum = std::numeric_limits<PixelId>::max();
                for (PixelId pixel : topology.nodeSupport(node))
                    expectedMinimum = std::min(expectedMinimum, pixel);
                requireEqual(adjustment.minPixel(tree, node), expectedMinimum, "component pixel minima agree with complete current supports");
            }
            std::vector<std::vector<NodeId>> expected(topology.numInternalNodeSlots());
            for (NodeId node : expectedXi) {
                const NodeId owner = topology.smallestNode(supports[node].front());
                require(topology.isAlive(owner), "pending support has a live component owner");
                for (PixelId pixel : supports[node])
                    requireEqual(topology.smallestNode(pixel), owner, "complete support satisfies the association relation");
                expected[owner].push_back(node);
            }
            for (NodeId owner = 0; owner < topology.numInternalNodeSlots(); ++owner) {
                const auto& associated = extTree.pendingSupportsIn(polarity, owner);
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
        for (PixelId pixel : context.pixels)
            requireEqual(dual.topology().smallestNode(pixel), context.dualProperPartNode, "complete flat zone lies in one dual proper part");
        require(best.first >= previousArea, "selected support cardinalities are non-decreasing");
        previousArea = best.first;
        const T nodeAltitude = primal.nodeAltitude(selected);
        const NodeId primalParent = primal.topology().parent(selected);
        const T parentAltitude = primal.nodeAltitude(primalParent);
        const auto signedResidualValue = static_cast<AltitudeDifference<T>>(nodeAltitude) -
                                         static_cast<AltitudeDifference<T>>(parentAltitude);
        require(polarity == Polarity::Maximum ? signedResidualValue > 0 : signedResidualValue < 0,
                "residual sign agrees with the selected polarity");
        for (NodeId node = 1; node < static_cast<NodeId>(pending.size()); ++node) {
            if (pending[node] && std::all_of(supports[node].begin(), supports[node].end(), [&](PixelId pixel) {
                    return std::find(context.pixels.begin(), context.pixels.end(), pixel) != context.pixels.end();
                })) {
                pending[node] = false;
                expectedParents[node] = static_cast<NodeId>(supports.size());
            }
        }
        supports.emplace_back(context.pixels);
        pending.push_back(true);
        expectedParents.push_back(0);
        const NodeId residualNode = extTree.addSelectedSupport(context, polarity, nodeAltitude, selected);
        if (polarity == Polarity::Maximum)
            adjustment.pruneMaxLeafAndUpdateMinTree(selected, context.pixels, context.dualProperPartNode, context.wholeDualProperPartNode);
        else
            adjustment.pruneMinLeafAndUpdateMaxTree(selected, context.pixels, context.dualProperPartNode, context.wholeDualProperPartNode);
        for (PixelId pixel : context.pixels) {
            (*leveledImage)[pixel] = parentAltitude;
            requireEqual(primal.topology().smallestNode(pixel), primalParent, "primal pruning establishes containment over the complete support");
            requireEqual(dual.topology().smallestNode(pixel), adjustment.getLastSelectedSupportNode(),
                         "dual transfer destination contains the complete selected support");
        }
        extTree.updatePendingSupports(residualNode, polarity, primalParent, adjustment.getLastSelectedSupportNode(), adjustment.getLastProperPartTransfers());
        requireRebuiltComponentTrees(minTree, maxTree, leveledImage, adjacency, "stepwise dual leaf update");
        verify();
    }
    const auto output = std::move(extTree).completeTree(image, adjacency, minTree.nodeAltitude(minTree.topology().root()));
    requireVectorEqual(output.parent, expectedParents, "tree completion attaches implicit Xi to the root and preserves assigned parents");
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
    checkResidualAssociations(makeTypedImage<std::int8_t>(2, 3, {-4, 2, -1, 2, -4, 0}), RegularGridAdjacency2D(2, 3, 1.0));
    checkResidualAssociations(makeTypedImage<double>(2, 3, {-1.5, 0.25, -0.5, 0.25, -1.5, 0.0}), RegularGridAdjacency2D(2, 3, 1.0));
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
    requireVectorEqual(context.pixels, std::vector<PixelId>{0}, "initial leaf support");
    adjustment.pruneMaxLeafAndUpdateMinTree(leaf, context.pixels, context.dualProperPartNode, context.wholeDualProperPartNode);
    requireRebuiltComponentTrees(minTree, maxTree, makeImage(1, 4, {1, 1, 1, 0}), adjacency,
                                "updater discovers the initial leaf boundary");
    require(maxTree.topology().isLeaf(parent), "primal contraction creates the next extremal leaf");
    BuilderAccess<std::uint8_t>::prepareLeaf(context, parent, maxTree, minTree);
    std::ranges::sort(context.pixels);
    requireVectorEqual(context.pixels, std::vector<PixelId>{0, 1, 2}, "merged leaf support from proper parts");
    adjustment.pruneMaxLeafAndUpdateMinTree(parent, context.pixels, context.dualProperPartNode, context.wholeDualProperPartNode);
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
        for (PixelId pixel : context.pixels)
            (*expectedImage)[pixel] = targetAltitude;
        detail::DualMinMaxTreeIncrementalFilterLeaf<T> adjustment(&minTree, &maxTree, adjacency);
        if (pruneMaximum)
            adjustment.pruneMaxLeafAndUpdateMinTree(leaf, context.pixels, context.dualProperPartNode, context.wholeDualProperPartNode);
        else
            adjustment.pruneMinLeafAndUpdateMaxTree(leaf, context.pixels, context.dualProperPartNode, context.wholeDualProperPartNode);
        reference.pruneNode(leaf);
        requireVectorEqual(parentBuffer(primal), parentBuffer(reference), "internal leaf contraction preserves public prune topology");
        requireEqual(primal.topology().getMutationVersion(), reference.topology().getMutationVersion(),
                     "internal leaf contraction preserves mutation versions");
        requireRebuiltComponentTrees(minTree, maxTree, expectedImage, adjacency, "certified leaf update");
    }
}

void testLeafBoundaryExcludesNonAdjacentDualBranches() {
    const RegularGridAdjacency2D adjacency(1, 4, 1.0);
    for (bool useWholePartCertificate : {false, true}) {
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
            requireEqual(context.wholeDualProperPartNode, context.dualProperPartNode, "whole-part certificate remains available");
            const NodeId certificate = useWholePartCertificate ? context.wholeDualProperPartNode : InvalidNode;
            if (invertContrast)
                adjustment.pruneMinLeafAndUpdateMaxTree(leaf, context.pixels, context.dualProperPartNode, certificate);
            else
                adjustment.pruneMaxLeafAndUpdateMinTree(leaf, context.pixels, context.dualProperPartNode, certificate);
            require(dual.topology().smallestNode(0) != dual.topology().smallestNode(3), "non-adjacent equal-level regions remain separate");
            const auto leveledImage = invertContrast ? makeImage(1, 4, {2, 2, 1, 2}) : makeImage(1, 4, {1, 1, 2, 1});
            requireRebuiltComponentTrees(minTree, maxTree, leveledImage, adjacency, "spatial boundary excludes the non-adjacent branch");
        }
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
        requireVectorEqual(context.pixels, std::vector<PixelId>{1}, "selected leaf support is one pixel");
        requireEqual(context.dualProperPartNode, dual.topology().smallestNode(3), "dual proper part also contains an unselected extremum");
        requireEqual(context.wholeDualProperPartNode, InvalidNode, "selected support is a strict subset of the dual proper part");

        detail::DualMinMaxTreeIncrementalFilterLeaf<std::uint8_t> adjustment(&minTree, &maxTree, adjacency);
        if (invertContrast) {
            adjustment.pruneMinLeafAndUpdateMaxTree(leaf, context.pixels, context.dualProperPartNode, context.wholeDualProperPartNode);
        } else {
            adjustment.pruneMaxLeafAndUpdateMinTree(leaf, context.pixels, context.dualProperPartNode, context.wholeDualProperPartNode);
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

void testBuilderAndFactoryEquivalence() {
    const auto image = makeImage(2, 3, {0, 2, 1, 2, 0, 0});
    const RegularGridAdjacency2D adjacency(2, 3, 1.0);

    const auto factoryTree = MorphologicalTreeFactory::createExtremumResidualTreeByCTs(image, adjacency);
    ExtremumResidualTreeByCTsBuilder<std::uint8_t> extremumBuilder(adjacency);
    const auto hierarchy = extremumBuilder.build(image);
    requireVectorEqual(parentBuffer(factoryTree), hierarchy.parent, "factory and builder parents");
    requireVectorEqual(smallestNodeMapBuffer(factoryTree), hierarchy.smallestNodeMap,
                       "factory and builder smallest-node map");
    requireVectorEqual(altitudeBuffer(factoryTree), hierarchy.nodeAltitudes, "factory and builder altitudes");

    const auto rebuilt = extremumBuilder.build(image);
    requireVectorEqual(hierarchy.parent, rebuilt.parent, "builder reuse preserves parents");
    requireVectorEqual(hierarchy.smallestNodeMap, rebuilt.smallestNodeMap, "builder reuse preserves smallest-node map");
    requireVectorEqual(hierarchy.nodeAltitudes, rebuilt.nodeAltitudes, "builder reuse preserves altitudes");
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
    testExtremumCandidatePriority();
    testResidualTreeConstructionInIsolation();
    testResidualAssociationsAgainstFullSupports();
    testResidualLeafPreparation();
    testLeafBoundaryExcludesNonAdjacentDualBranches();
    testPartialDualProperPartPreservesReattachments();
    testCertifiedLeafUpdate<std::uint8_t>();
    testCertifiedLeafUpdate<std::int32_t>();
    testCertifiedLeafUpdate<double>();
    testAltitudeTypesAndDegenerateDomain();
    testBuilderAndFactoryEquivalence();
    testSharedArbitrarySymmetricAdjacency();
    testOwnedBuildResultLifetime();
    return 0;
}
