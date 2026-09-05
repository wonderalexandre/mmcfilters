#include "support/TestSupport.hpp"

#include "mmcfilters/trees/TreeOfShapesProducer.hpp"
#include "mmcfilters/trees/detail/ComponentTreeUnionFind.hpp"
#include "mmcfilters/trees/detail/CommittedTreeAccess.hpp"
#include "mmcfilters/trees/TreeEditor.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <random>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

static_assert(sizeof(ToSFloodDepth) >= sizeof(std::uint32_t), "Tree-of-shapes flood depth must not wrap at 16 bits.");

// Independent parent-path oracle: no DFS intervals or RMQ data are consulted.
static void requireSharedLcaMatchesParentPaths(const MorphologicalTree& tree) {
    for (NodeId first : tree.aliveNodeIds()) {
        const auto firstAncestors = collectNodeIds(tree.ancestors(first));
        for (NodeId second : tree.aliveNodeIds()) {
            NodeId expected = InvalidNode;
            for (NodeId candidate : tree.ancestors(second)) {
                if (std::find(firstAncestors.begin(), firstAncestors.end(), candidate) != firstAncestors.end()) {
                    expected = candidate;
                    break;
                }
            }
            requireEqual(tree.lowestCommonAncestor(first, second), expected, "public LCA matches parent-path oracle");
            requireEqual(detail::CommittedTreeAccess::lowestCommonAncestor(tree, first, second), expected,
                         "internal LCA matches parent-path oracle");
        }
    }
    requireEqual(tree.lowestCommonAncestor(InvalidNode, tree.root()), InvalidNode, "invalid LCA endpoint");
    requireEqual(tree.lowestCommonAncestor(tree.root(), tree.numInternalNodeSlots()), InvalidNode, "out-of-range LCA endpoint");
}

int main() {
    auto requireThrows = [](auto&& fn, const std::string& label) {
        bool threw = false;
        try {
            fn();
        } catch (const std::exception&) {
            threw = true;
        }
        require(threw, label);
    };

    {
        auto image = makeComponentTreeFixture();
        RegularGridAdjacency2D adj(4, 4, 1.5);

        detail::ComponentTreeUnionFind maxBuilder(&adj, true);
        auto [parentMax, orderMax, numNodesMax] = maxBuilder.build(image);
        requireEqual(numNodesMax, 6, "ComponentTreeUnionFind max numNodes");
        requireEqual(static_cast<int>(orderMax.size()), 16, "ComponentTreeUnionFind max order size");
        requireVectorEqual(orderMax, std::vector<int>{15, 8, 12, 13, 2, 3, 7, 11, 0, 1, 4, 5, 6, 9, 10, 14}, "ComponentTreeUnionFind max ordered pixels");
        requireVectorEqual(parentMax, std::vector<int>{2, 0, 8, 2, 0, 0, 5, 2, 15, 5, 5, 2, 8, 8, 10, 15}, "ComponentTreeUnionFind max parent array");

        detail::ComponentTreeUnionFind minBuilder(&adj, false);
        auto [parentMin, orderMin, numNodesMin] = minBuilder.build(image);
        requireEqual(numNodesMin, 6, "ComponentTreeUnionFind min numNodes");
        requireVectorEqual(orderMin, std::vector<int>{10, 14, 5, 6, 9, 0, 1, 4, 2, 3, 7, 11, 8, 12, 13, 15}, "ComponentTreeUnionFind min ordered pixels");
        requireVectorEqual(parentMin, std::vector<int>{5, 0, 0, 2, 0, 10, 5, 2, 0, 5, 10, 2, 8, 8, 10, 2}, "ComponentTreeUnionFind min parent array");
    }

    {
        auto image = makeImage(2, 2, {0, 0, 0, 0});
        RegularGridAdjacency2D adj(2, 2, 1.5);

        detail::ComponentTreeUnionFind maxBuilder(&adj, true);
        auto [parentMax, orderMax, numNodesMax] = maxBuilder.build(image);
        requireEqual(numNodesMax, 1, "ComponentTreeUnionFind constant max numNodes");
        requireEqual(static_cast<int>(orderMax.size()), 4, "ComponentTreeUnionFind constant max order size");

        detail::ComponentTreeUnionFind minBuilder(&adj, false);
        auto [parentMin, orderMin, numNodesMin] = minBuilder.build(image);
        requireEqual(numNodesMin, 1, "ComponentTreeUnionFind constant min numNodes");
        requireEqual(static_cast<int>(orderMin.size()), 4, "ComponentTreeUnionFind constant min order size");

        auto maxTree = MorphologicalTreeFactory::createMaxTree(image);
        auto minTree = MorphologicalTreeFactory::createMinTree(image);
        requireEqual(maxTree.topology().numNodes(), 1, "constant max-tree node count");
        requireEqual(minTree.topology().numNodes(), 1, "constant min-tree node count");
        requireVectorEqual(collectNodeIds(maxTree.topology().ancestors(maxTree.topology().root())), std::vector<NodeId>{maxTree.topology().root()},
                           "singleton ancestors include self");
        requireVectorEqual(collectNodeIds(maxTree.topology().descendants(maxTree.topology().root())), std::vector<NodeId>{},
                           "singleton has no strict descendants");
        requireVectorEqual(collectNodeIds(maxTree.topology().subtreeNodes(maxTree.topology().root())), std::vector<NodeId>{maxTree.topology().root()},
                           "singleton subtree includes root");
        requireEqual(maxTree.topology().dfsEntryIndex(maxTree.topology().root()), 0, "singleton DFS entry index");
        requireEqual(maxTree.topology().dfsExitIndex(maxTree.topology().root()), 1, "singleton DFS exit index");
        requireVectorEqual(collectPixelIds(maxTree.topology().properPart(maxTree.topology().root())), std::vector<PixelId>{0, 1, 2, 3},
                           "constant max-tree proper parts");
        requireVectorEqual(collectPixelIds(minTree.topology().properPart(minTree.topology().root())), std::vector<PixelId>{0, 1, 2, 3},
                           "constant min-tree proper parts");

        auto valuedTree = MorphologicalTreeFactory::createMaxTree(image);
        requireVectorEqual(collectImageValues(valuedTree.reconstructFromNodeAltitudes()), collectImageValues(image), "constant valuedTree reconstruction");
    }

    {
        auto image = ImageFloat::create(2, 3);
        float values[] = {2.0f, 1.0f, 2.0f, 0.5f, 3.0f, 1.5f};
        for (int i = 0; i < 6; ++i) {
            (*image)[i] = values[i];
        }
        RegularGridAdjacency2D adj(2, 3, 1.5);
        detail::ComponentTreeUnionFind maxBuilder(&adj, true);
        auto ordered = maxBuilder.sort(image);
        requireEqual(static_cast<int>(ordered.size()), 6, "ComponentTreeUnionFind float sort size");
        std::vector<float> sortedValues;
        for (PixelId pixelId : ordered) {
            sortedValues.push_back((*image)[pixelId]);
        }
        require(std::is_sorted(sortedValues.begin(), sortedValues.end()), "ComponentTreeUnionFind float sort must be nondecreasing for max-tree");
    }

    {
        auto image = ImageFloat::create(1, 2);
        (*image)[0] = -0.0f;
        (*image)[1] = +0.0f;
        auto tree = MorphologicalTreeFactory::createMaxTree(image, 1.5);
        requireEqual(tree.topology().numNodes(), 1, "signed-zero plateau node count");
        require(std::signbit(tree.nodeAltitude(tree.topology().root())),
                "component-tree altitude must retain the first row-major proper-part representation");
    }

    {
        detail::AdjacencyUC adjacency(3, 3, true);
        requireVectorEqual(collectNodeIds(adjacency.getNeighborIndices(1, 1)), std::vector<NodeId>{1, 3, 7, 5}, "AdjacencyUC default 4-neighborhood");
        adjacency.setDiagonalConnection(1, 1, detail::DiagonalConnection::Sw);
        adjacency.setDiagonalConnection(1, 1, detail::DiagonalConnection::Ne);
        adjacency.setDiagonalConnection(1, 1, detail::DiagonalConnection::Se);
        adjacency.setDiagonalConnection(1, 1, detail::DiagonalConnection::Nw);
        requireVectorEqual(collectNodeIds(adjacency.getNeighborIndices(1, 1)), std::vector<NodeId>{1, 3, 7, 5, 6, 2, 8, 0},
                           "AdjacencyUC diagonal neighborhood");

        detail::PriorityQueueToS queue;
        queue.initial(10, 5);
        queue.priorityPush(20, 3, 4);
        queue.priorityPush(30, 7, 9);
        requireEqual(queue.priorityPop(), 10, "PriorityQueueToS first pop");
        requireEqual(queue.getCurrentPriority(), 5, "PriorityQueueToS first priority");
        requireEqual(queue.priorityPop(), 20, "PriorityQueueToS second pop");
        requireEqual(queue.getCurrentPriority(), 4, "PriorityQueueToS second priority");
        requireEqual(queue.priorityPop(), 30, "PriorityQueueToS third pop");
        requireEqual(queue.getCurrentPriority(), 7, "PriorityQueueToS third priority");
        require(queue.isEmpty(), "PriorityQueueToS must be empty after pops");

        detail::PriorityQueueToS tieQueue;
        tieQueue.initial(40, 5);
        tieQueue.priorityPush(50, 3, 3);
        tieQueue.priorityPush(60, 7, 7);
        requireEqual(tieQueue.priorityPop(), 40, "PriorityQueueToS tie first pop");
        requireEqual(tieQueue.priorityPop(), 50, "PriorityQueueToS must prefer the lower bucket on an equal-distance tie");
    }

    {
        auto image = makeImage(3, 3,
                               {
                                   1,
                                   2,
                                   1,
                                   2,
                                   3,
                                   2,
                                   1,
                                   2,
                                   1,
                               });

        TreeOfShapesProducer tosBuilder(makeTopographicConvention(image, TestTopographicImmersion::Min4Max8));
        auto [interpolationMin, interpolationMax, adjacency] = tosBuilder.interpolateImage4c8c(image);
        requireEqual(static_cast<int>(interpolationMin.size()), 49, "TreeOfShapesProducer interpolate min size");
        requireEqual(static_cast<int>(interpolationMax.size()), 49, "TreeOfShapesProducer interpolate max size");
        requireEqual(interpolationMin[ImageUtils::to1D(3, 3, 7)], 3, "TreeOfShapesProducer interpolated center min");
        requireEqual(interpolationMax[ImageUtils::to1D(3, 3, 7)], 3, "TreeOfShapesProducer interpolated center max");
        requireVectorEqual(collectNodeIds(adjacency.getNeighborIndices(ImageUtils::to1D(3, 3, 7))), std::vector<NodeId>{17, 23, 31, 25},
                           "TreeOfShapesProducer interpolated center adjacency");

        auto [imgU, orderInterpolated, _] = tosBuilder.sort(image);
        requireEqual(static_cast<int>(imgU.size()), 49, "TreeOfShapesProducer sorted interpolation size");
        requireEqual(static_cast<int>(orderInterpolated.size()), 49, "TreeOfShapesProducer sorted order size");
        requireEqual(orderInterpolated.front(), 0, "TreeOfShapesProducer sorted first interpolated pixel");
        auto sortedInterpolatedOrder = orderInterpolated;
        std::sort(sortedInterpolatedOrder.begin(), sortedInterpolatedOrder.end());
        for (PixelId pixelId = 0; pixelId < 49; ++pixelId) {
            requireEqual(sortedInterpolatedOrder[static_cast<std::size_t>(pixelId)], pixelId, "TreeOfShapesProducer sorted order must be a permutation");
        }

        auto result = tosBuilder.build<ToSGrayLevel>(image);
        requireEqual(static_cast<int>(result.parent.size()), 6, "TreeOfShapesProducer native node count");
        requireEqual(static_cast<int>(result.smallestNodeMap.size()), 9, "TreeOfShapesProducer native proper-part count");
        for (NodeId smallestNodeId : result.smallestNodeMap) {
            require(smallestNodeId >= 0 && smallestNodeId < static_cast<NodeId>(result.parent.size()), "TreeOfShapesProducer smallest node must be valid");
        }
    }

    {
        auto saddle = makeImage(2, 2,
                                {
                                    0,
                                    8,
                                    9,
                                    4,
                                });

        TreeOfShapesProducer min4max8(makeTopographicConvention(saddle, TestTopographicImmersion::Min4Max8));
        auto [min4max8Min, min4max8Max, min4max8Adj] = min4max8.interpolateImage4c8c(saddle);
        const NodeId center = ImageUtils::to1D(2, 2, 5);
        requireEqual(min4max8Min[center], static_cast<uint8_t>(8), "Min4cMax8c saddle center min");
        requireEqual(min4max8Max[center], static_cast<uint8_t>(9), "Min4cMax8c saddle center max");
        requireVectorEqual(collectNodeIds(min4max8Adj.getNeighborIndices(center)), std::vector<NodeId>{7, 11, 17, 13, 16, 8},
                           "Min4cMax8c saddle center adjacency");

        TreeOfShapesProducer min8max4(makeTopographicConvention(saddle, TestTopographicImmersion::Min8Max4));
        auto [min8max4Min, min8max4Max, min8max4Adj] = min8max4.interpolateImage4c8c(saddle);
        requireEqual(min8max4Min[center], static_cast<uint8_t>(0), "Min8cMax4c saddle center min");
        requireEqual(min8max4Max[center], static_cast<uint8_t>(4), "Min8cMax4c saddle center max");
        requireVectorEqual(collectNodeIds(min8max4Adj.getNeighborIndices(center)), std::vector<NodeId>{7, 11, 17, 13, 18, 6},
                           "Min8cMax4c saddle center adjacency");
    }

    {
        auto halfLevelBoundary = makeImage(1, 2, {0, 1});
        TreeOfShapesProducer selfDual(makeTopographicConvention(halfLevelBoundary));
        auto [interpolationMin, interpolationMax, adjacency] = selfDual.interpolateImage(halfLevelBoundary);
        const NodeId outerCorner = ImageUtils::to1D(0, 0, 5);
        const NodeId firstOriginal = ImageUtils::to1D(1, 1, 5);
        const NodeId secondOriginal = ImageUtils::to1D(1, 3, 5);
        requireEqual(interpolationMin[outerCorner], 1, "SelfDual outer boundary median must preserve half levels in Z/2");
        requireEqual(interpolationMax[outerCorner], 1, "SelfDual outer boundary median max must preserve half levels in Z/2");
        requireEqual(interpolationMin[firstOriginal], 0, "SelfDual first original pixel must be scaled by 2");
        requireEqual(interpolationMax[secondOriginal], 2, "SelfDual second original pixel must be scaled by 2");

        auto thinBoundary = makeImage(1, 3, {10, 20, 30});
        auto [thinMin, thinMax, thinAdjacency] = selfDual.interpolateImage(thinBoundary);
        (void)thinAdjacency;
        requireEqual(thinMin[ImageUtils::to1D(0, 0, 7)], static_cast<ToSGrayLevel>(40),
                     "SelfDual thin-image boundary median must use each boundary pixel exactly once");
        requireEqual(thinMax[ImageUtils::to1D(0, 0, 7)], static_cast<ToSGrayLevel>(40), "SelfDual thin-image boundary median max");

        auto [imgU, orderInterpolated, _] = selfDual.sort(halfLevelBoundary);
        requireEqual(orderInterpolated.front(), outerCorner, "SelfDual propagation must start from the outer boundary");
        requireEqual(imgU[outerCorner], static_cast<ToSFloodDepth>(1), "SelfDual propagated boundary level must preserve Z/2 half level");

        TreeOfShapesProducer lowerRightInfinity(makeTopographicConvention(halfLevelBoundary, TestTopographicImmersion::SelfDualSpan,
                                                                           TopographicDomainExtension::ExteriorRing, PixelId{14}));
        auto [customImgU, customOrderInterpolated, __] = lowerRightInfinity.sort(halfLevelBoundary);
        const NodeId lowerRightCorner = ImageUtils::to1D(2, 4, 5);
        requireEqual(customOrderInterpolated.front(), lowerRightCorner, "SelfDual propagation must honor a custom infinity pixel");
        requireEqual(customImgU[lowerRightCorner], static_cast<ToSFloodDepth>(1), "SelfDual custom infinity pixel must preserve Z/2 half level");
        TreeOfShapesProducer internalInfinity(makeTopographicConvention(halfLevelBoundary, TestTopographicImmersion::SelfDualSpan,
                                                                         TopographicDomainExtension::ExteriorRing, PixelId{6}));
        auto [internalImgU, internalOrderInterpolated, ___] = internalInfinity.sort(halfLevelBoundary);
        const NodeId firstOriginalSeed = ImageUtils::to1D(1, 1, 5);
        requireEqual(internalOrderInterpolated.front(), firstOriginalSeed, "SelfDual propagation must honor an internal custom infinity pixel");
        requireThrows([&]() {
            static_cast<void>(TreeOfShapesProducer(makeTopographicConvention(halfLevelBoundary, TestTopographicImmersion::SelfDualSpan,
                                                                             TopographicDomainExtension::ExteriorRing, PixelId{20}))
                                  .sort(halfLevelBoundary));
        },
                      "SelfDual infinity pixel must stay inside the interpolated domain");
    }

    {
        auto image = makeImage(1, 2, {0, 1});
        TreeOfShapesProducer producer(makeTopographicConvention(image, TestTopographicImmersion::SelfDualSpan,
                                                                 TopographicDomainExtension::None, PixelId{0}));
        const auto& convention = producer.convention();
        require(convention.domainExtension == TopographicDomainExtension::None, "TreeOfShapesProducer must retain its domain extension");

        auto [interpolationMin, interpolationMax, adjacency] = producer.interpolateImage(image);
        (void)adjacency;
        requireEqual(interpolationMin.size(), std::size_t{3}, "unpadded self-dual interpolation size");
        requireEqual(interpolationMax.size(), std::size_t{3}, "unpadded self-dual interpolation max size");
        requireEqual(interpolationMin[0], static_cast<ToSGrayLevel>(0), "unpadded first original level");
        requireEqual(interpolationMax[2], static_cast<ToSGrayLevel>(2), "unpadded second original level");

        auto [treeLevel, order, _] = producer.sort(image);
        requireEqual(treeLevel.size(), std::size_t{3}, "unpadded propagation domain size");
        requireEqual(order.front(), 0, "unpadded infinity pixel is expressed in the selected immersion domain");

        requireThrows(
            [&]() {
                TreeOfShapesProducer invalidSeed(makeTopographicConvention(image, TestTopographicImmersion::SelfDualSpan,
                                                                            TopographicDomainExtension::None, PixelId{3}));
                static_cast<void>(invalidSeed.sort(image));
            },
            "unpadded infinity pixel bounds");
    }

    {
        auto plateau = makeImage(2, 2,
                                 {
                                     0,
                                     0,
                                     2,
                                     2,
                                 });

        for (TestTopographicImmersion interpolation : {TestTopographicImmersion::SelfDualSpan, TestTopographicImmersion::Min4Max8, TestTopographicImmersion::Min8Max4}) {
            TreeOfShapesProducer builder(makeTopographicConvention(plateau, interpolation));
            auto result = builder.build<ToSGrayLevel>(plateau);
            int numRoots = 0;
            for (NodeId nodeId = 0; nodeId < static_cast<NodeId>(result.parent.size()); ++nodeId) {
                const NodeId parent = result.parent[static_cast<std::size_t>(nodeId)];
                require(parent != InvalidNode, "TreeOfShapesProducer native parent through interpolated plateau must be valid");
                if (parent == nodeId) {
                    ++numRoots;
                }
            }
            requireEqual(numRoots, 1, "TreeOfShapesProducer native plateau root count");
            requireEqual(static_cast<int>(result.smallestNodeMap.size()), 4, "TreeOfShapesProducer native plateau proper-part count");
            require(result.parent.size() >= 2, "TreeOfShapesProducer native plateau node count");
        }
    }

    {
        auto image = makeComponentTreeFixture();
        auto tree = makeComponentTree(image, true);
        auto parent = exportFlatHigraHierarchy(*tree).first;
        auto rebuilt = makeTreeFromHigraParent(parent, image->getNumRows(), image->getNumColumns(), true);
        const NodeId rebuiltRoot = rebuilt->root();

        requireEqual(rebuilt->parent(rebuiltRoot), rebuiltRoot, "Higra rebuild root parent");
        requireEqual(rebuilt->numNodes(), 6, "Higra rebuild numNodes");
        requireEqual(computeAreaViaAttributeFacade(*rebuilt, rebuiltRoot), 16, "Higra rebuild root area");
        const NodeId smallestNodeOfPixel10 = rebuilt->smallestNode(10);
        require(rebuilt->isAlive(smallestNodeOfPixel10), "Higra rebuild smallest node must be alive");
        require(rebuilt->parent(smallestNodeOfPixel10) != InvalidNode, "Higra rebuild smallestNodeId must remain attached");
        requireEqual(static_cast<int>(rebuilt->leaves().size()), 1, "Higra rebuild leaf count");
    }

    {
        auto image = makeComponentTreeFixture();
        auto maxTree = makeComponentTree(image, true);
        auto minTree = makeComponentTree(image, false);

        requireEqual(maxTree->lowestCommonAncestor(5, 5), 5, "LCA self");
        requireEqual(maxTree->lowestCommonAncestor(0, 5), 0, "LCA root and leaf");
        requireEqual(maxTree->lowestCommonAncestor(4, 5), 4, "LCA ancestor and descendant");
        requireEqual(maxTree->lowestCommonAncestor(5, 4), 4, "LCA symmetry on chain");
        requireEqual(minTree->lowestCommonAncestor(3, 4), 2, "LCA siblings in min-tree");
        requireEqual(minTree->lowestCommonAncestor(4, 5), 2, "LCA cousins in min-tree");
        requireSharedLcaMatchesParentPaths(*maxTree);
        requireSharedLcaMatchesParentPaths(*minTree);

        maxTree->mergeNodeIntoParent(5);
        requireEqual(maxTree->lowestCommonAncestor(5, 4), InvalidNode, "removed LCA endpoint");
        requireSharedLcaMatchesParentPaths(*maxTree);
        auto editor = maxTree->edit();
        const NodeId reused = editor.createDetachedNode();
        editor.attach(4, reused);
        editor.reparent(reused, 2);
        requireEqual(maxTree->lowestCommonAncestor(4, reused), 2, "LCA after moveNode");
        requireEqual(maxTree->lowestCommonAncestor(reused, 3), 2, "LCA after reattachment");
        requireSharedLcaMatchesParentPaths(*maxTree);
    }

    {
        constexpr NodeId branchLength = 32768;
        constexpr int numPixels = 2;
        constexpr NodeId firstBranchLeaf = branchLength;
        constexpr NodeId secondBranchRoot = branchLength + 1;
        constexpr NodeId secondBranchLeaf = branchLength * 2;
        constexpr int numNodeSlots = secondBranchLeaf + 1;

        auto higraNodeId = [](NodeId slotId) { return numPixels + slotId; };

        std::vector<NodeId> parent(static_cast<size_t>(numPixels + numNodeSlots), InvalidNode);
        parent[0] = higraNodeId(firstBranchLeaf);
        parent[1] = higraNodeId(secondBranchLeaf);
        parent[static_cast<size_t>(higraNodeId(0))] = higraNodeId(0);
        for (NodeId nodeId = 1; nodeId <= firstBranchLeaf; ++nodeId) {
            parent[static_cast<size_t>(higraNodeId(nodeId))] = higraNodeId(nodeId - 1);
        }
        parent[static_cast<size_t>(higraNodeId(secondBranchRoot))] = higraNodeId(0);
        for (NodeId nodeId = secondBranchRoot + 1; nodeId <= secondBranchLeaf; ++nodeId) {
            parent[static_cast<size_t>(higraNodeId(nodeId))] = higraNodeId(nodeId - 1);
        }

        const auto altitude = makeStrictHigraAltitude<std::uint32_t>(parent, numPixels, true);
        auto deepValuedTree = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(parent), std::span<const std::uint32_t>(altitude), 1, 2,
                                                                                MorphologicalTreeKind::MaxTree, RegularGridAdjacency2D(1, 2, 1.5));
        const auto& deepTree = deepValuedTree.topology();
        requireEqual(deepTree.numDescendants(0), numNodeSlots - 1, "deep tree descendant count");
        requireEqual(deepTree.lowestCommonAncestor(firstBranchLeaf, secondBranchLeaf), 0, "deep tree LCA across branches");
    }

    {
        const std::vector<NodeId> nodeParent{0, 0, 0};
        const std::vector<NodeId> smallestNodeMap{1, 1, 2, 2};
        const std::vector<std::uint8_t> altitude{0, 1, 0};
        auto valuedTree = MorphologicalTreeFactory::createFromNativeTopology(
            std::span<const NodeId>(nodeParent), std::span<const NodeId>(smallestNodeMap), std::span<const std::uint8_t>(altitude), 0, 2, 2,
            MorphologicalTreeSemantics{
                MorphologicalTreeKind::TreeOfShapes, NodeAltitudeOrder::Unconstrained,
                TopographicConvention{ComplementaryGridImmersion{ComplementaryAdjacencies{
                                          RegularGridAdjacency2D(2, 2, 1.0), RegularGridAdjacency2D(2, 2, 1.5)}},
                                      TopographicDomainExtension::ExteriorRing, PixelId{0}}});
        const MorphologicalTree& tree = valuedTree.topology();
        requireEqual(tree.properPartCardinality(tree.root()), 0, "connected-subset root proper-part cardinality");
        requireEqual(tree.numChildren(tree.root()), 2, "connected-subset root children");
        const auto* topographicConvention = tree.topographicConvention();
        require(topographicConvention != nullptr, "connected-subset topographic convention");
        const auto* complementaryImmersion = std::get_if<ComplementaryGridImmersion>(&topographicConvention->immersion);
        require(complementaryImmersion != nullptr, "connected-subset complementary-grid immersion");
        require(tree.nodeAltitudeOrder() == NodeAltitudeOrder::Unconstrained, "connected-subset altitude order");
        require(tree.hasEmptyProperPart(tree.root()), "connected-subset root has an empty proper part");
        require(!tree.isTreeOfPartialPartitions(), "empty root proper part must reject the partial-partitions specialization");
        requireThrows([&tree] { tree.validateTreeOfPartialPartitions(); }, "partial-partitions validation must reject an empty proper part");
        requireEqual(complementaryImmersion->complementaryAdjacencies.minAdjacency.getSize(), 5, "connected-subset minimum adjacency");
        requireEqual(complementaryImmersion->complementaryAdjacencies.maxAdjacency.getSize(), 9, "connected-subset maximum adjacency");
        requireEqual(computeAreaViaAttributeFacade(tree, tree.root()), 4, "connected-subset root area");
        requireVectorEqual(collectImageValues(valuedTree.reconstructFromNodeAltitudes()), std::vector<std::uint8_t>{1, 1, 0, 0},
                           "connected-subset reconstruction");
    }

    {
        const std::vector<NodeId> nodeParent{0, 0, 1};
        const std::vector<NodeId> smallestNodeMap{2};
        const std::vector<std::uint8_t> altitude{0, 1, 2};
        auto valuedTree = MorphologicalTreeFactory::createFromNativeTopology(
            std::span<const NodeId>(nodeParent), std::span<const NodeId>(smallestNodeMap), std::span<const std::uint8_t>(altitude), 0, 1, 1,
            MorphologicalTreeSemantics{MorphologicalTreeKind::Generic, NodeAltitudeOrder::Increasing, NoConstructionContext{}});
        const MorphologicalTree& tree = valuedTree.topology();
        require(tree.kind() == MorphologicalTreeKind::Generic, "generic native topology descriptive kind");
        require(tree.nodeAltitudeOrder() == NodeAltitudeOrder::Increasing, "generic native topology altitude order");
        require(std::holds_alternative<NoConstructionContext>(tree.constructionContext()), "generic native topology has no construction context");
        require(tree.hasGridDomain2D() && tree.gridDomain2D()->rows == 1 && tree.gridDomain2D()->columns == 1, "generic native topology explicit 2D domain");
        require(tree.hasEmptyProperPart(0) && tree.hasEmptyProperPart(1) && !tree.hasEmptyProperPart(2),
                "generic native topology derives structural nodes from pixel assignments to smallest nodes");
        require(!tree.isTreeOfPartialPartitions(), "generic chain with empty proper parts is broader than a tree of partial partitions");
        requireThrows([&tree] { tree.validateTreeOfPartialPartitions(); }, "tree-of-partial-partitions validation must reject the generic chain");
        for (NodeId nodeId : tree.aliveNodeIds()) {
            require(!collectNodeIds(tree.nodeSupport(nodeId)).empty(), "every committed generic node must have non-empty subtree support");
        }
        valuedTree.setNodeAltitudes(std::vector<std::uint8_t>{0, 1, 3});
        requireThrows([&valuedTree] { valuedTree.setNodeAltitudes(std::vector<std::uint8_t>{0, 3, 2}); },
                      "generic hierarchy altitude validation must use the declared order");
        requireThrows([&valuedTree] { valuedTree.setNodeAltitudes(std::vector<std::uint8_t>{0, 1, 1}); },
                      "generic hierarchy altitude validation must reject equality");
    }

    {
        const std::vector<NodeId> nodeParent{0, 0, 0};
        const std::vector<NodeId> smallestNodeMap{1, 2};
        const std::vector<std::uint8_t> altitude{10, 3, 20};
        auto valuedTree = MorphologicalTreeFactory::createFromNativeTopology(std::span<const NodeId>(nodeParent), std::span<const NodeId>(smallestNodeMap),
                                                                           std::span<const std::uint8_t>(altitude), 0, MorphologicalTreeSemantics{});
        const MorphologicalTree& tree = valuedTree.topology();

        require(!tree.hasGridDomain2D(), "abstract native topology must not invent a 2D domain");
        require(!tree.gridDomain2D().has_value(), "abstract native topology optional grid");
        requireEqual(tree.numPixels(), 2, "abstract native topology proper-part cardinality");
        requireEqual(computeAreaViaAttributeFacade(tree, tree.root()), 2, "support attributes must work without grid metadata");

        auto [grayNames, grayBuffer] = AttributeComputation::computeSingleAttribute(valuedTree, GrayLevelHeight);
        requireNear(grayBuffer[grayNames.linearIndex(0, GrayLevelHeight)], 10.0f, 1.0e-6f, "unconstrained gray height must use the farthest subtree altitude");
        requireNear(grayBuffer[grayNames.linearIndex(1, GrayLevelHeight)], 0.0f, 1.0e-6f, "unconstrained gray height leaf one");
        requireNear(grayBuffer[grayNames.linearIndex(2, GrayLevelHeight)], 0.0f, 1.0e-6f, "unconstrained gray height leaf two");

        requireThrows([&valuedTree] { static_cast<void>(valuedTree.reconstructFromNodeAltitudes()); }, "image reconstruction must require explicit 2D domain metadata");
        if constexpr (contract::validationsEnabled) {
            requireThrows([&valuedTree] { static_cast<void>(AttributeComputation::computeSingleAttribute(valuedTree, BoxWidth)); },
                          "geometric attributes must require explicit 2D domain metadata");

            bool namedCapabilityDiagnostic = false;
            try {
                static_cast<void>(AttributeComputation::computeSingleAttribute(valuedTree, BoxWidth));
            } catch (const std::invalid_argument& error) {
                const std::string message = error.what();
                namedCapabilityDiagnostic =
                    message.find("BOX_WIDTH") != std::string::npos && message.find("regular 2D pixel domain") != std::string::npos;
            }
            require(namedCapabilityDiagnostic, "capability diagnostics must name the attribute and missing contract");
        }
    }

    {
        std::mt19937 random(0x5A17u);
        for (int trial = 0; trial < 32; ++trial) {
            const int numNodes = 2 + (trial * 7) % 127;
            std::vector<NodeId> nodeParent(static_cast<std::size_t>(numNodes), 0);
            std::vector<int> childCount(static_cast<std::size_t>(numNodes), 0);
            nodeParent[0] = 0;
            for (NodeId nodeId = 1; nodeId < numNodes; ++nodeId) {
                const NodeId parentId = static_cast<NodeId>(random() % static_cast<unsigned>(nodeId));
                nodeParent[static_cast<std::size_t>(nodeId)] = parentId;
                ++childCount[static_cast<std::size_t>(parentId)];
            }

            std::vector<NodeId> smallestNodeMap;
            for (NodeId nodeId = 0; nodeId < numNodes; ++nodeId) {
                if (childCount[static_cast<std::size_t>(nodeId)] == 0) {
                    smallestNodeMap.push_back(nodeId);
                }
            }

            std::vector<std::uint8_t> altitude(static_cast<std::size_t>(numNodes));
            std::vector<std::uint8_t> invertedAltitude(static_cast<std::size_t>(numNodes));
            for (NodeId nodeId = 0; nodeId < numNodes; ++nodeId) {
                const auto value = static_cast<std::uint8_t>(random() % 256u);
                altitude[static_cast<std::size_t>(nodeId)] = value;
                invertedAltitude[static_cast<std::size_t>(nodeId)] = static_cast<std::uint8_t>(255u - value);
            }

            auto valuedTree = MorphologicalTreeFactory::createFromNativeHierarchy(
                NativeHierarchyView<std::uint8_t>{nodeParent, smallestNodeMap, altitude, 0, std::nullopt, MorphologicalTreeSemantics{}});
            auto inverted = MorphologicalTreeFactory::createFromNativeHierarchy(
                NativeHierarchyView<std::uint8_t>{nodeParent, smallestNodeMap, invertedAltitude, 0, std::nullopt, MorphologicalTreeSemantics{}});

            const MorphologicalTree& tree = valuedTree.topology();
            requireEqual(computeAreaViaAttributeFacade(tree, tree.root()), static_cast<int>(smallestNodeMap.size()),
                         "random generic hierarchy root support");

            auto [grayNames, grayBuffer] = AttributeComputation::computeSingleAttribute(valuedTree, GrayLevelHeight);
            auto [invertedNames, invertedBuffer] = AttributeComputation::computeSingleAttribute(inverted, GrayLevelHeight);

            for (NodeId nodeId : tree.aliveNodeIds()) {
                int oracle = 0;
                for (NodeId descendant : tree.subtreeNodes(nodeId)) {
                    oracle = std::max(oracle, std::abs(static_cast<int>(altitude[static_cast<std::size_t>(nodeId)]) -
                                                       static_cast<int>(altitude[static_cast<std::size_t>(descendant)])));
                }
                requireNear(grayBuffer[grayNames.linearIndex(nodeId, GrayLevelHeight)], static_cast<float>(oracle), 1.0e-6f, "random generic GrayLevelHeight oracle");
                requireNear(invertedBuffer[invertedNames.linearIndex(nodeId, GrayLevelHeight)], static_cast<float>(oracle), 1.0e-6f,
                            "GrayLevelHeight contrast inversion");
            }
        }
    }

    {
        const std::vector<NodeId> nodeParent{0, 0};
        const std::vector<NodeId> smallestNodeMap{0, 1};
        const std::vector<std::uint8_t> altitude{0, 1};
        auto capabilityDriven = MorphologicalTreeFactory::createFromNativeTopology(
            std::span<const NodeId>(nodeParent), std::span<const NodeId>(smallestNodeMap), std::span<const std::uint8_t>(altitude), 0, 1, 2,
            MorphologicalTreeSemantics{MorphologicalTreeKind::Generic, NodeAltitudeOrder::Increasing,
                                       SharedAdjacencyContext{RegularGridAdjacency2D(1, 2, 1.5)}});
        require(capabilityDriven.topology().kind() == MorphologicalTreeKind::Generic, "capability test generic kind");
        auto [maxDistNames, maxDistBuffer] = AttributeComputation::computeSingleAttribute(capabilityDriven, MaxDist);
        require(std::isfinite(maxDistBuffer[maxDistNames.linearIndex(0, MaxDist)]) && std::isfinite(maxDistBuffer[maxDistNames.linearIndex(1, MaxDist)]),
                "MAX_DIST must depend on capabilities rather than descriptive kind");

        if constexpr (contract::validationsEnabled) {
            requireThrows(
                [&] {
                    static_cast<void>(MorphologicalTreeFactory::createFromNativeTopology(
                        std::span<const NodeId>(nodeParent), std::span<const NodeId>(smallestNodeMap), std::span<const std::uint8_t>(altitude), 0, 1, 2,
                        MorphologicalTreeSemantics{MorphologicalTreeKind::MaxTree, NodeAltitudeOrder::Unconstrained,
                                                   SharedAdjacencyContext{RegularGridAdjacency2D(1, 2, 1.5)}}));
                },
                "native import must reject a max-tree kind with an incompatible altitude order");
        }
    }

    {
        requireThrows(
            [] {
                const std::vector<NodeId> nodeParent{0, 0};
                const std::vector<NodeId> smallestNodeMap{1, 2};
                const std::vector<std::uint8_t> altitude{0, 1};
                static_cast<void>(MorphologicalTreeFactory::createFromNativeTopology(std::span<const NodeId>(nodeParent),
                                                                                     std::span<const NodeId>(smallestNodeMap),
                                                                                     std::span<const std::uint8_t>(altitude), 0, 1, 2, MorphologicalTreeSemantics{}));
            },
            "native topology must reject an smallestNodeId outside the node domain");
        requireThrows(
            [] {
                const std::vector<NodeId> nodeParent{0, 2, 1};
                const std::vector<NodeId> smallestNodeMap{1, 2};
                const std::vector<std::uint8_t> altitude{0, 1, 2};
                static_cast<void>(MorphologicalTreeFactory::createFromNativeTopology(std::span<const NodeId>(nodeParent),
                                                                                     std::span<const NodeId>(smallestNodeMap),
                                                                                     std::span<const std::uint8_t>(altitude), 0, 1, 2, MorphologicalTreeSemantics{}));
            },
            "native topology must reject a hierarchy disconnected from root");
        requireThrows(
            [] {
                const std::vector<NodeId> nodeParent{0, 0};
                const std::vector<NodeId> smallestNodeMap{0, 1};
                const std::vector<std::uint8_t> altitude{0};
                static_cast<void>(MorphologicalTreeFactory::createFromNativeTopology(std::span<const NodeId>(nodeParent),
                                                                                     std::span<const NodeId>(smallestNodeMap),
                                                                                     std::span<const std::uint8_t>(altitude), 0, 1, 2, MorphologicalTreeSemantics{}));
            },
            "native topology must reject an altitude shape mismatch");
        requireThrows(
            [] {
                const std::vector<NodeId> nodeParent{0, 0};
                const std::vector<NodeId> smallestNodeMap{0};
                const std::vector<std::uint8_t> altitude{0, 1};
                static_cast<void>(MorphologicalTreeFactory::createFromNativeTopology(std::span<const NodeId>(nodeParent),
                                                                                     std::span<const NodeId>(smallestNodeMap),
                                                                                     std::span<const std::uint8_t>(altitude), 0, 1, 1, MorphologicalTreeSemantics{}));
            },
            "native topology must reject a live leaf with empty subtree support");
        requireThrows(
            [] {
                const std::vector<NodeId> nodeParent{0, 0};
                const std::vector<NodeId> smallestNodeMap{0, 1};
                const std::vector<std::uint8_t> altitude{0, 1};
                static_cast<void>(MorphologicalTreeFactory::createFromNativeTopology(
                    std::span<const NodeId>(nodeParent), std::span<const NodeId>(smallestNodeMap), std::span<const std::uint8_t>(altitude), 0,
                    MorphologicalTreeSemantics{MorphologicalTreeKind::Generic, NodeAltitudeOrder::Increasing,
                                               SharedAdjacencyContext{RegularGridAdjacency2D(1, 2, 1.5)}}));
            },
            "abstract pixel domain must reject grid adjacency without grid metadata");
        requireThrows(
            [] {
                const std::vector<NodeId> nodeParent{0, 0, 1};
                const std::vector<NodeId> smallestNodeMap{2};
                const std::vector<std::uint8_t> altitude{0, 2, 1};
                static_cast<void>(MorphologicalTreeFactory::createFromNativeTopology(
                    std::span<const NodeId>(nodeParent), std::span<const NodeId>(smallestNodeMap), std::span<const std::uint8_t>(altitude), 0, 1, 1,
                    MorphologicalTreeSemantics{MorphologicalTreeKind::Generic, NodeAltitudeOrder::Increasing, NoConstructionContext{}}));
            },
            "native topology altitude input must satisfy the declared generic order");
    }

    {
        requireThrows(
            [] {
                std::vector<NodeId> noAdjacencyParent = {1, 1};
                std::vector<std::uint8_t> altitude(noAdjacencyParent.size(), std::uint8_t{});
                auto tree = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(noAdjacencyParent), std::span<const std::uint8_t>(altitude),
                                                                            1, 1, MorphologicalTreeKind::MaxTree);
                static_cast<void>(tree);
            },
            "Higra max-tree import must reject missing adjacency");

        {
            std::vector<NodeId> tosParent = {1, 1};
            std::vector<std::uint8_t> altitude(tosParent.size(), std::uint8_t{});
            auto tree = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(tosParent), std::span<const std::uint8_t>(altitude), 1, 1,
                                                                        MorphologicalTreeKind::TreeOfShapes);
            require(std::holds_alternative<NoConstructionContext>(tree.topology().constructionContext()),
                    "Higra tree-of-shapes import records unavailable construction context explicitly");
        }

        requireThrows(
            [] {
                std::vector<NodeId> badOwner = {0, 2, 2};
                std::vector<std::uint8_t> altitude(badOwner.size(), std::uint8_t{});
                auto tree = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(badOwner), std::span<const std::uint8_t>(altitude), 1, 2,
                                                                            MorphologicalTreeKind::MaxTree, RegularGridAdjacency2D(1, 2, 1.5));
                static_cast<void>(tree);
            },
            "Higra import must reject invalid smallest nodes");

        requireThrows(
            [] {
                std::vector<NodeId> cyclicParent = {1, 2, 1};
                std::vector<std::uint8_t> altitude(cyclicParent.size(), std::uint8_t{});
                auto tree = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(cyclicParent), std::span<const std::uint8_t>(altitude), 1,
                                                                            1, MorphologicalTreeKind::MaxTree, RegularGridAdjacency2D(1, 1, 1.5));
                static_cast<void>(tree);
            },
            "Higra import must reject disconnected cycles");

        requireThrows(
            [] {
                std::vector<NodeId> noInternalNodes = {0, 0, 0, 0};
                std::vector<std::uint8_t> altitude(noInternalNodes.size(), std::uint8_t{});
                auto tree = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(noInternalNodes), std::span<const std::uint8_t>(altitude),
                                                                            2, 2, MorphologicalTreeKind::MaxTree, RegularGridAdjacency2D(2, 2, 1.5));
                static_cast<void>(tree);
            },
            "Higra import must reject missing internal nodes");

        requireThrows(
            [] {
                std::vector<NodeId> badHigraParent = {1, 2, 3, 3};
                std::vector<std::uint8_t> altitude = {0, 0, 0, 0};
                auto valuedTree =
                    MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(badHigraParent), std::span<const std::uint8_t>(altitude), 1, 2,
                                                                    MorphologicalTreeKind::MaxTree, RegularGridAdjacency2D(1, 2, 1.5));
                static_cast<void>(valuedTree);
            },
            "Higra import must reject leaves pointing outside the internal-node domain");

        requireThrows(
            [] {
                auto valuedTree = makeValuedComponentTree(makeComponentTreeFixture(), true);
                auto [higraParent, higraAltitude] = valuedTree->exportHigraHierarchy();
                higraAltitude.pop_back();
                auto imported =
                    MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(higraParent), std::span<const std::uint8_t>(higraAltitude), 4, 4,
                                                                    MorphologicalTreeKind::MaxTree, RegularGridAdjacency2D(4, 4, 1.5));
                static_cast<void>(imported);
            },
            "Higra import must reject altitude buffers with the wrong size");

        if constexpr (contract::validationsEnabled) {
            requireThrows(
                [] {
                    auto valuedTree = MorphologicalTreeFactory::createMaxTree(makeComponentTreeFixture());
                    static_cast<void>(valuedTree.nodeAltitude(999));
                },
                "valuedTree altitude access must reject invalid node ids");
        }
    }

    return 0;
}
