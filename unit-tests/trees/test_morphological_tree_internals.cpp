#include "support/TestSupport.hpp"

#include "mmcfilters/trees/TreeOfShapesProducer.hpp"
#include "mmcfilters/trees/detail/ComponentTreeUnionFind.hpp"
#include "mmcfilters/trees/TreeEditor.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <random>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

static_assert(sizeof(ToSFloodDepth) >= sizeof(std::uint32_t), "Tree-of-shapes flood depth must not wrap at 16 bits.");

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
        requireEqual(maxTree.topology().getNumNodes(), 1, "constant max-tree node count");
        requireEqual(minTree.topology().getNumNodes(), 1, "constant min-tree node count");
        requireVectorEqual(collectNodeIds(maxTree.topology().getProperParts(maxTree.topology().getRoot())), std::vector<NodeId>{0, 1, 2, 3},
                           "constant max-tree proper parts");
        requireVectorEqual(collectNodeIds(minTree.topology().getProperParts(minTree.topology().getRoot())), std::vector<NodeId>{0, 1, 2, 3},
                           "constant min-tree proper parts");

        auto weighted = MorphologicalTreeFactory::createMaxTree(image);
        requireVectorEqual(collectImageValues(weighted.reconstructionImage()), collectImageValues(image), "constant weighted reconstruction");
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
        for (int pixelId : ordered) {
            sortedValues.push_back((*image)[pixelId]);
        }
        require(std::is_sorted(sortedValues.begin(), sortedValues.end()), "ComponentTreeUnionFind float sort must be nondecreasing for max-tree");
    }

    {
        auto image = ImageFloat::create(1, 2);
        (*image)[0] = -0.0f;
        (*image)[1] = +0.0f;
        auto tree = MorphologicalTreeFactory::createMaxTree(image, 1.5);
        requireEqual(tree.topology().getNumNodes(), 1, "signed-zero plateau node count");
        require(std::signbit(tree.getAltitude(tree.topology().getRoot())),
                "component-tree altitude must retain the first row-major proper-part representation");
    }

    {
        detail::AdjacencyUC adjacency(3, 3, true);
        requireVectorEqual(collectNodeIds(adjacency.getNeighborIndices(1, 1)), std::vector<NodeId>{1, 3, 7, 5}, "AdjacencyUC default 4-neighborhood");
        adjacency.setDiagonalConnection(1, 1, detail::DiagonalConnection::SW);
        adjacency.setDiagonalConnection(1, 1, detail::DiagonalConnection::NE);
        adjacency.setDiagonalConnection(1, 1, detail::DiagonalConnection::SE);
        adjacency.setDiagonalConnection(1, 1, detail::DiagonalConnection::NW);
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

        TreeOfShapesProducer tosBuilder(ToSInterpolation::Min4cMax8c);
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
        for (int pixelId = 0; pixelId < 49; ++pixelId) {
            requireEqual(sortedInterpolatedOrder[static_cast<std::size_t>(pixelId)], pixelId, "TreeOfShapesProducer sorted order must be a permutation");
        }

        auto result = tosBuilder.build(image);
        requireEqual(static_cast<int>(result.nodeParent.size()), 6, "TreeOfShapesProducer native node count");
        requireEqual(static_cast<int>(result.properPartOwner.size()), 9, "TreeOfShapesProducer native proper-part count");
        for (NodeId owner : result.properPartOwner) {
            require(owner >= 0 && owner < static_cast<NodeId>(result.nodeParent.size()), "TreeOfShapesProducer proper-part owner must be valid");
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

        TreeOfShapesProducer min4max8(ToSInterpolation::Min4cMax8c);
        auto [min4max8Min, min4max8Max, min4max8Adj] = min4max8.interpolateImage4c8c(saddle);
        const NodeId center = ImageUtils::to1D(2, 2, 5);
        requireEqual(min4max8Min[center], static_cast<uint8_t>(8), "Min4cMax8c saddle center min");
        requireEqual(min4max8Max[center], static_cast<uint8_t>(9), "Min4cMax8c saddle center max");
        requireVectorEqual(collectNodeIds(min4max8Adj.getNeighborIndices(center)), std::vector<NodeId>{7, 11, 17, 13, 16, 8},
                           "Min4cMax8c saddle center adjacency");

        TreeOfShapesProducer min8max4(ToSInterpolation::Min8cMax4c);
        auto [min8max4Min, min8max4Max, min8max4Adj] = min8max4.interpolateImage4c8c(saddle);
        requireEqual(min8max4Min[center], static_cast<uint8_t>(0), "Min8cMax4c saddle center min");
        requireEqual(min8max4Max[center], static_cast<uint8_t>(4), "Min8cMax4c saddle center max");
        requireVectorEqual(collectNodeIds(min8max4Adj.getNeighborIndices(center)), std::vector<NodeId>{7, 11, 17, 13, 18, 6},
                           "Min8cMax4c saddle center adjacency");
    }

    {
        auto halfLevelBoundary = makeImage(1, 2, {0, 1});
        TreeOfShapesProducer selfDual(ToSInterpolation::SelfDual);
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

        TreeOfShapesProducer lowerRightInfinity(ToSInterpolation::SelfDual, 2, 4);
        auto [customImgU, customOrderInterpolated, __] = lowerRightInfinity.sort(halfLevelBoundary);
        const NodeId lowerRightCorner = ImageUtils::to1D(2, 4, 5);
        requireEqual(customOrderInterpolated.front(), lowerRightCorner, "SelfDual propagation must honor a custom infinity seed");
        requireEqual(customImgU[lowerRightCorner], static_cast<ToSFloodDepth>(1), "SelfDual custom infinity seed must preserve Z/2 half level");
        TreeOfShapesProducer internalInfinity(ToSInterpolation::SelfDual, 1, 1);
        auto [internalImgU, internalOrderInterpolated, ___] = internalInfinity.sort(halfLevelBoundary);
        const NodeId firstOriginalSeed = ImageUtils::to1D(1, 1, 5);
        requireEqual(internalOrderInterpolated.front(), firstOriginalSeed, "SelfDual propagation must honor an internal custom infinity seed");
        requireThrows([&]() { static_cast<void>(TreeOfShapesProducer(ToSInterpolation::SelfDual, 3, 5).sort(halfLevelBoundary)); },
                      "SelfDual infinity seed must stay inside the interpolated domain");
    }

    {
        auto image = makeImage(1, 2, {0, 1});
        TreeOfShapesProducer producer(TreeOfShapesProducerOptions{ToSInterpolation::SelfDual, ToSPaddingPolicy::NoPadding, 0, 0});
        const auto& options = producer.options();
        require(options.padding == ToSPaddingPolicy::NoPadding, "TreeOfShapesProducer must retain its padding policy");

        auto [interpolationMin, interpolationMax, adjacency] = producer.interpolateImage(image);
        (void)adjacency;
        requireEqual(interpolationMin.size(), std::size_t{3}, "unpadded self-dual interpolation size");
        requireEqual(interpolationMax.size(), std::size_t{3}, "unpadded self-dual interpolation max size");
        requireEqual(interpolationMin[0], static_cast<ToSGrayLevel>(0), "unpadded first original level");
        requireEqual(interpolationMax[2], static_cast<ToSGrayLevel>(2), "unpadded second original level");

        auto [treeLevel, order, _] = producer.sort(image);
        requireEqual(treeLevel.size(), std::size_t{3}, "unpadded propagation domain size");
        requireEqual(order.front(), 0, "unpadded infinity seed is expressed in the selected immersion domain");

        requireThrows(
            [&]() {
                TreeOfShapesProducer invalidSeed(TreeOfShapesProducerOptions{ToSInterpolation::SelfDual, ToSPaddingPolicy::NoPadding, 1, 0});
                static_cast<void>(invalidSeed.sort(image));
            },
            "unpadded infinity seed bounds");
    }

    {
        auto plateau = makeImage(2, 2,
                                 {
                                     0,
                                     0,
                                     2,
                                     2,
                                 });

        for (ToSInterpolation interpolation : {ToSInterpolation::SelfDual, ToSInterpolation::Min4cMax8c, ToSInterpolation::Min8cMax4c}) {
            TreeOfShapesProducer builder(interpolation);
            auto result = builder.build(plateau);
            int numRoots = 0;
            for (NodeId nodeId = 0; nodeId < static_cast<NodeId>(result.nodeParent.size()); ++nodeId) {
                const NodeId parent = result.nodeParent[static_cast<std::size_t>(nodeId)];
                require(parent != InvalidNode, "TreeOfShapesProducer native parent through interpolated plateau must be valid");
                if (parent == nodeId) {
                    ++numRoots;
                }
            }
            requireEqual(numRoots, 1, "TreeOfShapesProducer native plateau root count");
            requireEqual(static_cast<int>(result.properPartOwner.size()), 4, "TreeOfShapesProducer native plateau proper-part count");
            require(result.nodeParent.size() >= 2, "TreeOfShapesProducer native plateau node count");
        }
    }

    {
        auto image = makeComponentTreeFixture();
        auto tree = makeComponentTree(image, true);
        auto parent = exportFlatHigraHierarchy(*tree).first;
        auto rebuilt = makeTreeFromHigraParent(parent, image->getNumRows(), image->getNumCols(), true);
        const NodeId rebuiltRoot = rebuilt->getRoot();

        requireEqual(rebuilt->getNodeParent(rebuiltRoot), rebuiltRoot, "Higra rebuild root parent");
        requireEqual(rebuilt->getNumNodes(), 6, "Higra rebuild numNodes");
        requireEqual(computeAreaViaAttributeFacade(*rebuilt, rebuiltRoot), 16, "Higra rebuild root area");
        const NodeId ownerOfPixel10 = rebuilt->getProperPartOwner(10);
        require(rebuilt->isAlive(ownerOfPixel10), "Higra rebuild proper-part owner must be alive");
        require(rebuilt->getNodeParent(ownerOfPixel10) != InvalidNode, "Higra rebuild owner must remain attached");
        requireEqual(static_cast<int>(rebuilt->getLeaves().size()), 1, "Higra rebuild leaf count");
    }

    {
        auto image = makeComponentTreeFixture();
        auto maxTree = makeComponentTree(image, true);
        auto minTree = makeComponentTree(image, false);

        requireEqual(maxTree->getLowestCommonAncestor(5, 5), 5, "LCA self");
        requireEqual(maxTree->getLowestCommonAncestor(0, 5), 0, "LCA root and leaf");
        requireEqual(maxTree->getLowestCommonAncestor(4, 5), 4, "LCA ancestor and descendant");
        requireEqual(maxTree->getLowestCommonAncestor(5, 4), 4, "LCA symmetry on chain");
        requireEqual(minTree->getLowestCommonAncestor(3, 4), 2, "LCA siblings in min-tree");
        requireEqual(minTree->getLowestCommonAncestor(4, 5), 2, "LCA cousins in min-tree");

        maxTree->mergeNodeIntoParent(5);
        auto editor = maxTree->edit();
        const NodeId reused = editor.createDetachedNode();
        editor.attach(4, reused);
        editor.reparent(reused, 2);
        requireEqual(maxTree->getLowestCommonAncestor(4, reused), 2, "LCA after moveNode");
        requireEqual(maxTree->getLowestCommonAncestor(reused, 3), 2, "LCA after reattachment");
    }

    {
        constexpr NodeId branchLength = 2048;
        constexpr NodeId numProperParts = 2;
        constexpr NodeId firstBranchLeaf = branchLength;
        constexpr NodeId secondBranchRoot = branchLength + 1;
        constexpr NodeId secondBranchLeaf = branchLength * 2;
        constexpr NodeId numNodeSlots = secondBranchLeaf + 1;

        auto higraNodeId = [](NodeId slotId) { return numProperParts + slotId; };

        std::vector<NodeId> parent(static_cast<size_t>(numProperParts + numNodeSlots), InvalidNode);
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

        const auto altitude = makeStrictHigraAltitude<std::uint32_t>(parent, numProperParts, true);
        auto deepTreeWeighted = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(parent), std::span<const std::uint32_t>(altitude), 1, 2,
                                                                                MorphologicalTreeKind::MAX_TREE, RegularGridAdjacency2D(1, 2, 1.5));
        const auto& deepTree = deepTreeWeighted.topology();
        requireEqual(deepTree.getNodeNumDescendants(0), numNodeSlots - 1, "deep tree descendant count");
        requireEqual(deepTree.getLowestCommonAncestor(firstBranchLeaf, secondBranchLeaf), 0, "deep tree LCA across branches");
    }

    {
        const std::vector<NodeId> nodeParent{0, 0, 0};
        const std::vector<NodeId> properPartOwner{1, 1, 2, 2};
        const std::vector<std::uint8_t> altitude{0, 1, 0};
        auto weighted = MorphologicalTreeFactory::createFromNativeTopology(
            std::span<const NodeId>(nodeParent), std::span<const NodeId>(properPartOwner), std::span<const std::uint8_t>(altitude), 0, 2, 2,
            HierarchySemantics{MorphologicalTreeKind::GENERIC, AltitudeOrder::UNCONSTRAINED,
                               DirectionalGridAdjacency2D{RegularGridAdjacency2D(2, 2, 1.0), RegularGridAdjacency2D(2, 2, 1.5)}});
        const MorphologicalTree& tree = weighted.topology();
        requireEqual(tree.getNumProperParts(tree.getRoot()), 0, "native partial-partition root direct proper parts");
        requireEqual(tree.getNumChildren(tree.getRoot()), 2, "native partial-partition root children");
        require(tree.hasDirectionalGridAdjacency2D(), "native partial-partition directional adjacency capability");
        require(tree.getAdjacencyMode() == AdjacencyMode::DIRECTIONAL, "native partial-partition adjacency mode");
        require(tree.getAltitudeOrder() == AltitudeOrder::UNCONSTRAINED, "native partial-partition altitude order");
        require(tree.isStructuralNode(tree.getRoot()), "native partial-partition root structural role");
        requireEqual(tree.getDecreasingGridAdjacency2D()->getSize(), 5, "native partial-partition decreasing adjacency");
        requireEqual(tree.getIncreasingGridAdjacency2D()->getSize(), 9, "native partial-partition increasing adjacency");
        requireEqual(computeAreaViaAttributeFacade(tree, tree.getRoot()), 4, "native partial-partition root area");
        requireVectorEqual(collectImageValues(weighted.reconstructionImage()), std::vector<std::uint8_t>{1, 1, 0, 0},
                           "native partial-partition reconstruction");
    }

    {
        const std::vector<NodeId> nodeParent{0, 0, 1};
        const std::vector<NodeId> properPartOwner{2};
        const std::vector<std::uint8_t> altitude{0, 1, 2};
        auto weighted = MorphologicalTreeFactory::createFromNativeTopology(
            std::span<const NodeId>(nodeParent), std::span<const NodeId>(properPartOwner), std::span<const std::uint8_t>(altitude), 0, 1, 1,
            HierarchySemantics{MorphologicalTreeKind::GENERIC, AltitudeOrder::INCREASING_FROM_ROOT, NoAdjacency{}});
        const MorphologicalTree& tree = weighted.topology();
        require(tree.getDescriptiveKind() == MorphologicalTreeKind::GENERIC, "generic native topology descriptive kind");
        require(tree.getAltitudeOrder() == AltitudeOrder::INCREASING_FROM_ROOT, "generic native topology altitude order");
        require(tree.getAdjacencyMode() == AdjacencyMode::NONE, "generic native topology adjacency mode");
        require(tree.hasGridDomain2D() && tree.getGridDomain2D()->rows == 1 && tree.getGridDomain2D()->cols == 1, "generic native topology explicit 2D domain");
        require(tree.isStructuralNode(0) && tree.isStructuralNode(1) && !tree.isStructuralNode(2),
                "generic native topology derives structural nodes from direct ownership");
        for (NodeId nodeId : tree.getAliveNodeIds()) {
            require(!collectNodeIds(tree.getConnectedComponent(nodeId)).empty(), "every committed generic node must have non-empty subtree support");
        }
        weighted.setAltitudeBuffer(std::vector<std::uint8_t>{0, 1, 3});
        requireThrows([&weighted] { weighted.setAltitudeBuffer(std::vector<std::uint8_t>{0, 3, 2}); },
                      "generic hierarchy altitude validation must use the declared order");
        requireThrows([&weighted] { weighted.setAltitudeBuffer(std::vector<std::uint8_t>{0, 1, 1}); },
                      "generic hierarchy altitude validation must reject equality");
    }

    {
        const std::vector<NodeId> nodeParent{0, 0, 0};
        const std::vector<NodeId> properPartOwner{1, 2};
        const std::vector<std::uint8_t> altitude{10, 3, 20};
        auto weighted = MorphologicalTreeFactory::createFromNativeTopology(std::span<const NodeId>(nodeParent), std::span<const NodeId>(properPartOwner),
                                                                           std::span<const std::uint8_t>(altitude), 0, HierarchySemantics{});
        const MorphologicalTree& tree = weighted.topology();

        require(!tree.hasGridDomain2D(), "abstract native topology must not invent a 2D domain");
        require(!tree.getGridDomain2D().has_value(), "abstract native topology optional grid");
        requireEqual(tree.getNumTotalProperParts(), 2, "abstract native topology proper-part cardinality");
        requireEqual(computeAreaViaAttributeFacade(tree, tree.getRoot()), 2, "support attributes must work without grid metadata");

        auto [grayNames, grayBuffer] = AttributeComputation::computeSingleAttribute(weighted, GRAY_HEIGHT);
        requireNear(grayBuffer[grayNames.linearIndex(0, GRAY_HEIGHT)], 10.0f, 1.0e-6f, "unconstrained gray height must use the farthest subtree altitude");
        requireNear(grayBuffer[grayNames.linearIndex(1, GRAY_HEIGHT)], 0.0f, 1.0e-6f, "unconstrained gray height leaf one");
        requireNear(grayBuffer[grayNames.linearIndex(2, GRAY_HEIGHT)], 0.0f, 1.0e-6f, "unconstrained gray height leaf two");

        requireThrows([&weighted] { static_cast<void>(weighted.reconstructionImage()); }, "image reconstruction must require explicit 2D domain metadata");
        requireThrows([&weighted] { static_cast<void>(AttributeComputation::computeSingleAttribute(weighted, BOX_WIDTH)); },
                      "geometric attributes must require explicit 2D domain metadata");

        bool namedCapabilityDiagnostic = false;
        try {
            static_cast<void>(AttributeComputation::computeSingleAttribute(weighted, BOX_WIDTH));
        } catch (const std::invalid_argument& error) {
            const std::string message = error.what();
            namedCapabilityDiagnostic = message.find("BOX_WIDTH") != std::string::npos && message.find("regular 2D proper-part domain") != std::string::npos;
        }
        require(namedCapabilityDiagnostic, "capability diagnostics must name the attribute and missing contract");
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

            std::vector<NodeId> properPartOwner;
            for (NodeId nodeId = 0; nodeId < numNodes; ++nodeId) {
                if (childCount[static_cast<std::size_t>(nodeId)] == 0) {
                    properPartOwner.push_back(nodeId);
                }
            }

            std::vector<std::uint8_t> altitude(static_cast<std::size_t>(numNodes));
            std::vector<std::uint8_t> invertedAltitude(static_cast<std::size_t>(numNodes));
            for (NodeId nodeId = 0; nodeId < numNodes; ++nodeId) {
                const auto value = static_cast<std::uint8_t>(random() % 256u);
                altitude[static_cast<std::size_t>(nodeId)] = value;
                invertedAltitude[static_cast<std::size_t>(nodeId)] = static_cast<std::uint8_t>(255u - value);
            }

            auto weighted = MorphologicalTreeFactory::createFromNativeHierarchy(
                NativeHierarchyView<std::uint8_t>{nodeParent, properPartOwner, altitude, 0, std::nullopt, HierarchySemantics{}});
            auto inverted = MorphologicalTreeFactory::createFromNativeHierarchy(
                NativeHierarchyView<std::uint8_t>{nodeParent, properPartOwner, invertedAltitude, 0, std::nullopt, HierarchySemantics{}});

            const MorphologicalTree& tree = weighted.topology();
            requireEqual(computeAreaViaAttributeFacade(tree, tree.getRoot()), static_cast<int>(properPartOwner.size()),
                         "random generic hierarchy root support");

            auto [grayNames, grayBuffer] = AttributeComputation::computeSingleAttribute(weighted, GRAY_HEIGHT);
            auto [invertedNames, invertedBuffer] = AttributeComputation::computeSingleAttribute(inverted, GRAY_HEIGHT);

            for (NodeId nodeId : tree.getAliveNodeIds()) {
                int oracle = 0;
                for (NodeId descendant : tree.getNodeSubtree(nodeId)) {
                    oracle = std::max(oracle, std::abs(static_cast<int>(altitude[static_cast<std::size_t>(nodeId)]) -
                                                       static_cast<int>(altitude[static_cast<std::size_t>(descendant)])));
                }
                requireNear(grayBuffer[grayNames.linearIndex(nodeId, GRAY_HEIGHT)], static_cast<float>(oracle), 1.0e-6f, "random generic GRAY_HEIGHT oracle");
                requireNear(invertedBuffer[invertedNames.linearIndex(nodeId, GRAY_HEIGHT)], static_cast<float>(oracle), 1.0e-6f,
                            "GRAY_HEIGHT contrast inversion");
            }
        }
    }

    {
        const std::vector<NodeId> nodeParent{0, 0};
        const std::vector<NodeId> properPartOwner{0, 1};
        const std::vector<std::uint8_t> altitude{0, 1};
        auto capabilityDriven = MorphologicalTreeFactory::createFromNativeTopology(
            std::span<const NodeId>(nodeParent), std::span<const NodeId>(properPartOwner), std::span<const std::uint8_t>(altitude), 0, 1, 2,
            HierarchySemantics{MorphologicalTreeKind::TREE_OF_SHAPES, AltitudeOrder::INCREASING_FROM_ROOT,
                               UniformGridAdjacency2D{RegularGridAdjacency2D(1, 2, 1.5)}});
        require(capabilityDriven.topology().getDescriptiveKind() == MorphologicalTreeKind::TREE_OF_SHAPES, "capability test descriptive kind");
        auto [maxDistNames, maxDistBuffer] = AttributeComputation::computeSingleAttribute(capabilityDriven, MAX_DIST);
        require(std::isfinite(maxDistBuffer[maxDistNames.linearIndex(0, MAX_DIST)]) && std::isfinite(maxDistBuffer[maxDistNames.linearIndex(1, MAX_DIST)]),
                "MAX_DIST must depend on capabilities rather than descriptive kind");

        auto misleadingDescriptor = MorphologicalTreeFactory::createFromNativeTopology(
            std::span<const NodeId>(nodeParent), std::span<const NodeId>(properPartOwner), std::span<const std::uint8_t>(altitude), 0, 1, 2,
            HierarchySemantics{MorphologicalTreeKind::MAX_TREE, AltitudeOrder::UNCONSTRAINED, UniformGridAdjacency2D{RegularGridAdjacency2D(1, 2, 1.5)}});
        requireThrows([&misleadingDescriptor] { static_cast<void>(AttributeComputation::computeSingleAttribute(misleadingDescriptor, MAX_DIST)); },
                      "MAX_DIST must reject a missing altitude-order capability even when the descriptor says max-tree");
    }

    {
        requireThrows(
            [] {
                const std::vector<NodeId> nodeParent{0, 0};
                const std::vector<NodeId> properPartOwner{1, 2};
                const std::vector<std::uint8_t> altitude{0, 1};
                static_cast<void>(MorphologicalTreeFactory::createFromNativeTopology(std::span<const NodeId>(nodeParent),
                                                                                     std::span<const NodeId>(properPartOwner),
                                                                                     std::span<const std::uint8_t>(altitude), 0, 1, 2, HierarchySemantics{}));
            },
            "native topology must reject an owner outside the node domain");
        requireThrows(
            [] {
                const std::vector<NodeId> nodeParent{0, 2, 1};
                const std::vector<NodeId> properPartOwner{1, 2};
                const std::vector<std::uint8_t> altitude{0, 1, 2};
                static_cast<void>(MorphologicalTreeFactory::createFromNativeTopology(std::span<const NodeId>(nodeParent),
                                                                                     std::span<const NodeId>(properPartOwner),
                                                                                     std::span<const std::uint8_t>(altitude), 0, 1, 2, HierarchySemantics{}));
            },
            "native topology must reject a hierarchy disconnected from root");
        requireThrows(
            [] {
                const std::vector<NodeId> nodeParent{0, 0};
                const std::vector<NodeId> properPartOwner{0, 1};
                const std::vector<std::uint8_t> altitude{0};
                static_cast<void>(MorphologicalTreeFactory::createFromNativeTopology(std::span<const NodeId>(nodeParent),
                                                                                     std::span<const NodeId>(properPartOwner),
                                                                                     std::span<const std::uint8_t>(altitude), 0, 1, 2, HierarchySemantics{}));
            },
            "native topology must reject an altitude shape mismatch");
        requireThrows(
            [] {
                const std::vector<NodeId> nodeParent{0, 0};
                const std::vector<NodeId> properPartOwner{0};
                const std::vector<std::uint8_t> altitude{0, 1};
                static_cast<void>(MorphologicalTreeFactory::createFromNativeTopology(std::span<const NodeId>(nodeParent),
                                                                                     std::span<const NodeId>(properPartOwner),
                                                                                     std::span<const std::uint8_t>(altitude), 0, 1, 1, HierarchySemantics{}));
            },
            "native topology must reject a live leaf with empty subtree support");
        requireThrows(
            [] {
                const std::vector<NodeId> nodeParent{0, 0};
                const std::vector<NodeId> properPartOwner{0, 1};
                const std::vector<std::uint8_t> altitude{0, 1};
                static_cast<void>(MorphologicalTreeFactory::createFromNativeTopology(
                    std::span<const NodeId>(nodeParent), std::span<const NodeId>(properPartOwner), std::span<const std::uint8_t>(altitude), 0,
                    HierarchySemantics{MorphologicalTreeKind::GENERIC, AltitudeOrder::INCREASING_FROM_ROOT,
                                       UniformGridAdjacency2D{RegularGridAdjacency2D(1, 2, 1.5)}}));
            },
            "abstract proper-part domain must reject grid adjacency without grid metadata");
        requireThrows(
            [] {
                const std::vector<NodeId> nodeParent{0, 0, 1};
                const std::vector<NodeId> properPartOwner{2};
                const std::vector<std::uint8_t> altitude{0, 2, 1};
                static_cast<void>(MorphologicalTreeFactory::createFromNativeTopology(
                    std::span<const NodeId>(nodeParent), std::span<const NodeId>(properPartOwner), std::span<const std::uint8_t>(altitude), 0, 1, 1,
                    HierarchySemantics{MorphologicalTreeKind::GENERIC, AltitudeOrder::INCREASING_FROM_ROOT, NoAdjacency{}}));
            },
            "native topology altitude input must satisfy the declared generic order");
    }

    {
        requireThrows(
            [] {
                std::vector<NodeId> noAdjacencyParent = {1, 1};
                std::vector<std::uint8_t> altitude(noAdjacencyParent.size(), std::uint8_t{});
                auto tree = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(noAdjacencyParent), std::span<const std::uint8_t>(altitude),
                                                                            1, 1, MorphologicalTreeKind::MAX_TREE);
                static_cast<void>(tree);
            },
            "Higra max-tree import must reject missing adjacency");

        {
            std::vector<NodeId> tosParent = {1, 1};
            std::vector<std::uint8_t> altitude(tosParent.size(), std::uint8_t{});
            auto tree = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(tosParent), std::span<const std::uint8_t>(altitude), 1, 1,
                                                                        MorphologicalTreeKind::TREE_OF_SHAPES);
            require(!tree.topology().hasUniformGridAdjacency2D(), "Higra tree-of-shapes import may omit adjacency");
        }

        requireThrows(
            [] {
                std::vector<NodeId> badOwner = {0, 2, 2};
                std::vector<std::uint8_t> altitude(badOwner.size(), std::uint8_t{});
                auto tree = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(badOwner), std::span<const std::uint8_t>(altitude), 1, 2,
                                                                            MorphologicalTreeKind::MAX_TREE, RegularGridAdjacency2D(1, 2, 1.5));
                static_cast<void>(tree);
            },
            "Higra import must reject invalid proper-part owners");

        requireThrows(
            [] {
                std::vector<NodeId> cyclicParent = {1, 2, 1};
                std::vector<std::uint8_t> altitude(cyclicParent.size(), std::uint8_t{});
                auto tree = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(cyclicParent), std::span<const std::uint8_t>(altitude), 1,
                                                                            1, MorphologicalTreeKind::MAX_TREE, RegularGridAdjacency2D(1, 1, 1.5));
                static_cast<void>(tree);
            },
            "Higra import must reject disconnected cycles");

        requireThrows(
            [] {
                std::vector<NodeId> noInternalNodes = {0, 0, 0, 0};
                std::vector<std::uint8_t> altitude(noInternalNodes.size(), std::uint8_t{});
                auto tree = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(noInternalNodes), std::span<const std::uint8_t>(altitude),
                                                                            2, 2, MorphologicalTreeKind::MAX_TREE, RegularGridAdjacency2D(2, 2, 1.5));
                static_cast<void>(tree);
            },
            "Higra import must reject missing internal nodes");

        requireThrows(
            [] {
                std::vector<NodeId> badHigraParent = {1, 2, 3, 3};
                std::vector<std::uint8_t> altitude = {0, 0, 0, 0};
                auto weighted =
                    MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(badHigraParent), std::span<const std::uint8_t>(altitude), 1, 2,
                                                                    MorphologicalTreeKind::MAX_TREE, RegularGridAdjacency2D(1, 2, 1.5));
                static_cast<void>(weighted);
            },
            "Higra import must reject leaves pointing outside the internal-node domain");

        requireThrows(
            [] {
                auto weighted = makeWeightedComponentTree(makeComponentTreeFixture(), true);
                auto [higraParent, higraAltitude] = weighted->exportHigraHierarchy();
                higraAltitude.pop_back();
                auto imported =
                    MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(higraParent), std::span<const std::uint8_t>(higraAltitude), 4, 4,
                                                                    MorphologicalTreeKind::MAX_TREE, RegularGridAdjacency2D(4, 4, 1.5));
                static_cast<void>(imported);
            },
            "Higra import must reject altitude buffers with the wrong size");

        requireThrows(
            [] {
                auto weighted = MorphologicalTreeFactory::createMaxTree(makeComponentTreeFixture());
                static_cast<void>(weighted.getAltitude(999));
            },
            "weighted altitude access must reject invalid node ids");
    }

    return 0;
}
