#include "support/TestSupport.hpp"

#include "mmcfilters/trees/BuilderMorphologicalTreeByUnionFind.hpp"
#include "mmcfilters/trees/TreeEditor.hpp"

#include <algorithm>
#include <memory>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

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
        AdjacencyRelation adj(4, 4, 1.5);

        BuilderComponentTree maxBuilder(&adj, true);
        auto [parentMax, orderMax, numNodesMax] = maxBuilder.createTreeByUnionFind(image);
        requireEqual(numNodesMax, 6, "BuilderComponentTree max numNodes");
        requireEqual(static_cast<int>(orderMax.size()), 16, "BuilderComponentTree max order size");
        requireVectorEqual(
            orderMax,
            std::vector<int>{15, 8, 12, 13, 2, 3, 7, 11, 0, 1, 4, 5, 6, 9, 10, 14},
            "BuilderComponentTree max ordered pixels"
        );
        requireVectorEqual(
            parentMax,
            std::vector<int>{2, 0, 8, 2, 0, 0, 5, 2, 15, 5, 5, 2, 8, 8, 10, 15},
            "BuilderComponentTree max parent array"
        );

        BuilderComponentTree minBuilder(&adj, false);
        auto [parentMin, orderMin, numNodesMin] = minBuilder.createTreeByUnionFind(image);
        requireEqual(numNodesMin, 6, "BuilderComponentTree min numNodes");
        requireVectorEqual(
            orderMin,
            std::vector<int>{10, 14, 5, 6, 9, 0, 1, 4, 2, 3, 7, 11, 8, 12, 13, 15},
            "BuilderComponentTree min ordered pixels"
        );
        requireVectorEqual(
            parentMin,
            std::vector<int>{5, 0, 0, 2, 0, 10, 5, 2, 0, 5, 10, 2, 8, 8, 10, 2},
            "BuilderComponentTree min parent array"
        );
    }

    {
        auto image = makeImage(2, 2, {0, 0, 0, 0});
        AdjacencyRelation adj(2, 2, 1.5);

        BuilderComponentTree maxBuilder(&adj, true);
        auto [parentMax, orderMax, numNodesMax] = maxBuilder.createTreeByUnionFind(image);
        requireEqual(numNodesMax, 1, "BuilderComponentTree constant max numNodes");
        requireEqual(static_cast<int>(orderMax.size()), 4, "BuilderComponentTree constant max order size");

        BuilderComponentTree minBuilder(&adj, false);
        auto [parentMin, orderMin, numNodesMin] = minBuilder.createTreeByUnionFind(image);
        requireEqual(numNodesMin, 1, "BuilderComponentTree constant min numNodes");
        requireEqual(static_cast<int>(orderMin.size()), 4, "BuilderComponentTree constant min order size");

        auto maxTree = MorphologicalTree::createComponentTree(image, true);
        auto minTree = MorphologicalTree::createComponentTree(image, false);
        requireEqual(maxTree.getNumNodes(), 1, "constant max-tree node count");
        requireEqual(minTree.getNumNodes(), 1, "constant min-tree node count");
        requireVectorEqual(collectNodeIds(maxTree.getProperParts(maxTree.getRoot())), std::vector<NodeId>{0, 1, 2, 3}, "constant max-tree proper parts");
        requireVectorEqual(collectNodeIds(minTree.getProperParts(minTree.getRoot())), std::vector<NodeId>{0, 1, 2, 3}, "constant min-tree proper parts");

        auto weighted = WeightedMorphologicalTree::createComponentTree(image, true);
        requireVectorEqual(collectImageValues(weighted.reconstructionImage()), collectImageValues(image), "constant weighted reconstruction");
    }

    {
        auto image = ImageFloat::create(2, 3);
        float values[] = {2.0f, 1.0f, 2.0f, 0.5f, 3.0f, 1.5f};
        for (int i = 0; i < 6; ++i) {
            (*image)[i] = values[i];
        }
        AdjacencyRelation adj(2, 3, 1.5);
        BuilderComponentTree maxBuilder(&adj, true);
        auto ordered = maxBuilder.sort(image);
        requireEqual(static_cast<int>(ordered.size()), 6, "BuilderComponentTree float sort size");
        std::vector<float> sortedValues;
        for (int pixelId : ordered) {
            sortedValues.push_back((*image)[pixelId]);
        }
        require(std::is_sorted(sortedValues.begin(), sortedValues.end()), "BuilderComponentTree float sort must be nondecreasing for max-tree");
    }

    {
        AdjacencyUC adjacency(3, 3, true);
        requireVectorEqual(collectNodeIds(adjacency.getNeighborPixels(1, 1)), std::vector<NodeId>{1, 3, 7, 5}, "AdjacencyUC default 4-neighborhood");
        adjacency.setDiagonalConnection(1, 1, DiagonalConnection::SW);
        adjacency.setDiagonalConnection(1, 1, DiagonalConnection::NE);
        adjacency.setDiagonalConnection(1, 1, DiagonalConnection::SE);
        adjacency.setDiagonalConnection(1, 1, DiagonalConnection::NW);
        requireVectorEqual(collectNodeIds(adjacency.getNeighborPixels(1, 1)), std::vector<NodeId>{1, 3, 7, 5, 6, 2, 8, 0}, "AdjacencyUC diagonal neighborhood");

        PriorityQueueToS queue;
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
    }

    {
        auto image = makeImage(
            3,
            3,
            {
                1, 2, 1,
                2, 3, 2,
                1, 2, 1,
            }
        );

        BuilderTreeOfShape tosBuilder(true);
        auto [interpolationMin, interpolationMax, adjacency] = tosBuilder.interpolateImage4c8c(image);
        requireEqual(static_cast<int>(interpolationMin.size()), 49, "BuilderTreeOfShape interpolate min size");
        requireEqual(static_cast<int>(interpolationMax.size()), 49, "BuilderTreeOfShape interpolate max size");
        requireEqual(interpolationMin[ImageUtils::to1D(3, 3, 7)], 3, "BuilderTreeOfShape interpolated center min");
        requireEqual(interpolationMax[ImageUtils::to1D(3, 3, 7)], 3, "BuilderTreeOfShape interpolated center max");
        requireVectorEqual(
            collectNodeIds(adjacency.getNeighborPixels(ImageUtils::to1D(3, 3, 7))),
            std::vector<NodeId>{17, 23, 31, 25, 32, 16},
            "BuilderTreeOfShape interpolated center adjacency"
        );

        auto [imgU, orderInterpolated, _] = tosBuilder.sort(image);
        requireEqual(static_cast<int>(imgU.size()), 49, "BuilderTreeOfShape sorted interpolation size");
        requireEqual(static_cast<int>(orderInterpolated.size()), 49, "BuilderTreeOfShape sorted order size");
        requireEqual(orderInterpolated.front(), 8, "BuilderTreeOfShape sorted first interpolated pixel");
        requireEqual(orderInterpolated.back(), 24, "BuilderTreeOfShape sorted last interpolated pixel");

        auto [parent, orderedPixels, numNodes] = tosBuilder.createTreeByUnionFind(image);
        requireEqual(numNodes, 6, "BuilderTreeOfShape numNodes");
        requireVectorEqual(orderedPixels, std::vector<int>{0, 3, 1, 7, 5, 6, 2, 8, 4}, "BuilderTreeOfShape ordered original pixels");
        requireVectorEqual(parent, std::vector<int>{0, 3, 3, 0, 3, 3, 3, 3, 3}, "BuilderTreeOfShape parent array");
    }

    {
        auto image = makeComponentTreeFixture();
        auto tree = makeComponentTree(image, true);
        auto parent = exportFlatHigraHierarchy(*tree).first;
        auto rebuilt = makeTreeFromHigraParent(parent, image->getNumRows(), image->getNumCols(), true);
        const NodeId rebuiltRoot = rebuilt->getRoot();

        requireEqual(rebuilt->getNodeParent(rebuiltRoot), rebuiltRoot, "Higra rebuild root parent");
        requireEqual(rebuilt->getNumNodes(), 6, "Higra rebuild numNodes");
        requireEqual(computeAreaAttribute(*rebuilt, rebuiltRoot), 16, "Higra rebuild root area");
        const NodeId ownerOfPixel10 = rebuilt->getSmallestComponent(10);
        require(rebuilt->isAlive(ownerOfPixel10), "Higra rebuild proper-part owner must be alive");
        require(rebuilt->getNodeParent(ownerOfPixel10) != InvalidNode, "Higra rebuild owner must remain attached");
        requireEqual(static_cast<int>(rebuilt->getLeaves().size()), 1, "Higra rebuild leaf count");
    }

    {
        auto image = makeComponentTreeFixture();
        auto maxTree = makeComponentTree(image, true);
        auto minTree = makeComponentTree(image, false);
        TreeEditor editor(*maxTree);

        requireEqual(maxTree->getLowestCommonAncestor(5, 5), 5, "LCA self");
        requireEqual(maxTree->getLowestCommonAncestor(0, 5), 0, "LCA root and leaf");
        requireEqual(maxTree->getLowestCommonAncestor(4, 5), 4, "LCA ancestor and descendant");
        requireEqual(maxTree->getLowestCommonAncestor(5, 4), 4, "LCA symmetry on chain");
        requireEqual(minTree->getLowestCommonAncestor(3, 4), 2, "LCA siblings in min-tree");
        requireEqual(minTree->getLowestCommonAncestor(4, 5), 2, "LCA cousins in min-tree");

        maxTree->mergeNodeIntoParent(5);
        const NodeId reused = editor.createDetachedNode();
        editor.attach(4, reused);
        editor.reparent(reused, 2);
        requireEqual(maxTree->getLowestCommonAncestor(4, reused), 2, "LCA after moveNode");
        requireEqual(maxTree->getLowestCommonAncestor(reused, 3), 2, "LCA after reattachment");
    }

    {
        constexpr NodeId branchLength = 2048;
        constexpr NodeId numProperParts = 1;
        constexpr NodeId firstBranchLeaf = branchLength;
        constexpr NodeId secondBranchRoot = branchLength + 1;
        constexpr NodeId secondBranchLeaf = branchLength * 2;
        constexpr NodeId numNodeSlots = secondBranchLeaf + 1;

        auto higraNodeId = [](NodeId slotId) {
            return numProperParts + slotId;
        };

        std::vector<NodeId> parent(static_cast<size_t>(numProperParts + numNodeSlots), InvalidNode);
        parent[0] = higraNodeId(firstBranchLeaf);
        parent[static_cast<size_t>(higraNodeId(0))] = higraNodeId(0);
        for (NodeId nodeId = 1; nodeId <= firstBranchLeaf; ++nodeId) {
            parent[static_cast<size_t>(higraNodeId(nodeId))] = higraNodeId(nodeId - 1);
        }
        parent[static_cast<size_t>(higraNodeId(secondBranchRoot))] = higraNodeId(0);
        for (NodeId nodeId = secondBranchRoot + 1; nodeId <= secondBranchLeaf; ++nodeId) {
            parent[static_cast<size_t>(higraNodeId(nodeId))] = higraNodeId(nodeId - 1);
        }

        auto deepTree = MorphologicalTree::createFromHigraParent(
            parent,
            1,
            1,
            MorphologicalTree::MAX_TREE,
            AdjacencyRelation(1, 1, 1.5));
        requireEqual(deepTree.getNodeNumDescendants(0), numNodeSlots - 1, "deep tree descendant count");
        requireEqual(deepTree.getLowestCommonAncestor(firstBranchLeaf, secondBranchLeaf), 0, "deep tree LCA across branches");
    }

    {
        requireThrows(
            [] {
                std::vector<NodeId> noAdjacencyParent = {1, 1};
                auto tree = MorphologicalTree::createFromHigraParent(
                    noAdjacencyParent,
                    1,
                    1,
                    MorphologicalTree::MAX_TREE);
                static_cast<void>(tree);
            },
            "Higra max-tree import must reject missing adjacency"
        );

        {
            std::vector<NodeId> tosParent = {1, 1};
            auto tree = MorphologicalTree::createFromHigraParent(
                tosParent,
                1,
                1,
                MorphologicalTree::TREE_OF_SHAPES);
            require(!tree.hasAdjacencyRelation(), "Higra tree-of-shapes import may omit adjacency");
        }

        requireThrows(
            [] {
                std::vector<NodeId> badOwner = {0, 2, 2};
                auto tree = MorphologicalTree::createFromHigraParent(
                    badOwner,
                    1,
                    2,
                    MorphologicalTree::MAX_TREE,
                    AdjacencyRelation(1, 2, 1.5));
                static_cast<void>(tree);
            },
            "Higra import must reject invalid proper-part owners"
        );

        requireThrows(
            [] {
                std::vector<NodeId> cyclicParent = {1, 2, 1};
                auto tree = MorphologicalTree::createFromHigraParent(
                    cyclicParent,
                    1,
                    1,
                    MorphologicalTree::MAX_TREE,
                    AdjacencyRelation(1, 1, 1.5));
                static_cast<void>(tree);
            },
            "Higra import must reject disconnected cycles"
        );

        requireThrows(
            [] {
                std::vector<NodeId> noInternalNodes = {0, 0, 0, 0};
                auto tree = MorphologicalTree::createFromHigraParent(
                    noInternalNodes,
                    2,
                    2,
                    MorphologicalTree::MAX_TREE,
                    AdjacencyRelation(2, 2, 1.5));
                static_cast<void>(tree);
            },
            "Higra import must reject missing internal nodes"
        );

        requireThrows(
            [] {
                std::vector<NodeId> badHigraParent = {1, 2, 3, 3};
                std::vector<AltitudeType> altitude = {0, 0, 0, 0};
                auto weighted = WeightedMorphologicalTree::createFromHigraParent(
                    badHigraParent,
                    altitude,
                    1,
                    2,
                    MorphologicalTree::MAX_TREE,
                    AdjacencyRelation(1, 2, 1.5));
                static_cast<void>(weighted);
            },
            "Higra import must reject leaves pointing outside the internal-node domain"
        );

        requireThrows(
            [] {
                auto weighted = makeWeightedComponentTree(makeComponentTreeFixture(), true);
                auto [higraParent, higraAltitude] = weighted->exportHigraHierarchy();
                higraAltitude.pop_back();
                auto imported = WeightedMorphologicalTree::createFromHigraParent(
                    higraParent,
                    higraAltitude,
                    4,
                    4,
                    MorphologicalTree::MAX_TREE,
                    AdjacencyRelation(4, 4, 1.5));
                static_cast<void>(imported);
            },
            "Higra import must reject altitude buffers with the wrong size"
        );

        requireThrows(
            [] {
                auto weighted = WeightedMorphologicalTree::createComponentTree(makeComponentTreeFixture(), true);
                static_cast<void>(weighted.getAltitude(999));
            },
            "weighted altitude access must reject invalid node ids"
        );
    }

    return 0;
}
