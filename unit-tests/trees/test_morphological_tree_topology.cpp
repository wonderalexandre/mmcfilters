#include "support/TestSupport.hpp"
#include "mmcfilters/attributes/AttributeComputation.hpp"

#include <algorithm>
#include <memory>
#include <numeric>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

int main() {
    auto computeArea = [](const MorphologicalTree& tree, NodeId nodeId) {
        auto [attrNames, buffer] = AttributeComputation::computeSingleTopologyAttribute(tree, Area);
        return static_cast<int>(buffer[attrNames.linearIndex(nodeId, Area)]);
    };

    auto image = makeComponentTreeFixture();
    auto maxTree = makeComponentTree(image, true);
    auto minTree = makeComponentTree(image, false);
    auto valuedMaxTree = makeValuedComponentTree(image, true);
    auto valuedMinTree = makeValuedComponentTree(image, false);
    const auto maxHigraParent = exportFlatHigraHierarchy(*maxTree).first;

    auto requireDfsIntervals = [](const MorphologicalTree& tree, const std::string& label) {
        const auto preOrderNodes = collectNodeIds(tree.subtreeNodes(tree.root()));
        const auto postOrderNodes = collectNodeIds(tree.postOrder());
        requireEqual(static_cast<int>(preOrderNodes.size()), tree.numNodes(), label + " pre-order coverage");
        requireEqual(static_cast<int>(postOrderNodes.size()), tree.numNodes(), label + " post-order coverage");

        auto entryOrderedNodes = collectNodeIds(tree.aliveNodeIds());
        auto exitOrderedNodes = entryOrderedNodes;
        std::sort(entryOrderedNodes.begin(), entryOrderedNodes.end(), [&](NodeId lhs, NodeId rhs) { return tree.dfsEntryIndex(lhs) < tree.dfsEntryIndex(rhs); });
        std::sort(exitOrderedNodes.begin(), exitOrderedNodes.end(), [&](NodeId lhs, NodeId rhs) { return tree.dfsExitIndex(lhs) < tree.dfsExitIndex(rhs); });
        requireVectorEqual(entryOrderedNodes, preOrderNodes, label + " DFS entry order equals pre-order");
        requireVectorEqual(exitOrderedNodes, postOrderNodes, label + " DFS exit order equals post-order");

        std::vector<int> eventIndices;
        eventIndices.reserve(static_cast<std::size_t>(2 * tree.numNodes()));
        for (NodeId nodeId : tree.aliveNodeIds()) {
            const int entryIndex = tree.dfsEntryIndex(nodeId);
            const int exitIndex = tree.dfsExitIndex(nodeId);
            require(entryIndex < exitIndex, label + " DFS entry precedes exit");
            requireEqual(tree.numDescendants(nodeId), (exitIndex - entryIndex - 1) / 2, label + " descendant count from DFS interval");
            eventIndices.push_back(entryIndex);
            eventIndices.push_back(exitIndex);
        }
        std::sort(eventIndices.begin(), eventIndices.end());
        std::vector<int> expected(static_cast<std::size_t>(2 * tree.numNodes()));
        std::iota(expected.begin(), expected.end(), 0);
        requireVectorEqual(eventIndices, expected, label + " DFS events form one interleaved permutation");
        requireEqual(tree.dfsEntryIndex(tree.root()), 0, label + " root DFS entry index");
        requireEqual(tree.dfsExitIndex(tree.root()), 2 * tree.numNodes() - 1, label + " root DFS exit index");

        for (NodeId ancestor : tree.aliveNodeIds()) {
            for (NodeId node : tree.aliveNodeIds()) {
                const bool intervalContains = tree.dfsEntryIndex(ancestor) <= tree.dfsEntryIndex(node) &&
                                              tree.dfsExitIndex(ancestor) >= tree.dfsExitIndex(node);
                requireEqual(tree.isAncestor(ancestor, node), intervalContains, label + " ancestry equals DFS interval containment");
                requireEqual(tree.isDescendant(node, ancestor), intervalContains, label + " descendancy equals reversed DFS containment");
            }
        }
    };

    requireEqual(maxTree->numPixels(), 16, "max-tree num proper parts");
    requireEqual(maxTree->numInternalNodeSlots(), 6, "max-tree internal slots");
    requireEqual(maxTree->numInternalNodeSlots(), 6, "max-tree node-id space size");
    requireEqual(maxTree->root(), 0, "max-tree root alias");
    requireEqual(maxTree->root(), 0, "max-tree root dense node id");
    requireEqual(maxTree->getNumFreeNodeSlots(), 0, "max-tree free slots");
    requireEqual(maxTree->numLeafNodes(), 1, "max-tree leaf nodes");

    requireVectorEqual(collectNodeIds(maxTree->aliveNodeIds()), {0, 1, 2, 3, 4, 5}, "max-tree alive node ids");
    requireVectorEqual(collectNodeIds(maxTree->children(0)), {1}, "max-tree root children");
    requireVectorEqual(collectNodeIds(maxTree->children(3)), {4}, "max-tree node 3 children");
    requireVectorEqual(collectNodeIds(maxTree->postOrder()), {5, 4, 3, 2, 1, 0}, "max-tree post-order");
    requireVectorEqual(collectNodeIds(maxTree->breadthFirstTraversal()), {0, 1, 2, 3, 4, 5}, "max-tree breadth-first");
    requireVectorEqual(collectNodeIds(maxTree->ancestors(5)), {5, 4, 3, 2, 1, 0}, "max-tree path to root");
    requireVectorEqual(collectNodeIds(maxTree->subtreeNodes(2)), {2, 3, 4, 5}, "max-tree subtree");
    requireVectorEqual(collectNodeIds(maxTree->descendants(2)), {3, 4, 5}, "max-tree descendants");
    requireEqual(computeArea(*maxTree, 3), 8, "max-tree node area by node id");
    requireEqual(maxTree->numDescendants(2), 3, "max-tree descendants count by node id");
    requireEqual(maxTree->numSiblings(4), 0, "max-tree siblings count by node id");
    requireEqual(maxTree->properPartCardinality(3), 3, "max-tree direct proper-part count by node id");
    requireDfsIntervals(*maxTree, "max-tree");
    requireEqual(valuedMaxTree->nodeAltitude(0), 0, "valuedTree max-tree root altitude by node id");
    requireEqual(valuedMaxTree->nodeResidue(5), 1, "valuedTree max-tree node residue by node id");

    requireEqual(maxTree->parent(5), 4, "max-tree node parent");
    requireEqual(maxTree->parent(0), 0, "max-tree root parent must point to itself");
    require(maxTree->isNode(5), "5 must be a node in max-tree");
    require(!maxTree->isNode(10), "10 must not be an internal node in max-tree");
    require(maxTree->isPixel(10), "10 must be a proper part in max-tree");
    require(maxTree->isAlive(5), "5 must be alive in max-tree");
    require(maxTree->isRoot(0), "0 must be root in max-tree");
    require(!maxTree->isLeaf(4), "4 must not be leaf in max-tree");
    require(maxTree->isLeaf(5), "5 must be leaf in max-tree");
    require(!maxTree->isAlive(InvalidNode), "invalid node id must not be alive");
    if constexpr (contract::validationsEnabled) {
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(maxTree->parent(InvalidNode)); }, "invalid parent must throw");
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(collectNodeIds(maxTree->children(InvalidNode))); }, "invalid children must throw");
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(maxTree->numChildren(999)); }, "invalid numChildren must throw");
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(maxTree->dfsEntryIndex(999)); }, "invalid dfsEntryIndex must throw");
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(collectPixelIds(maxTree->properPart(999))); }, "invalid properPart must throw");
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(maxTree->isAncestor(999, 0)); }, "invalid isAncestor source must throw");
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(maxTree->isAncestor(0, 999)); }, "invalid isAncestor target must throw");
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(maxTree->isDescendant(999, 0)); }, "invalid isDescendant source must throw");
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(maxTree->isDescendant(0, 999)); }, "invalid isDescendant target must throw");
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(valuedMaxTree->nodeAltitude(InvalidNode)); }, "invalid valuedTree nodeAltitude must throw");
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(valuedMaxTree->nodeResidue(999)); }, "invalid valuedTree nodeResidue must throw");
    }

    auto sparseTree = makeComponentTree(image, true);
    sparseTree->mergeNodeIntoParent(4);
    require(!sparseTree->isAlive(4), "merged node slot must no longer be alive");
    if constexpr (contract::validationsEnabled) {
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(sparseTree->parent(4)); }, "dead-slot parent must throw");
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(collectNodeIds(sparseTree->children(4))); }, "dead-slot children must throw");
    }
    requireDfsIntervals(*sparseTree, "sparse max-tree");
    auto sparseValuedTree = makeValuedComponentTree(image, true);
    sparseValuedTree->mergeNodeIntoParent(4);
    if constexpr (contract::validationsEnabled) {
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(sparseValuedTree->nodeAltitude(4)); }, "dead-slot valuedTree nodeAltitude must throw");
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(sparseValuedTree->nodeResidue(4)); }, "dead-slot valuedTree nodeResidue must throw");
    }
    requireThrows<std::invalid_argument>([&]() { sparseValuedTree->setNodeAltitude(4, 7); }, "dead-slot valuedTree setNodeAltitude must throw");

    auto rebuiltFromHigra = makeTreeFromHigraParent(maxHigraParent, 4, 4, true);
    const NodeId rebuiltRoot = rebuiltFromHigra->root();
    requireVectorEqual(collectNodeIds(rebuiltFromHigra->aliveNodeIds()), {0, 1, 2, 3, 4, 5}, "Higra import alive nodes");
    requireEqual(rebuiltFromHigra->parent(rebuiltRoot), rebuiltRoot, "Higra import root parent must point to itself");
    requireEqual(computeArea(*rebuiltFromHigra, rebuiltRoot), 16, "Higra import root area");
    const auto rebuiltLeaves = rebuiltFromHigra->leaves();
    requireEqual(static_cast<int>(rebuiltLeaves.size()), 1, "Higra import leaf count");
    auto rebuiltLeafPath = collectNodeIds(rebuiltFromHigra->ancestors(rebuiltLeaves.front()));
    requireEqual(static_cast<int>(rebuiltLeafPath.size()), 6, "Higra import path to root length");
    requireEqual(rebuiltLeafPath.back(), rebuiltRoot, "Higra import path must end at root");
    const NodeId rebuiltOwnerOfPixel10 = rebuiltFromHigra->smallestNode(10);
    require(rebuiltFromHigra->isAlive(rebuiltOwnerOfPixel10), "Higra import smallest node must be alive");
    requireDfsIntervals(*rebuiltFromHigra, "rebuilt max-tree");

    requireEqual(minTree->numPixels(), 16, "min-tree num proper parts");
    requireEqual(minTree->numInternalNodeSlots(), 6, "min-tree internal slots");
    requireEqual(minTree->numInternalNodeSlots(), 6, "min-tree node-id space size");
    requireEqual(minTree->root(), 0, "min-tree root node id");
    requireEqual(minTree->numLeafNodes(), 2, "min-tree leaf nodes");

    requireVectorEqual(collectNodeIds(minTree->aliveNodeIds()), {0, 1, 2, 3, 4, 5}, "min-tree alive node ids");
    requireVectorEqual(collectNodeIds(minTree->children(2)), {3, 4}, "min-tree node 2 children");
    requireVectorEqual(collectNodeIds(minTree->ancestors(5)), {5, 3, 2, 1, 0}, "min-tree path to root");
    requireVectorEqual(collectNodeIds(minTree->descendants(2)), {3, 5, 4}, "min-tree descendants");
    requireEqual(minTree->numSiblings(3), 1, "min-tree siblings count by node id");
    requireEqual(minTree->numDescendants(2), 3, "min-tree descendants count by node id");
    requireDfsIntervals(*minTree, "min-tree");
    requireEqual(valuedMinTree->nodeAltitude(0), 5, "valuedTree min-tree root altitude by node id");

    return 0;
}
