#include "support/TestSupport.hpp"
#include "mmcfilters/attributes/AttributeComputedIncrementally.hpp"

#include <memory>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

int main() {
    auto computeArea = [](MorphologicalTree& tree, NodeId nodeId) {
        auto [attrNames, buffer] = AttributeComputedIncrementally::computeSingleAttribute(tree, AREA);
        return static_cast<int>(buffer[attrNames.linearIndex(nodeId, AREA)]);
    };

    auto image = makeComponentTreeFixture();
    auto maxTree = std::make_shared<MorphologicalTree>(image, true);
    auto minTree = std::make_shared<MorphologicalTree>(image, false);
    auto weightedMaxTree = std::make_shared<WeightedMorphologicalTree>(image, true);
    auto weightedMinTree = std::make_shared<WeightedMorphologicalTree>(image, false);
    const auto maxParent = exportParentArray(*maxTree);

    requireEqual(maxTree->getNumTotalProperParts(), 16, "max-tree num proper parts");
    requireEqual(maxTree->getNumInternalNodeSlots(), 6, "max-tree internal slots");
    requireEqual(maxTree->getNumInternalNodeSlots(), 6, "max-tree node-id space size");
    requireEqual(maxTree->getRoot(), 0, "max-tree root alias");
    requireEqual(maxTree->getRoot(), 0, "max-tree root dense node id");
    requireEqual(maxTree->getNumFreeNodeSlots(), 0, "max-tree free slots");
    requireEqual(maxTree->getNumLeafNodes(), 1, "max-tree leaf nodes");

    requireVectorEqual(collectNodeIds(maxTree->getAliveNodeIds()), {0, 1, 2, 3, 4, 5}, "max-tree alive node ids");
    requireVectorEqual(collectNodeIds(maxTree->getChildren(0)), {1}, "max-tree root children");
    requireVectorEqual(collectNodeIds(maxTree->getChildren(3)), {4}, "max-tree node 3 children");
    requireVectorEqual(collectNodeIds(maxTree->getPostOrderNodes()), {5, 4, 3, 2, 1, 0}, "max-tree post-order");
    requireVectorEqual(collectNodeIds(maxTree->getIteratorBreadthFirstTraversal()), {0, 1, 2, 3, 4, 5}, "max-tree breadth-first");
    requireVectorEqual(collectNodeIds(maxTree->getPathToRootNodes(5)), {5, 4, 3, 2, 1, 0}, "max-tree path to root");
    requireVectorEqual(collectNodeIds(maxTree->getNodeSubtree(2)), {2, 3, 4, 5}, "max-tree subtree");
    requireVectorEqual(collectNodeIds(maxTree->getDescendants(2)), {3, 4, 5}, "max-tree descendants");
    requireEqual(computeArea(*maxTree, 3), 8, "max-tree node area by node id");
    requireEqual(maxTree->getNodeNumDescendants(2), 3, "max-tree descendants count by node id");
    requireEqual(maxTree->getNodeNumSiblings(4), 0, "max-tree siblings count by node id");
    requireEqual(maxTree->getNumProperParts(3), 3, "max-tree direct proper-part count by node id");
    require(maxTree->getNodeTimePreOrder(0) < maxTree->getNodeTimePostOrder(0), "max-tree preorder must precede postorder");
    requireEqual(weightedMaxTree->getAltitude(0), 0, "weighted max-tree root altitude by node id");
    requireEqual(weightedMaxTree->getNodeResidue(5), 1, "weighted max-tree node residue by node id");

    requireEqual(maxTree->getNodeParent(5), 4, "max-tree node parent");
    requireEqual(maxTree->getNodeParent(0), 0, "max-tree root parent must point to itself");
    require(maxTree->isNode(5), "5 must be a node in max-tree");
    require(!maxTree->isNode(10), "10 must not be an internal node in max-tree");
    require(maxTree->isProperPart(10), "10 must be a proper part in max-tree");
    require(maxTree->isAlive(5), "5 must be alive in max-tree");
    require(maxTree->isRoot(0), "0 must be root in max-tree");
    require(!maxTree->isLeaf(4), "4 must not be leaf in max-tree");
    require(maxTree->isLeaf(5), "5 must be leaf in max-tree");

    auto rebuiltFromParent = std::make_shared<MorphologicalTree>(std::span<const NodeId>(maxParent), 4, 4, true);
    std::cerr << "after parent-array constructor\n";
    requireEqual(rebuiltFromParent->getRoot(), 0, "parent-array constructor root alias");
    requireVectorEqual(collectNodeIds(rebuiltFromParent->getAliveNodeIds()), {0, 1, 2, 3, 4, 5}, "parent-array constructor alive nodes");
    requireVectorEqual(collectNodeIds(rebuiltFromParent->getChildren(3)), {4}, "parent-array constructor children");
    requireVectorEqual(collectNodeIds(rebuiltFromParent->getPathToRootNodes(5)), {5, 4, 3, 2, 1, 0}, "parent-array constructor path to root");
    requireEqual(computeArea(*rebuiltFromParent, 3), 8, "parent-array constructor area");
    requireEqual(rebuiltFromParent->getSmallestComponent(10), 5, "parent-array constructor smallest component");
    requireEqual(rebuiltFromParent->getNodeParent(0), 0, "parent-array constructor root parent must point to itself");

    auto resetFromParent = MorphologicalTree::create(4, 4, true, AdjacencyRelation(4, 4, 1.5));
    resetFromParent.reset(maxParent);
    requireEqual(resetFromParent.getRoot(), 0, "parent-array reset root alias");
    requireVectorEqual(collectNodeIds(resetFromParent.getAliveNodeIds()), {0, 1, 2, 3, 4, 5}, "parent-array reset alive nodes");
    requireVectorEqual(collectNodeIds(resetFromParent.getChildren(3)), {4}, "parent-array reset children");
    requireVectorEqual(collectNodeIds(resetFromParent.getPathToRootNodes(5)), {5, 4, 3, 2, 1, 0}, "parent-array reset path to root");
    requireEqual(computeArea(resetFromParent, 3), 8, "parent-array reset area");
    requireEqual(resetFromParent.getSmallestComponent(10), 5, "parent-array reset smallest component");
    requireEqual(resetFromParent.getNodeParent(0), 0, "parent-array reset root parent must point to itself");

    requireEqual(minTree->getNumTotalProperParts(), 16, "min-tree num proper parts");
    requireEqual(minTree->getNumInternalNodeSlots(), 6, "min-tree internal slots");
    requireEqual(minTree->getNumInternalNodeSlots(), 6, "min-tree node-id space size");
    requireEqual(minTree->getRoot(), 0, "min-tree root node id");
    requireEqual(minTree->getNumLeafNodes(), 2, "min-tree leaf nodes");

    requireVectorEqual(collectNodeIds(minTree->getAliveNodeIds()), {0, 1, 2, 3, 4, 5}, "min-tree alive node ids");
    requireVectorEqual(collectNodeIds(minTree->getChildren(2)), {3, 4}, "min-tree node 2 children");
    requireVectorEqual(collectNodeIds(minTree->getPathToRootNodes(5)), {5, 3, 2, 1, 0}, "min-tree path to root");
    requireVectorEqual(collectNodeIds(minTree->getDescendants(2)), {3, 5, 4}, "min-tree descendants");
    requireEqual(minTree->getNodeNumSiblings(3), 1, "min-tree siblings count by node id");
    requireEqual(minTree->getNodeNumDescendants(2), 3, "min-tree descendants count by node id");
    requireEqual(weightedMinTree->getAltitude(0), 5, "weighted min-tree root altitude by node id");

    return 0;
}
