#include "support/TestSupport.hpp"
#include "mmcfilters/attributes/AttributeComputation.hpp"

#include <memory>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

int main() {
    auto computeArea = [](const MorphologicalTree& tree, NodeId nodeId) {
        auto [attrNames, buffer] = AttributeComputation::computeSingleTopologyAttribute(tree, AREA);
        return static_cast<int>(buffer[attrNames.linearIndex(nodeId, AREA)]);
    };

    auto image = makeComponentTreeFixture();
    auto maxTree = makeComponentTree(image, true);
    auto minTree = makeComponentTree(image, false);
    auto weightedMaxTree = makeWeightedComponentTree(image, true);
    auto weightedMinTree = makeWeightedComponentTree(image, false);
    const auto maxHigraParent = exportFlatHigraHierarchy(*maxTree).first;

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
    require(!maxTree->isAlive(InvalidNode), "invalid node id must not be alive");
    requireThrows<std::invalid_argument>([&]() { static_cast<void>(maxTree->getNodeParent(InvalidNode)); }, "invalid getNodeParent must throw");
    requireThrows<std::invalid_argument>([&]() { static_cast<void>(collectNodeIds(maxTree->getChildren(InvalidNode))); }, "invalid getChildren must throw");
    requireThrows<std::invalid_argument>([&]() { static_cast<void>(maxTree->getNumChildren(999)); }, "invalid getNumChildren must throw");
    requireThrows<std::invalid_argument>([&]() { static_cast<void>(maxTree->getNodeTimePreOrder(999)); }, "invalid getNodeTimePreOrder must throw");
    requireThrows<std::invalid_argument>([&]() { static_cast<void>(collectNodeIds(maxTree->getProperParts(999))); }, "invalid getProperParts must throw");
    requireThrows<std::invalid_argument>([&]() { static_cast<void>(maxTree->isAncestor(999, 0)); }, "invalid isAncestor source must throw");
    requireThrows<std::invalid_argument>([&]() { static_cast<void>(maxTree->isAncestor(0, 999)); }, "invalid isAncestor target must throw");
    requireThrows<std::invalid_argument>([&]() { static_cast<void>(maxTree->isDescendant(999, 0)); }, "invalid isDescendant source must throw");
    requireThrows<std::invalid_argument>([&]() { static_cast<void>(maxTree->isDescendant(0, 999)); }, "invalid isDescendant target must throw");
    requireThrows<std::invalid_argument>([&]() { static_cast<void>(weightedMaxTree->getAltitude(InvalidNode)); }, "invalid weighted getAltitude must throw");
    requireThrows<std::invalid_argument>([&]() { static_cast<void>(weightedMaxTree->getNodeResidue(999)); }, "invalid weighted getNodeResidue must throw");

    auto sparseTree = makeComponentTree(image, true);
    sparseTree->mergeNodeIntoParent(4);
    require(!sparseTree->isAlive(4), "merged node slot must no longer be alive");
    requireThrows<std::invalid_argument>([&]() { static_cast<void>(sparseTree->getNodeParent(4)); }, "dead-slot getNodeParent must throw");
    requireThrows<std::invalid_argument>([&]() { static_cast<void>(collectNodeIds(sparseTree->getChildren(4))); }, "dead-slot getChildren must throw");
    auto sparseWeightedTree = makeWeightedComponentTree(image, true);
    sparseWeightedTree->mergeNodeIntoParent(4);
    requireThrows<std::invalid_argument>([&]() { static_cast<void>(sparseWeightedTree->getAltitude(4)); }, "dead-slot weighted getAltitude must throw");
    requireThrows<std::invalid_argument>([&]() { static_cast<void>(sparseWeightedTree->getNodeResidue(4)); }, "dead-slot weighted getNodeResidue must throw");
    requireThrows<std::invalid_argument>([&]() { sparseWeightedTree->setAltitude(4, 7); }, "dead-slot weighted setAltitude must throw");

    auto rebuiltFromHigra = makeTreeFromHigraParent(maxHigraParent, 4, 4, true);
    const NodeId rebuiltRoot = rebuiltFromHigra->getRoot();
    requireVectorEqual(collectNodeIds(rebuiltFromHigra->getAliveNodeIds()), {0, 1, 2, 3, 4, 5}, "Higra import alive nodes");
    requireEqual(rebuiltFromHigra->getNodeParent(rebuiltRoot), rebuiltRoot, "Higra import root parent must point to itself");
    requireEqual(computeArea(*rebuiltFromHigra, rebuiltRoot), 16, "Higra import root area");
    const auto rebuiltLeaves = rebuiltFromHigra->getLeaves();
    requireEqual(static_cast<int>(rebuiltLeaves.size()), 1, "Higra import leaf count");
    auto rebuiltLeafPath = collectNodeIds(rebuiltFromHigra->getPathToRootNodes(rebuiltLeaves.front()));
    requireEqual(static_cast<int>(rebuiltLeafPath.size()), 6, "Higra import path to root length");
    requireEqual(rebuiltLeafPath.back(), rebuiltRoot, "Higra import path must end at root");
    const NodeId rebuiltOwnerOfPixel10 = rebuiltFromHigra->getProperPartOwner(10);
    require(rebuiltFromHigra->isAlive(rebuiltOwnerOfPixel10), "Higra import proper-part owner must be alive");

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
