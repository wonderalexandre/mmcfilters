#include "support/TestSupport.hpp"

#include <memory>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

namespace {

std::vector<int> pixelsOfConnectedComponent(const MorphologicalTree& tree, NodeId nodeId) {
    std::vector<int> pixels;
    for (NodeId subtreeNodeId : tree.getNodeSubtree(nodeId)) {
        for (int properPart : tree.getProperParts(subtreeNodeId)) {
            pixels.push_back(properPart);
        }
    }
    return pixels;
}

} // namespace

int main() {
    auto image = makeComponentTreeFixture();
    auto tree = makeComponentTree(image, true);

    requireEqual(tree->getNumProperParts(0), 1, "root direct proper parts");
    requireEqual(tree->getNumProperParts(1), 3, "node 1 direct proper parts");
    requireEqual(tree->getNumProperParts(2), 4, "node 2 direct proper parts");
    requireEqual(tree->getNumProperParts(3), 3, "node 3 direct proper parts");
    requireEqual(tree->getNumProperParts(4), 3, "node 4 direct proper parts");
    requireEqual(tree->getNumProperParts(5), 2, "node 5 direct proper parts");

    requireVectorEqual(collectNodeIds(tree->getProperParts(0)), {15}, "root proper parts");
    requireVectorEqual(collectNodeIds(tree->getProperParts(1)), {8, 12, 13}, "node 1 proper parts");
    requireVectorEqual(collectNodeIds(tree->getProperParts(2)), {2, 3, 7, 11}, "node 2 proper parts");
    requireVectorEqual(collectNodeIds(tree->getProperParts(3)), {0, 1, 4}, "node 3 proper parts");
    requireVectorEqual(collectNodeIds(tree->getProperParts(4)), {5, 6, 9}, "node 4 proper parts");
    requireVectorEqual(collectNodeIds(tree->getProperParts(5)), {10, 14}, "node 5 proper parts");

    requireEqual(tree->getProperPartOwner(15), 0, "pixel 15 proper-part owner");
    requireEqual(tree->getProperPartOwner(8), 1, "pixel 8 proper-part owner");
    requireEqual(tree->getProperPartOwner(2), 2, "pixel 2 proper-part owner");
    requireEqual(tree->getProperPartOwner(0), 3, "pixel 0 proper-part owner");
    requireEqual(tree->getProperPartOwner(5), 4, "pixel 5 proper-part owner");
    requireEqual(tree->getProperPartOwner(10), 5, "pixel 10 proper-part owner");

    requireVectorEqual(collectNodeIds(tree->getProperParts(3)), {0, 1, 4}, "node 3 direct proper parts");
    requireVectorEqual(collectNodeIds(tree->getProperParts(4)), {5, 6, 9}, "node 4 direct proper parts");
    requireVectorEqual(pixelsOfConnectedComponent(*tree, 3), {0, 1, 4, 5, 6, 9, 10, 14}, "node 3 pixels of CC");
    requireVectorEqual(pixelsOfConnectedComponent(*tree, 4), {5, 6, 9, 10, 14}, "node 4 pixels of CC");

    return 0;
}
