#include "support/TestSupport.hpp"

#include <memory>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

namespace {

std::vector<int> referenceNodeSupport(const MorphologicalTree& tree, NodeId nodeId) {
    std::vector<PixelId> pixels;
    for (NodeId subtreeNodeId : tree.subtreeNodes(nodeId)) {
        for (PixelId pixel : tree.properPart(subtreeNodeId)) {
            pixels.push_back(pixel);
        }
    }
    return pixels;
}

} // namespace

int main() {
    auto image = makeComponentTreeFixture();
    auto tree = makeComponentTree(image, true);

    require(tree->isTreeOfPartialPartitions(), "component tree must satisfy the tree-of-partial-partitions specialization");
    tree->validateTreeOfPartialPartitions();
    requireEqual(static_cast<int>(tree->smallestNodeMap().size()), tree->numPixels(), "smallest-node map size");

    requireEqual(tree->properPartCardinality(0), 1, "root direct proper parts");
    requireEqual(tree->properPartCardinality(1), 3, "node 1 direct proper parts");
    requireEqual(tree->properPartCardinality(2), 4, "node 2 direct proper parts");
    requireEqual(tree->properPartCardinality(3), 3, "node 3 direct proper parts");
    requireEqual(tree->properPartCardinality(4), 3, "node 4 direct proper parts");
    requireEqual(tree->properPartCardinality(5), 2, "node 5 direct proper parts");

    requireVectorEqual(collectPixelIds(tree->properPart(0)), {15}, "root proper parts");
    requireVectorEqual(collectPixelIds(tree->properPart(1)), {8, 12, 13}, "node 1 proper parts");
    requireVectorEqual(collectPixelIds(tree->properPart(2)), {2, 3, 7, 11}, "node 2 proper parts");
    requireVectorEqual(collectPixelIds(tree->properPart(3)), {0, 1, 4}, "node 3 proper parts");
    requireVectorEqual(collectPixelIds(tree->properPart(4)), {5, 6, 9}, "node 4 proper parts");
    requireVectorEqual(collectPixelIds(tree->properPart(5)), {10, 14}, "node 5 proper parts");

    requireEqual(tree->smallestNode(15), 0, "smallest node containing pixel 15");
    requireEqual(tree->smallestNode(8), 1, "smallest node containing pixel 8");
    requireEqual(tree->smallestNode(2), 2, "smallest node containing pixel 2");
    requireEqual(tree->smallestNode(0), 3, "smallest node containing pixel 0");
    requireEqual(tree->smallestNode(5), 4, "smallest node containing pixel 5");
    requireEqual(tree->smallestNode(10), 5, "smallest node containing pixel 10");

    for (PixelId pixel = 0; pixel < tree->numPixels(); ++pixel) {
        const NodeId smallest = tree->smallestNode(pixel);
        const auto pixels = collectPixelIds(tree->properPart(smallest));
        require(std::find(pixels.begin(), pixels.end(), pixel) != pixels.end(), "pixel must belong to the proper part of its smallest node");
    }

    requireVectorEqual(collectPixelIds(tree->properPart(3)), {0, 1, 4}, "node 3 direct proper parts");
    requireVectorEqual(collectPixelIds(tree->properPart(4)), {5, 6, 9}, "node 4 direct proper parts");
    requireVectorEqual(referenceNodeSupport(*tree, 3), collectNodeIds(tree->nodeSupport(3)), "node support must equal the union of subtree proper parts");
    requireVectorEqual(referenceNodeSupport(*tree, 4), collectNodeIds(tree->nodeSupport(4)), "descendant node support must agree with ancestry");

    return 0;
}
