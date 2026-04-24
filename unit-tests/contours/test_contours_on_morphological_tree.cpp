#include "support/TestSupport.hpp"

#include "mmcfilters/contours/ContoursComputedIncrementally.hpp"

#include <algorithm>
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

}

int main() {
    auto image = makeComponentTreeFixture();

    for (bool isMaxtree : {true, false}) {
        auto tree = makeComponentTree(image, isMaxtree);
        auto contours = ContoursComputedIncrementally::extractCompactContours(*tree);
        AdjacencyRelation adj4(tree->getNumRowsOfImage(), tree->getNumColsOfImage(), 1);
        std::vector<uint8_t> mask(static_cast<std::size_t>(image->getSize()), 0);

        for (NodeId nodeId : tree->getAliveNodeIds()) {
            std::vector<int> ccPixels;
            for (int pixel : pixelsOfConnectedComponent(*tree, nodeId)) {
                ccPixels.push_back(pixel);
                mask[static_cast<std::size_t>(pixel)] = 1;
            }

            std::vector<int> expectedContour;
            for (int pixel : ccPixels) {
                auto [row, col] = ImageUtils::to2D(pixel, tree->getNumColsOfImage());
                bool isContour = row == 0 || col == 0 || row == tree->getNumRowsOfImage() - 1 || col == tree->getNumColsOfImage() - 1;
                if (!isContour) {
                    for (int q : adj4.getAdjPixels(pixel)) {
                        if (!mask[static_cast<std::size_t>(q)]) {
                            isContour = true;
                            break;
                        }
                    }
                }
                if (isContour) {
                    expectedContour.push_back(pixel);
                }
            }

            std::vector<int> actualContour = contours.buildContourVector(nodeId);
            std::sort(expectedContour.begin(), expectedContour.end());
            std::sort(actualContour.begin(), actualContour.end());
            requireVectorEqual(actualContour, expectedContour, isMaxtree ? "max-tree contour regression" : "min-tree contour regression");

            for (int pixel : ccPixels) {
                mask[static_cast<std::size_t>(pixel)] = 0;
            }
        }
    }

    {
        auto tree = makeComponentTree(image, true);
        tree->mergeNodeIntoParent(4);

        auto contours = ContoursComputedIncrementally::extractCompactContours(*tree);
        requireVectorEqual(collectNodeIds(tree->getAliveNodeIds()), {0, 1, 2, 3, 5}, "alive node ids after contour middle-slot merge");

        std::vector<NodeId> contourNodes;
        for (auto&& [nodeId, contour] : contours.contoursLazy()) {
            contourNodes.push_back(nodeId);
            std::vector<int> values;
            for (int pixel : contour) {
                values.push_back(pixel);
            }
            std::vector<int> eagerValues = contours.buildContourVector(nodeId);
            std::sort(values.begin(), values.end());
            std::sort(eagerValues.begin(), eagerValues.end());
            requireVectorEqual(values, eagerValues, "lazy contour proxy must match eager contour after middle-slot merge");
        }

        requireVectorEqual(contourNodes, {0, 1, 2, 3, 5}, "lazy contour node ids after middle-slot merge");
        require(!contours.buildContourVector(5).empty(), "leaf contour must stay readable after middle-slot merge");
    }

    return 0;
}
