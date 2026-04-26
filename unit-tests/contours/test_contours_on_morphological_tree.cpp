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

std::vector<int> contourVector(
    const ContoursComputedIncrementally::IncrementalContours& contours,
    NodeId nodeId) {
    std::vector<int> values;
    for (int pixel : contours.getContour(nodeId)) {
        values.push_back(pixel);
    }
    return values;
}

}

int main() {
    auto image = makeComponentTreeFixture();

    for (bool isMaxtree : {true, false}) {
        auto tree = makeComponentTree(image, isMaxtree);
        auto contours = ContoursComputedIncrementally::extractCompactContours(*tree);
        require(!contours.isMaterialized(), "contours must start without global materialization");
        AdjacencyRelation adj4(
            tree->getNumRowsOfImage(),
            tree->getNumColsOfImage(),
            ContoursComputedIncrementally::ContourSideAdjacencyRadius);
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

            std::vector<int> actualContour = contourVector(contours, nodeId);
            std::vector<int> rangeContour = contourVector(contours, nodeId);
            std::sort(expectedContour.begin(), expectedContour.end());
            std::sort(actualContour.begin(), actualContour.end());
            std::sort(rangeContour.begin(), rangeContour.end());
            requireVectorEqual(actualContour, expectedContour, isMaxtree ? "max-tree contour regression" : "min-tree contour regression");
            requireVectorEqual(rangeContour, actualContour, isMaxtree ? "max-tree range contour regression" : "min-tree range contour regression");

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
        require(!contours.isMaterialized(), "contours must not materialize during extraction");
        std::vector<int> leafContourFromIncrementalRead = contourVector(contours, 5);
        require(contours.isContourMaterialized(5), "contour iteration must cache the requested leaf");
        require(!contours.isContourMaterialized(tree->getRoot()), "leaf contour iteration must not materialize the root");
        require(!contours.isMaterialized(), "contour iteration must not materialize all contours");

        auto contoursByNode = contours.contoursByNode();
        require(!contours.isMaterialized(), "creating the contours-by-node range must not materialize contours");

        std::vector<NodeId> contourNodes;
        for (auto&& [nodeId, contour] : contoursByNode) {
            contourNodes.push_back(nodeId);
            std::vector<int> values;
            for (int pixel : contour) {
                values.push_back(pixel);
            }
            std::vector<int> eagerValues = contourVector(contours, nodeId);
            std::sort(values.begin(), values.end());
            std::sort(eagerValues.begin(), eagerValues.end());
            requireVectorEqual(values, eagerValues, "contours-by-node range must match eager contour after middle-slot merge");
        }

        requireVectorEqual(contourNodes, {0, 1, 2, 3, 5}, "contours-by-node node ids after middle-slot merge");
        require(contours.isMaterialized(), "iterating all contour ranges must materialize all contours incrementally");

        contours.materializeAll();
        require(contours.isMaterialized(), "explicit global materialization must materialize all contours");
        require(!contourVector(contours, 5).empty(), "leaf contour must stay readable after middle-slot merge");
        auto leafMaterialized = contourVector(contours, 5);
        std::sort(leafContourFromIncrementalRead.begin(), leafContourFromIncrementalRead.end());
        std::sort(leafMaterialized.begin(), leafMaterialized.end());
        requireVectorEqual(leafContourFromIncrementalRead, leafMaterialized, "incremental contour read must match materialized contour");
        requireThrows<std::invalid_argument>(
            [&]() { static_cast<void>(contourVector(contours, 4)); },
            "contour access must reject dead slots");
        requireThrows<std::invalid_argument>(
            [&]() { static_cast<void>(contours.getContour(InvalidNode)); },
            "contour access must reject invalid node ids");
    }

    return 0;
}
