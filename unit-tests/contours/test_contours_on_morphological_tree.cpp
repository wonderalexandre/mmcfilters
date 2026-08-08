#include "support/TestSupport.hpp"

#include "mmcfilters/contours/ContoursComputedIncrementally.hpp"

#include <algorithm>
#include <cstdint>
#include <memory>
#include <span>
#include <string>

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

std::vector<int> contourVector(const ContoursComputedIncrementally::IncrementalContours& contours, NodeId nodeId) {
    std::vector<int> values;
    for (int pixel : contours.getContour(nodeId)) {
        values.push_back(pixel);
    }
    return values;
}

void appendContourFromDeltas(const MorphologicalTree& tree, const ContoursComputedIncrementally::LocalContourDeltas& deltas, NodeId nodeId,
                             std::vector<int>& values, std::vector<uint8_t>& pixelMark) {
    for (NodeId child : tree.getChildren(nodeId)) {
        appendContourFromDeltas(tree, deltas, child, values, pixelMark);
    }

    for (int pixel : deltas.additions(nodeId)) {
        if (!pixelMark[static_cast<std::size_t>(pixel)]) {
            pixelMark[static_cast<std::size_t>(pixel)] = 1;
            values.push_back(pixel);
        }
    }

    for (int pixel : deltas.removals(nodeId)) {
        pixelMark[static_cast<std::size_t>(pixel)] = 0;
    }

    std::size_t writeIndex = 0;
    for (int pixel : values) {
        if (pixelMark[static_cast<std::size_t>(pixel)]) {
            values[writeIndex++] = pixel;
        }
    }
    values.resize(writeIndex);
}

std::vector<int> contourVectorFromDeltas(const MorphologicalTree& tree, const ContoursComputedIncrementally::LocalContourDeltas& deltas, NodeId nodeId) {
    std::vector<int> values;
    std::vector<uint8_t> pixelMark(static_cast<std::size_t>(tree.getNumRowsOfGridDomain2D() * tree.getNumColsOfGridDomain2D()), 0);
    appendContourFromDeltas(tree, deltas, nodeId, values, pixelMark);
    return values;
}

std::vector<int> expectedContourForNode(const MorphologicalTree& tree, NodeId nodeId) {
    RegularGridAdjacency2D adj4(tree.getNumRowsOfGridDomain2D(), tree.getNumColsOfGridDomain2D(), ContoursComputedIncrementally::ContourSideAdjacencyRadius);
    std::vector<uint8_t> mask(static_cast<std::size_t>(tree.getNumRowsOfGridDomain2D() * tree.getNumColsOfGridDomain2D()), 0);
    std::vector<int> ccPixels = pixelsOfConnectedComponent(tree, nodeId);
    for (int pixel : ccPixels) {
        mask[static_cast<std::size_t>(pixel)] = 1;
    }

    std::vector<int> expectedContour;
    for (int pixel : ccPixels) {
        auto [row, col] = ImageUtils::to2D(pixel, tree.getNumColsOfGridDomain2D());
        bool isContour = row == 0 || col == 0 || row == tree.getNumRowsOfGridDomain2D() - 1 || col == tree.getNumColsOfGridDomain2D() - 1;
        if (!isContour) {
            for (int q : adj4.getAdjacentIndices(pixel)) {
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

    std::sort(expectedContour.begin(), expectedContour.end());
    return expectedContour;
}

void verifyContoursAgainstSupportMasks(const MorphologicalTree& tree, const std::string& label) {
    auto contours = ContoursComputedIncrementally::extractCompactContours(tree);
    auto deltas = ContoursComputedIncrementally::extractContourDeltas(tree);
    require(!contours.isMaterialized(), label + " contours must start without global materialization");

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        std::vector<int> expectedContour = expectedContourForNode(tree, nodeId);
        std::vector<int> actualContour = contourVector(contours, nodeId);
        std::vector<int> rangeContour = contourVector(contours, nodeId);
        std::vector<int> deltaContour = contourVectorFromDeltas(tree, deltas, nodeId);
        std::sort(actualContour.begin(), actualContour.end());
        std::sort(rangeContour.begin(), rangeContour.end());
        std::sort(deltaContour.begin(), deltaContour.end());
        requireVectorEqual(actualContour, expectedContour, label + " contour regression node " + std::to_string(nodeId));
        requireVectorEqual(rangeContour, actualContour, label + " range contour regression node " + std::to_string(nodeId));
        requireVectorEqual(deltaContour, actualContour, label + " contour deltas regression node " + std::to_string(nodeId));
    }
}

MorphologicalTree makeTwoBranchInteriorTreeOfShapes() {
    std::vector<NodeId> parent = {
        9, 11, 10, 9, 9, 10, 11, 11, 10, 11, 11, 11,
    };
    std::vector<std::uint8_t> altitude(parent.size(), std::uint8_t{});
    auto weighted = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(parent), std::span<const std::uint8_t>(altitude), 3, 3,
                                                                    MorphologicalTreeKind::TREE_OF_SHAPES);
    return weighted.topology().clone();
}

} // namespace

int main() {
    auto image = makeComponentTreeFixture();

    for (bool isMaxtree : {true, false}) {
        auto tree = makeComponentTree(image, isMaxtree);
        verifyContoursAgainstSupportMasks(*tree, isMaxtree ? "max-tree" : "min-tree");

        if (isMaxtree) {
            auto staleContours = ContoursComputedIncrementally::extractCompactContours(*tree);
            tree->mergeNodeIntoParent(4);
            requireThrows<std::logic_error>([&]() { static_cast<void>(staleContours.isMaterialized()); },
                                            "contours must reject materialization status after topology mutation");
            requireThrows<std::logic_error>([&]() { static_cast<void>(staleContours.getContour(tree->getRoot())); },
                                            "contours must reject contour access after topology mutation");
            requireThrows<std::logic_error>([&]() { staleContours.materializeAll(); }, "contours must reject materializeAll after topology mutation");
        }

        auto weighted = makeWeightedComponentTree(image, isMaxtree);
        std::vector<std::int16_t> int16Altitude;
        int16Altitude.reserve(weighted->getAltitudeBuffer().size());
        for (std::uint8_t level : weighted->getAltitudeBuffer()) {
            int16Altitude.push_back(static_cast<std::int16_t>(level));
        }
        const WeightedTreeView<std::int16_t> int16View(weighted->topology(), std::span<const std::int16_t>(int16Altitude.data(), int16Altitude.size()));
        auto topologyContours = ContoursComputedIncrementally::extractCompactContours(weighted->topology());
        auto viewContours = ContoursComputedIncrementally::extractCompactContours(weighted->asView());
        auto int16ViewContours = ContoursComputedIncrementally::extractCompactContours(int16View);
        if (isMaxtree) {
            auto staleWeighted = makeWeightedComponentTree(image, true);
            const auto staleView = staleWeighted->asView();
            staleWeighted->mergeNodeIntoParent(4);
            requireThrows<std::logic_error>([&]() { static_cast<void>(ContoursComputedIncrementally::extractCompactContours(staleView)); },
                                            "contour extraction must reject stale WeightedTreeView");
            requireThrows<std::logic_error>([&]() { static_cast<void>(ContoursComputedIncrementally::extractContourDeltas(staleView)); },
                                            "contour delta extraction must reject stale WeightedTreeView");
        }
        for (NodeId nodeId : weighted->topology().getAliveNodeIds()) {
            auto topologyContour = contourVector(topologyContours, nodeId);
            auto viewContour = contourVector(viewContours, nodeId);
            auto int16ViewContour = contourVector(int16ViewContours, nodeId);
            std::sort(topologyContour.begin(), topologyContour.end());
            std::sort(viewContour.begin(), viewContour.end());
            std::sort(int16ViewContour.begin(), int16ViewContour.end());
            requireVectorEqual(viewContour, topologyContour, isMaxtree ? "max-tree contours via view" : "min-tree contours via view");
            requireVectorEqual(int16ViewContour, topologyContour, isMaxtree ? "max-tree contours via int16 view" : "min-tree contours via int16 view");
        }
    }

    {
        auto tree = makeTreeOfShapes(image, ToSInterpolation::SelfDual);
        verifyContoursAgainstSupportMasks(*tree, "tree of shapes");
    }

    {
        auto tree = makeTwoBranchInteriorTreeOfShapes();
        verifyContoursAgainstSupportMasks(tree, "two-branch interior tree of shapes");
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
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(contourVector(contours, 4)); }, "contour access must reject dead slots");
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(contours.getContour(InvalidNode)); }, "contour access must reject invalid node ids");
    }

    return 0;
}
