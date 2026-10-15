#include "support/TestSupport.hpp"

#include "mmcfilters/contours/ContourComputation.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <memory>
#include <ranges>
#include <span>
#include <string>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

namespace {

std::vector<PixelId> pixelsOfConnectedComponent(const MorphologicalTree& tree, NodeId nodeId) {
    std::vector<PixelId> pixels;
    for (NodeId subtreeNodeId : tree.subtreeNodes(nodeId)) {
        for (PixelId pixel : tree.properPart(subtreeNodeId)) {
            pixels.push_back(pixel);
        }
    }
    return pixels;
}

std::vector<PixelId> contourVector(const ContourComputation& contours, NodeId nodeId) {
    std::vector<PixelId> values;
    for (PixelId pixel : contours.contour(nodeId)) {
        values.push_back(pixel);
    }
    return values;
}

std::vector<PixelId> expectedContourForNode(const MorphologicalTree& tree, NodeId nodeId) {
    RegularGridAdjacency2D adj4(tree.numRows(), tree.numColumns(), 1.0);
    std::vector<uint8_t> mask(static_cast<std::size_t>(tree.numRows() * tree.numColumns()), 0);
    std::vector<PixelId> ccPixels = pixelsOfConnectedComponent(tree, nodeId);
    for (PixelId pixel : ccPixels) {
        mask[static_cast<std::size_t>(pixel)] = 1;
    }

    std::vector<PixelId> expectedContour;
    for (PixelId pixel : ccPixels) {
        auto [row, column] = ImageUtils::to2D(pixel, tree.numColumns());
        bool isContour = row == 0 || column == 0 || row == tree.numRows() - 1 || column == tree.numColumns() - 1;
        if (!isContour) {
            for (PixelId neighbor : adj4.getAdjacentIndices(pixel)) {
                if (!mask[static_cast<std::size_t>(neighbor)]) {
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

static_assert(std::input_iterator<ContourComputation::iterator>);
static_assert(std::ranges::input_range<const ContourComputation>);
static_assert(!std::ranges::forward_range<ContourComputation>);

void verifyContoursAgainstSupportMasks(const MorphologicalTree& tree, const std::string& label) {
    ContourComputation contours(tree);
    std::vector<bool> seen(tree.numInternalNodeSlots(), false);
    std::size_t count = 0;
    auto check = [&](NodeId node, std::span<const PixelId> pixels) {
        require(!seen[node], label + " visits each node once");
        for (NodeId child : tree.children(node)) {
            require(seen[child], label + " children precede their parent");
        }
        seen[node] = true;
        ++count;
        auto expected = expectedContourForNode(tree, node);
        std::vector<PixelId> actual(pixels.begin(), pixels.end());
        std::sort(actual.begin(), actual.end());
        requireVectorEqual(actual, expected, label + " streamed contour node " + std::to_string(node));
        auto point = contours.contour(node);
        std::sort(point.begin(), point.end());
        requireVectorEqual(point, expected, label + " point contour node " + std::to_string(node));
    };
    for (auto [node, pixels] : contours) {
        check(node, pixels);
    }
    requireEqual(count, static_cast<std::size_t>(tree.numNodes()), label + " visits all live nodes");
    std::fill(seen.begin(), seen.end(), false);
    count = 0;
    contours.forEachContour(check);
    requireEqual(count, static_cast<std::size_t>(tree.numNodes()), label + " callback visits all live nodes");
}

void verifyIteratorLifetimeAndRestart(const MorphologicalTree& tree) {
    ContourComputation contours(tree);
    auto first = contours.begin();
    auto second = contours.begin();
    auto shared = first;
    const auto firstNode = (*first).first;
    const auto saved = contours.contour(firstNode);
    ++first;
    requireEqual((*shared).first, (*first).first, "iterator copies share their single-pass position");
    requireEqual((*second).first, firstNode, "begin creates independent traversals");
    requireVectorEqual(saved, contours.contour(firstNode), "on-demand query owns its data across advancement");
    auto detached = ContourComputation(tree).begin();
    requireEqual((*detached).first, firstNode, "iterator keeps indexes alive after range destruction");
    std::size_t calls = 0;
    requireThrows<std::runtime_error>([&]() {
        contours.forEachContour([&](NodeId, std::span<const PixelId>) {
            ++calls;
            throw std::runtime_error("stop traversal");
        });
    }, "consumer exceptions propagate");
    requireEqual(calls, std::size_t{1}, "throwing callback stops immediately");
    verifyContoursAgainstSupportMasks(tree, "new traversal after early interruption");
    while (first != contours.end()) {
        ++first;
    }
    require(shared == contours.end(), "copied iterator observes exhaustion");
    requireThrows<std::out_of_range>([&]() { static_cast<void>(*first); }, "exhausted iterator rejects dereference");
}

MorphologicalTree makeTwoBranchInteriorTreeOfShapes() {
    std::vector<NodeId> parent = {
        9, 11, 10, 9, 9, 10, 11, 11, 10, 11, 11, 11,
    };
    std::vector<std::uint8_t> altitude(parent.size(), std::uint8_t{});
    auto valuedTree = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(parent), std::span<const std::uint8_t>(altitude), 3, 3,
                                                                    MorphologicalTreeKind::TreeOfShapes);
    return valuedTree.topology().clone();
}

void verifyDeepChainTraversal() {
    constexpr NodeId chainLeaf = 65539;
    constexpr NodeId stopNode = 5;
    constexpr PixelId centerPixel = 4;
    std::vector<NodeId> parent(static_cast<std::size_t>(chainLeaf + 1), NodeId{0});
    for (NodeId node = 1; node <= chainLeaf; ++node) {
        parent[static_cast<std::size_t>(node)] = node - 1;
    }
    std::array<NodeId, 9> smallestNodes{};
    smallestNodes.fill(stopNode);
    smallestNodes[centerPixel] = chainLeaf;
    const std::vector<uint8_t> altitude(parent.size(), uint8_t{0});
    const auto valuedTree = MorphologicalTreeFactory::createFromNativeTopology(
        std::span<const NodeId>(parent), std::span<const NodeId>(smallestNodes), std::span<const uint8_t>(altitude), NodeId{0}, 3, 3,
        MorphologicalTreeSemantics{MorphologicalTreeKind::Generic, NodeAltitudeOrder::Unconstrained, NoConstructionContext{}});
    auto contours = ContourComputation(valuedTree.topology());

    const std::vector<PixelId> exteriorPixels{0, 1, 2, 3, 5, 6, 7, 8};
    requireVectorEqual(contours.contour(chainLeaf), {centerPixel}, "on-demand query on deep leaf");
    std::size_t count = 0;
    for (auto [node, pixels] : contours) {
        std::vector<PixelId> actual(pixels.begin(), pixels.end());
        std::sort(actual.begin(), actual.end());
        requireVectorEqual(actual, node > stopNode ? std::vector<PixelId>{centerPixel} : exteriorPixels,
                           "deep chain applies the stop event without recursion or generation wrap");
        ++count;
    }
    requireEqual(count, static_cast<std::size_t>(chainLeaf + 1), "deep chain visits every node");
    auto root = contours.contour(0);
    std::sort(root.begin(), root.end());
    requireVectorEqual(root, exteriorPixels, "root query after full traversal");
}

} // namespace

int main() {
    verifyDeepChainTraversal();
    auto image = makeComponentTreeFixture();

    for (bool isMaxtree : {true, false}) {
        auto tree = makeComponentTree(image, isMaxtree);
        verifyContoursAgainstSupportMasks(*tree, isMaxtree ? "max-tree" : "min-tree");

        if (isMaxtree && contract::validationsEnabled) {
            auto staleContours = ContourComputation(*tree);
            tree->mergeNodeIntoParent(4);
            requireThrows<std::logic_error>([&]() { static_cast<void>(staleContours.contour(tree->root())); },
                                            "contours reject point access after topology mutation");
            requireThrows<std::logic_error>([&]() { static_cast<void>(staleContours.begin()); },
                                            "contours reject iteration after topology mutation");
        }

        auto valuedTree = makeValuedComponentTree(image, isMaxtree);
        std::vector<std::int16_t> int16Altitude;
        int16Altitude.reserve(valuedTree->nodeAltitudes().size());
        for (std::uint8_t level : valuedTree->nodeAltitudes()) {
            int16Altitude.push_back(static_cast<std::int16_t>(level));
        }
        const ValuedMorphologicalTreeView<std::int16_t> int16View(valuedTree->topology(), std::span<const std::int16_t>(int16Altitude.data(), int16Altitude.size()));
        auto topologyContours = ContourComputation(valuedTree->topology());
        auto viewContours = ContourComputation(valuedTree->asView());
        auto int16ViewContours = ContourComputation(int16View);
        if (isMaxtree && contract::validationsEnabled) {
            auto staleValuedTree = makeValuedComponentTree(image, true);
            const auto staleView = staleValuedTree->asView();
            staleValuedTree->mergeNodeIntoParent(4);
            requireThrows<std::logic_error>([&]() { static_cast<void>(ContourComputation(staleView)); },
                                            "contour extraction must reject stale ValuedMorphologicalTreeView");
        }
        for (NodeId nodeId : valuedTree->topology().aliveNodeIds()) {
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

    // Branching, diagonal contacts, holes and border pixels, independently
    // checked against masks for both component-tree construction adjacencies.
    for (int pattern = 0; pattern < 512; ++pattern) {
        auto binary = ImageUInt8::create(3, 3);
        for (PixelId pixel = 0; pixel < 9; ++pixel) {
            (*binary)[pixel] = static_cast<uint8_t>((pattern >> pixel) & 1);
        }
        for (double radius : {1.0, 1.5}) {
            auto maxTree = MorphologicalTreeFactory::createMaxTree(binary, radius);
            auto minTree = MorphologicalTreeFactory::createMinTree(binary, radius);
            verifyContoursAgainstSupportMasks(maxTree.topology(), "binary max-tree");
            verifyContoursAgainstSupportMasks(minTree.topology(), "binary min-tree");
        }
    }

    {
        auto tree = makeTreeOfShapes(image, TestTopographicImmersion::SelfDualSpan);
        verifyContoursAgainstSupportMasks(*tree, "tree of shapes");
    }

    {
        auto tree = makeTwoBranchInteriorTreeOfShapes();
        verifyContoursAgainstSupportMasks(tree, "two-branch interior tree of shapes");
    }

    {
        auto tree = makeComponentTree(image, true);
        tree->mergeNodeIntoParent(4);

        auto contours = ContourComputation(*tree);
        requireVectorEqual(collectNodeIds(tree->aliveNodeIds()), {0, 1, 2, 3, 5}, "alive node ids after contour middle-slot merge");
        verifyContoursAgainstSupportMasks(*tree, "merged tree with dead slot");
        verifyIteratorLifetimeAndRestart(*tree);
        requireThrows<std::out_of_range>([&]() { static_cast<void>(contours.contour(4)); }, "on-demand query rejects dead slots");
        requireThrows<std::out_of_range>([&]() { static_cast<void>(contours.contour(InvalidNode)); }, "on-demand query rejects invalid nodes");
    }

    if (contract::validationsEnabled) {
        auto tree = makeComponentTree(image, true);
        ContourComputation contours(*tree);
        auto it = contours.begin();
        tree->mergeNodeIntoParent(4);
        requireThrows<std::logic_error>([&]() { ++it; }, "iterator rejects mutation on advance");
        requireThrows<std::logic_error>([&]() { static_cast<void>(*it); }, "iterator rejects mutation on dereference");

        auto callbackTree = makeComponentTree(image, true);
        requireThrows<std::logic_error>([&]() {
            ContourComputation(*callbackTree).forEachContour([&](NodeId node, std::span<const PixelId>) {
                if (node == callbackTree->root()) {
                    callbackTree->mergeNodeIntoParent(4);
                }
            });
        }, "callback mutation is detected even on the last contour");
    }

    return 0;
}
