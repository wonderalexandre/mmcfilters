#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/AttributeComputation.hpp"
#include "mmcfilters/contours/ContourTraceComputation.hpp"
#include "mmcfilters/contours/ContourComputation.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <memory>
#include <span>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

namespace {

static_assert(std::is_same_v<decltype(std::declval<const ContourTraceComputation&>().trace(NodeId{})), ContourTrace>);
static_assert(std::input_iterator<ContourTraceComputation::iterator>);

std::vector<PixelId> pixelsOfConnectedComponent(const MorphologicalTree& tree, NodeId nodeId) {
    std::vector<PixelId> pixels;
    for (NodeId subtreeNodeId : tree.subtreeNodes(nodeId)) {
        for (PixelId pixel : tree.properPart(subtreeNodeId)) {
            pixels.push_back(pixel);
        }
    }
    return pixels;
}

int sideIndex(ContourSide side) { return static_cast<int>(side); }

int packEdge(const ContourEdge& edge) { return ContourTraceComputation::packEdge(edge.pixel, edge.side); }

PixelId neighborPixel(const MorphologicalTree& tree, PixelId pixel, ContourSide side) {
    const int rows = tree.numRows();
    const int columns = tree.numColumns();
    const auto [row, column] = ImageUtils::to2D(pixel, columns);
    switch (side) {
    case ContourSide::North:
        return row == 0 ? InvalidPixel : ImageUtils::to1D(row - 1, column, columns);
    case ContourSide::West:
        return column == 0 ? InvalidPixel : ImageUtils::to1D(row, column - 1, columns);
    case ContourSide::East:
        return column == columns - 1 ? InvalidPixel : ImageUtils::to1D(row, column + 1, columns);
    case ContourSide::South:
        return row == rows - 1 ? InvalidPixel : ImageUtils::to1D(row + 1, column, columns);
    }
    return InvalidPixel;
}

std::vector<int> expectedEdgesForNode(const MorphologicalTree& tree, NodeId nodeId) {
    static constexpr std::array<ContourSide, 4> sides{ContourSide::North, ContourSide::West, ContourSide::East, ContourSide::South};

    std::vector<uint8_t> mask(static_cast<std::size_t>(tree.numPixels()), 0);
    for (PixelId pixel : pixelsOfConnectedComponent(tree, nodeId)) {
        mask[static_cast<std::size_t>(pixel)] = 1;
    }

    std::vector<int> expected;
    for (PixelId pixel = 0; pixel < tree.numPixels(); ++pixel) {
        if (!mask[static_cast<std::size_t>(pixel)]) {
            continue;
        }
        for (ContourSide side : sides) {
            const PixelId neighbor = neighborPixel(tree, pixel, side);
            if (neighbor == InvalidPixel || !mask[static_cast<std::size_t>(neighbor)]) {
                expected.push_back(ContourTraceComputation::packEdge(pixel, side));
            }
        }
    }
    std::sort(expected.begin(), expected.end());
    return expected;
}

template <class Trace> std::vector<int> edgeVector(const Trace& trace) {
    std::vector<int> values;
    for (ContourEdge edge : trace.edges()) {
        values.push_back(packEdge(edge));
    }
    std::sort(values.begin(), values.end());
    return values;
}

std::vector<int> edgeVector(const ContourTraceComputation& traces, NodeId nodeId) {
    return edgeVector(traces.trace(nodeId));
}

std::vector<PixelId> contourPixelProjection(std::span<const int> packedEdges) {
    std::vector<PixelId> pixels;
    pixels.reserve(packedEdges.size());
    for (int packedEdge : packedEdges) {
        pixels.push_back(ContourTraceComputation::unpackEdge(packedEdge).pixel);
    }
    std::sort(pixels.begin(), pixels.end());
    pixels.erase(std::unique(pixels.begin(), pixels.end()), pixels.end());
    return pixels;
}

std::vector<PixelId> contourVector(const ContourComputation& contours, NodeId nodeId) {
    std::vector<PixelId> values;
    for (PixelId pixel : contours.contour(nodeId)) {
        values.push_back(pixel);
    }
    std::sort(values.begin(), values.end());
    return values;
}

template <class Trace> std::vector<int> boundarySignature(const Trace& trace) {
    std::vector<std::array<int, 4>> boundarySummaries;
    for (const ContourBoundary& boundary : trace.boundaries()) {
        int edgeSum = 0;
        for (ContourEdge edge : trace.boundaryEdges(boundary)) {
            edgeSum += packEdge(edge) + 1;
        }
        boundarySummaries.push_back(std::array<int, 4>{static_cast<int>(boundary.kind), static_cast<int>(boundary.edgeCount), boundary.doubledSignedArea, edgeSum});
    }
    std::sort(boundarySummaries.begin(), boundarySummaries.end());

    std::vector<int> signature;
    signature.reserve(boundarySummaries.size() * 4);
    for (const auto& summary : boundarySummaries) {
        signature.insert(signature.end(), summary.begin(), summary.end());
    }
    return signature;
}

std::array<int, 4> directionalCounts(std::span<const int> packedEdges) {
    std::array<int, 4> counts{0, 0, 0, 0};
    for (int packedEdge : packedEdges) {
        ++counts[static_cast<std::size_t>(sideIndex(ContourTraceComputation::unpackEdge(packedEdge).side))];
    }
    return counts;
}

NodeId findNodeByArea(const MorphologicalTree& tree, int area) {
    auto [names, values] = AttributeComputation::computeSingleTopologyAttribute(tree, Area);
    for (NodeId nodeId : tree.aliveNodeIds()) {
        if (static_cast<int>(values[static_cast<std::size_t>(names.linearIndex(nodeId, Area))]) == area) {
            return nodeId;
        }
    }
    return InvalidNode;
}

void verifyTraceEdgesAgainstSupportMasks(const MorphologicalTree& tree, const std::string& label) {
    auto traces = ContourTraceComputation(tree);
    auto contours = ContourComputation(tree);
    std::vector<AttributeOrGroup> requests{ContourPerimeter, ContourSideNorth, ContourSideWest, ContourSideEast, ContourSideSouth};
    auto attributeResult = AttributeComputation::computeTopologyAttributes(tree, requests);

    for (NodeId nodeId : tree.aliveNodeIds()) {
        std::vector<int> actualEdges = edgeVector(traces, nodeId);
        std::vector<int> expectedEdges = expectedEdgesForNode(tree, nodeId);
        requireVectorEqual(actualEdges, expectedEdges, label + " trace edge regression node " + std::to_string(nodeId));

        std::vector<PixelId> projectedPixels = contourPixelProjection(std::span<const int>(actualEdges));
        requireVectorEqual(projectedPixels, contourVector(contours, nodeId), label + " trace projection node " + std::to_string(nodeId));

        const std::array<int, 4> counts = directionalCounts(std::span<const int>(actualEdges));
        const auto& names = attributeResult.attributeNames();
        const auto& values = attributeResult.values();
        requireEqual(static_cast<int>(values[static_cast<std::size_t>(names.linearIndex(nodeId, ContourPerimeter))]), static_cast<int>(actualEdges.size()),
                     label + " perimeter from traced edges node " + std::to_string(nodeId));
        requireEqual(static_cast<int>(values[static_cast<std::size_t>(names.linearIndex(nodeId, ContourSideNorth))]),
                     counts[static_cast<std::size_t>(sideIndex(ContourSide::North))], label + " north side count node " + std::to_string(nodeId));
        requireEqual(static_cast<int>(values[static_cast<std::size_t>(names.linearIndex(nodeId, ContourSideWest))]),
                     counts[static_cast<std::size_t>(sideIndex(ContourSide::West))], label + " west side count node " + std::to_string(nodeId));
        requireEqual(static_cast<int>(values[static_cast<std::size_t>(names.linearIndex(nodeId, ContourSideEast))]),
                     counts[static_cast<std::size_t>(sideIndex(ContourSide::East))], label + " east side count node " + std::to_string(nodeId));
        requireEqual(static_cast<int>(values[static_cast<std::size_t>(names.linearIndex(nodeId, ContourSideSouth))]),
                     counts[static_cast<std::size_t>(sideIndex(ContourSide::South))], label + " south side count node " + std::to_string(nodeId));
    }
}

void verifyOnePixelBoundary() {
    auto image = makeImage(3, 3,
                           {
                               1,
                               1,
                               1,
                               1,
                               2,
                               1,
                               1,
                               1,
                               1,
                           });
    auto tree = makeComponentTree(image, true);
    const NodeId nodeId = findNodeByArea(*tree, 1);
    require(nodeId != InvalidNode, "one-pixel fixture must have an area-1 node");

    auto trace = ContourTraceComputation(*tree).trace(nodeId);
    auto boundaries = trace.boundaries();
    requireEqual(boundaries.size(), std::size_t{1}, "one-pixel support boundary count");
    requireEqual(static_cast<int>(boundaries[0].kind), static_cast<int>(ContourBoundaryKind::External), "one-pixel support boundary kind");
    requireEqual(boundaries[0].edgeCount, uint32_t{4}, "one-pixel support boundary edge count");
    requireEqual(boundaries[0].doubledSignedArea, 2, "one-pixel support signed area");
    requireEqual(trace.boundaryEdges(boundaries[0]).size(), std::size_t{4}, "one-pixel support boundary edge range");
}

void verifyRingBoundarySeparation() {
    auto image = makeImage(3, 3,
                           {
                               2,
                               2,
                               2,
                               2,
                               1,
                               2,
                               2,
                               2,
                               2,
                           });
    auto tree = makeComponentTree(image, true);
    const NodeId ringNode = findNodeByArea(*tree, 8);
    require(ringNode != InvalidNode, "ring fixture must have an area-8 node");

    auto trace = ContourTraceComputation(*tree).trace(ringNode);
    auto boundaries = trace.boundaries();
    requireEqual(boundaries.size(), std::size_t{2}, "ring support boundary count");

    int externalBoundaries = 0;
    int internalBoundaries = 0;
    int externalEdges = 0;
    int internalEdges = 0;
    int totalEdges = 0;
    for (const ContourBoundary& boundary : boundaries) {
        totalEdges += static_cast<int>(boundary.edgeCount);
        if (boundary.kind == ContourBoundaryKind::External) {
            ++externalBoundaries;
            externalEdges += static_cast<int>(boundary.edgeCount);
            require(boundary.doubledSignedArea > 0, "ring external boundary must have positive signed area");
        } else {
            ++internalBoundaries;
            internalEdges += static_cast<int>(boundary.edgeCount);
            require(boundary.doubledSignedArea < 0, "ring internal boundary must have negative signed area");
        }
    }

    requireEqual(externalBoundaries, 1, "ring external boundary count");
    requireEqual(internalBoundaries, 1, "ring internal boundary count");
    requireEqual(externalEdges, 12, "ring external edge count");
    requireEqual(internalEdges, 4, "ring internal edge count");
    requireEqual(totalEdges, 16, "ring total edge count");

    const ContourBoundary external = trace.externalBoundary();
    requireEqual(static_cast<int>(external.kind), static_cast<int>(ContourBoundaryKind::External), "direct external boundary kind");
    requireEqual(external.edgeCount, uint32_t{12}, "direct external boundary edge count");

    std::vector<ContourEdge> orderedEdges;
    for (ContourEdge edge : trace.boundaryEdges(external)) {
        orderedEdges.push_back(edge);
    }
    std::vector<PixelId> orderedPixels;
    for (PixelId pixel : trace.boundaryPixels(external)) {
        orderedPixels.push_back(pixel);
    }
    requireEqual(orderedPixels.size(), orderedEdges.size(), "boundary pixel projection size");
    for (std::size_t index = 0; index < orderedEdges.size(); ++index) {
        requireEqual(orderedPixels[index], orderedEdges[index].pixel, "boundary pixel projection preserves edge order");
    }
    auto distinctPixels = orderedPixels;
    std::sort(distinctPixels.begin(), distinctPixels.end());
    distinctPixels.erase(std::unique(distinctPixels.begin(), distinctPixels.end()), distinctPixels.end());
    require(distinctPixels.size() < orderedPixels.size(), "boundary pixel projection must retain repetitions from distinct pixel sides");
}

void requireBoundarySummary(const MorphologicalTree& tree, NodeId nodeId, int expectedExternalBoundaries, int expectedInternalBoundaries, int expectedExternalEdges,
                        int expectedInternalEdges, const std::string& label) {
    auto trace = ContourTraceComputation(tree).trace(nodeId);
    auto boundaries = trace.boundaries();

    int externalBoundaries = 0;
    int internalBoundaries = 0;
    int externalEdges = 0;
    int internalEdges = 0;
    for (const ContourBoundary& boundary : boundaries) {
        if (boundary.kind == ContourBoundaryKind::External) {
            ++externalBoundaries;
            externalEdges += static_cast<int>(boundary.edgeCount);
            require(boundary.doubledSignedArea > 0, label + " external boundary signed area");
        } else {
            ++internalBoundaries;
            internalEdges += static_cast<int>(boundary.edgeCount);
            require(boundary.doubledSignedArea < 0, label + " internal boundary signed area");
        }
    }

    requireEqual(externalBoundaries, expectedExternalBoundaries, label + " external boundary count");
    requireEqual(internalBoundaries, expectedInternalBoundaries, label + " internal boundary count");
    requireEqual(externalEdges, expectedExternalEdges, label + " external edge count");
    requireEqual(internalEdges, expectedInternalEdges, label + " internal edge count");
    requireEqual(externalEdges + internalEdges, static_cast<int>(trace.edges().size()), label + " boundaries must cover every trace edge");
}

void verifyBorderTouchingBoundary() {
    auto image = makeImage(3, 3,
                           {
                               2,
                               2,
                               1,
                               2,
                               2,
                               1,
                               1,
                               1,
                               1,
                           });
    auto tree = makeComponentTree(image, true);
    const NodeId nodeId = findNodeByArea(*tree, 4);
    require(nodeId != InvalidNode, "border-touching fixture must have an area-4 node");

    requireBoundarySummary(*tree, nodeId, 1, 0, 8, 0, "border-touching support");
}

void verifyMultipleInternalBoundaries() {
    auto image = makeImage(5, 5,
                           {
                               2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2,
                           });
    auto tree = makeComponentTree(image, true);
    const NodeId nodeId = findNodeByArea(*tree, 23);
    require(nodeId != InvalidNode, "two-hole fixture must have an area-23 node");

    requireBoundarySummary(*tree, nodeId, 1, 2, 20, 8, "two-hole support");
}

void verifyDiagonalTouchingBoundariesAreDeterministic() {
    auto image = makeImage(2, 2,
                           {
                               2,
                               1,
                               1,
                               2,
                           });
    auto tree = makeComponentTree(image, true);
    const NodeId nodeId = findNodeByArea(*tree, 2);
    require(nodeId != InvalidNode, "diagonal-touch fixture must have an area-2 node");

    requireBoundarySummary(*tree, nodeId, 1, 0, 8, 0, "eight-connected diagonal-touching support");
}

void verifyClosureAndComplementaryPairing() {
    // The center background pixel reaches the exterior diagonally. Closure must
    // occur at the initial directed edge even when another vertex is revisited.
    auto image = makeImage(3, 3, {1, 0, 2, 0, 2, 0, 0, 0, 0});
    for (double radius : {1.0, 1.5}) {
        auto tree = makeComponentTree(image, false, radius);
        const NodeId node = findNodeByArea(*tree, 7);
        require(node != InvalidNode, "diagonal background fixture must expose the area-7 support");
        if (radius == 1.0) {
            requireBoundarySummary(*tree, node, 1, 0, 16, 0, "4/8 diagonal background");
        } else {
            requireBoundarySummary(*tree, node, 1, 1, 12, 4, "8/4 diagonal background");
        }
    }
}

// Independent topology oracle: count foreground components and bounded
// complementary components directly on a framed support bitmap.
int countMaskComponents(const std::vector<uint8_t>& mask, int rows, int columns, uint8_t value, bool fourConnected) {
    std::vector<uint8_t> seen(mask.size(), 0);
    std::vector<PixelId> pending;
    int components = 0;
    for (PixelId seed = 0; seed < rows * columns; ++seed) {
        if (seen[seed] || mask[seed] != value) {
            continue;
        }
        ++components;
        seen[seed] = 1;
        pending.push_back(seed);
        while (!pending.empty()) {
            const PixelId pixel = pending.back();
            pending.pop_back();
            for (int dr = -1; dr <= 1; ++dr) {
                for (int dc = -1; dc <= 1; ++dc) {
                    if ((dr == 0 && dc == 0) || (fourConnected && dr != 0 && dc != 0)) {
                        continue;
                    }
                    const int row = pixel / columns + dr;
                    const int column = pixel % columns + dc;
                    if (row < 0 || row >= rows || column < 0 || column >= columns) {
                        continue;
                    }
                    const PixelId next = row * columns + column;
                    if (!seen[next] && mask[next] == value) {
                        seen[next] = 1;
                        pending.push_back(next);
                    }
                }
            }
        }
    }
    return components;
}

std::pair<int, int> edgeVertices(ContourEdge edge, int columns) {
    const int stride = columns + 1;
    const int topLeft = (edge.pixel / columns) * stride + edge.pixel % columns;
    switch (edge.side) {
    case ContourSide::North: return {topLeft, topLeft + 1};
    case ContourSide::East: return {topLeft + 1, topLeft + stride + 1};
    case ContourSide::South: return {topLeft + stride + 1, topLeft + stride};
    case ContourSide::West: return {topLeft + stride, topLeft};
    }
    throw std::runtime_error("invalid test side");
}

template <class Trace>
void verifyBoundaryTopology(const MorphologicalTree& tree, const Trace& trace, NodeId node, bool fourConnected) {
    const int rows = tree.numRows() + 2;
    const int columns = tree.numColumns() + 2;
    const auto support = pixelsOfConnectedComponent(tree, node);
    std::vector<uint8_t> mask(static_cast<std::size_t>(rows * columns), 0);
    for (PixelId pixel : support) {
        mask[(pixel / tree.numColumns() + 1) * columns + pixel % tree.numColumns() + 1] = 1;
    }
    int external = 0, internal = 0, area2 = 0;
    std::vector<int> visitedEdges;
    for (const auto& boundary : trace.boundaries()) {
        external += boundary.kind == ContourBoundaryKind::External;
        internal += boundary.kind == ContourBoundaryKind::Internal;
        area2 += boundary.doubledSignedArea;
        std::vector<ContourEdge> edges;
        for (auto edge : trace.boundaryEdges(boundary)) {
            edges.push_back(edge);
            visitedEdges.push_back(packEdge(edge));
        }
        require(!edges.empty(), "a traced boundary must be nonempty");
        for (std::size_t i = 0; i < edges.size(); ++i) {
            requireEqual(edgeVertices(edges[i], tree.numColumns()).second,
                         edgeVertices(edges[(i + 1) % edges.size()], tree.numColumns()).first,
                         "ordered boundary must be continuous and closed");
        }
    }
    std::sort(visitedEdges.begin(), visitedEdges.end());
    requireVectorEqual(visitedEdges, expectedEdgesForNode(tree, node), "boundaries must partition the exact boundary without repetitions");
    requireEqual(area2, 2 * static_cast<int>(support.size()), "boundary areas must recover support area");
    requireEqual(external, countMaskComponents(mask, rows, columns, 1, fourConnected), "external boundaries must match foreground components");
    requireEqual(internal, countMaskComponents(mask, rows, columns, 0, !fourConnected) - 1, "internal boundaries must match complementary holes");
}

void verifyExhaustiveTernaryBoundaryTopology() {
    for (int pattern = 0; pattern < 19683; ++pattern) {
        auto image = ImageUInt8::create(3, 3, uint8_t{0});
        int code = pattern;
        for (PixelId pixel = 0; pixel < 9; ++pixel) {
            (*image)[pixel] = static_cast<uint8_t>(code % 3);
            code /= 3;
        }
        for (bool isMaxtree : {false, true}) {
            for (double radius : {1.0, 1.5}) {
                auto tree = makeComponentTree(image, isMaxtree, radius);
                auto traces = ContourTraceComputation(*tree);
                for (auto [node, trace] : traces) {
                    verifyBoundaryTopology(*tree, trace, node, radius == 1.0);
                }
            }
        }
    }
}

void verifyShapePolarityPairing() {
    for (auto immersion : {TestTopographicImmersion::Min4Max8, TestTopographicImmersion::Min8Max4}) {
        for (bool upper : {false, true}) {
            auto image = ImageUInt8::create(5, 5, uint8_t{1});
            (*image)[6] = (*image)[12] = upper ? 2 : 0;
            auto tree = MorphologicalTreeFactory::createTreeOfShapes<ToSGrayLevel>(image, makeTopographicConvention(image, immersion));
            const bool eightConnected = upper == (immersion == TestTopographicImmersion::Min4Max8);
            require((findNodeByArea(tree.topology(), 2) != InvalidNode) == eightConnected,
                    "shape fixture must join diagonal pixels exactly for the eight-connected polarity");
            auto traces = ContourTraceComputation(tree.asView());
            for (auto [node, trace] : traces) {
                const bool lower = !tree.topology().isRoot(node) && tree.asView().nodeAltitude(node) < tree.asView().nodeAltitude(tree.topology().parent(node));
                const bool four = immersion == TestTopographicImmersion::Min4Max8 ? lower : !lower;
                verifyBoundaryTopology(tree.topology(), trace, node, four);
                if (pixelsOfConnectedComponent(tree.topology(), node).size() == 2) {
                    auto topologyOnly = ContourTraceComputation(tree.topology());
                    requireThrows<std::invalid_argument>([&]() { static_cast<void>(topologyOnly.trace(node)); },
                                                        "ambiguous complementary shape requires a valued view");
                }
            }
        }
    }
}

void verifyTraversalMatchesNodeQueries() {
    auto image = makeComponentTreeFixture();

    for (bool isMaxtree : {true, false}) {
        auto tree = makeComponentTree(image, isMaxtree);
        auto traces = ContourTraceComputation(*tree);
        std::vector<std::vector<int>> traversalSignatures(static_cast<std::size_t>(tree->numInternalNodeSlots()));
        std::vector<NodeId> traversalOrder;
        for (auto [node, trace] : traces) {
            traversalSignatures[static_cast<std::size_t>(node)] = boundarySignature(trace);
            traversalOrder.push_back(node);
        }

        const std::string label = isMaxtree ? "max-tree" : "min-tree";
        requireVectorEqual(traversalOrder, collectNodeIds(tree->postOrder()), label + " trace traversal order");
        for (NodeId node : tree->aliveNodeIds()) {
            requireVectorEqual(traversalSignatures[static_cast<std::size_t>(node)], boundarySignature(traces.trace(node)),
                               label + " traversal and node query node " + std::to_string(node));
        }
    }
}

void verifyOwnedTraceRetainsItsStorage() {
    auto image = makeComponentTreeFixture();
    auto tree = makeComponentTree(image, true);
    auto traces = ContourTraceComputation(*tree);
    auto iterator = traces.begin();
    require(iterator != std::default_sentinel, "owned trace fixture must have a first trace");
    const auto [node, borrowedTrace] = *iterator;
    ContourTrace ownedTrace(borrowedTrace);
    const auto expectedSignature = boundarySignature(ownedTrace);

    while (iterator != std::default_sentinel) {
        ++iterator;
    }
    requireVectorEqual(boundarySignature(ownedTrace), expectedSignature, "owned trace after traversal completion");
    requireVectorEqual(edgeVector(ownedTrace), expectedEdgesForNode(*tree, node), "owned trace edges after traversal completion");
}

void verifyIndependentTraversalPositions() {
    auto image = makeComponentTreeFixture();
    auto tree = makeComponentTree(image, true);
    auto traces = ContourTraceComputation(*tree);
    auto first = traces.begin();
    auto independent = traces.begin();
    require(first != std::default_sentinel && independent != std::default_sentinel, "independent traversal fixture must be nonempty");
    requireEqual((*first).first, (*independent).first, "independent traversals begin at the same node");

    auto sharedPosition = first;
    ++first;
    requireEqual((*sharedPosition).first, (*first).first, "iterator copies share one traversal position");
    require((*independent).first != (*first).first, "independent begin calls have separate positions");
}

void verifyCallbackMatchesIteration() {
    auto image = makeComponentTreeFixture();
    auto tree = makeComponentTree(image, true);
    auto traces = ContourTraceComputation(*tree);
    std::vector<std::vector<int>> expected(static_cast<std::size_t>(tree->numInternalNodeSlots()));
    for (auto [node, trace] : traces) {
        expected[static_cast<std::size_t>(node)] = boundarySignature(trace);
    }
    int callbackCount = 0;
    traces.forEachTrace([&](NodeId node, ContourTraceView trace) {
        requireVectorEqual(boundarySignature(trace), expected[static_cast<std::size_t>(node)], "callback and iterator trace");
        ++callbackCount;
    });
    requireEqual(callbackCount, tree->numNodes(), "callback trace count");
}

} // namespace

int main() {
    auto image = makeComponentTreeFixture();

    for (bool isMaxtree : {true, false}) {
        auto tree = makeComponentTree(image, isMaxtree);
        verifyTraceEdgesAgainstSupportMasks(*tree, isMaxtree ? "max-tree" : "min-tree");

        if (isMaxtree && contract::validationsEnabled) {
            auto staleTraces = ContourTraceComputation(*tree);
            tree->mergeNodeIntoParent(4);
            requireThrows<std::logic_error>([&]() { static_cast<void>(staleTraces.begin()); },
                                            "trace traversal must reject topology mutation");
            requireThrows<std::logic_error>([&]() { static_cast<void>(staleTraces.trace(tree->root())); },
                                            "node trace must reject topology mutation");
        }

        auto valuedTree = makeValuedComponentTree(image, isMaxtree);
        auto topologyTraces = ContourTraceComputation(valuedTree->topology());
        auto viewTraces = ContourTraceComputation(valuedTree->asView());
        for (NodeId nodeId : valuedTree->topology().aliveNodeIds()) {
            requireVectorEqual(edgeVector(viewTraces, nodeId), edgeVector(topologyTraces, nodeId),
                               isMaxtree ? "max-tree traces via view" : "min-tree traces via view");
        }

        if (isMaxtree && contract::validationsEnabled) {
            auto staleValuedTree = makeValuedComponentTree(image, true);
            const auto staleView = staleValuedTree->asView();
            staleValuedTree->mergeNodeIntoParent(4);
            requireThrows<std::logic_error>([&]() { static_cast<void>(ContourTraceComputation(staleView)); },
                                            "trace extraction must reject stale ValuedMorphologicalTreeView");
        }
    }

    {
        auto tree = makeTreeOfShapes(image, TestTopographicImmersion::SelfDualSpan);
        verifyTraceEdgesAgainstSupportMasks(*tree, "tree of shapes");
    }

    verifyOnePixelBoundary();
    verifyRingBoundarySeparation();
    verifyBorderTouchingBoundary();
    verifyMultipleInternalBoundaries();
    verifyDiagonalTouchingBoundariesAreDeterministic();
    verifyClosureAndComplementaryPairing();
    verifyExhaustiveTernaryBoundaryTopology();
    verifyShapePolarityPairing();
    verifyTraversalMatchesNodeQueries();
    verifyOwnedTraceRetainsItsStorage();
    verifyIndependentTraversalPositions();
    verifyCallbackMatchesIteration();

    return 0;
}
