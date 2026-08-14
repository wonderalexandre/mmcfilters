#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/AttributeComputation.hpp"
#include "mmcfilters/contours/ContourTraceComputation.hpp"
#include "mmcfilters/contours/ContoursComputedIncrementally.hpp"

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

using IncrementalContourTraces = ContourTraceComputation::IncrementalContourTraces;

static_assert(std::is_same_v<decltype(std::declval<const IncrementalContourTraces&>().getLoops(NodeId{})), std::vector<ContourTraceLoop>>);

std::vector<PixelId> pixelsOfConnectedComponent(const MorphologicalTree& tree, NodeId nodeId) {
    std::vector<PixelId> pixels;
    for (NodeId subtreeNodeId : tree.subtreeNodes(nodeId)) {
        for (PixelId pixel : tree.properPart(subtreeNodeId)) {
            pixels.push_back(pixel);
        }
    }
    return pixels;
}

int sideIndex(ContourTraceSide side) { return static_cast<int>(side); }

int packEdge(const ContourTraceEdge& edge) { return ContourTraceComputation::packEdge(edge.pixel, edge.side); }

PixelId neighborPixel(const MorphologicalTree& tree, PixelId pixel, ContourTraceSide side) {
    const int rows = tree.numRows();
    const int columns = tree.numColumns();
    const auto [row, column] = ImageUtils::to2D(pixel, columns);
    switch (side) {
    case ContourTraceSide::North:
        return row == 0 ? InvalidPixel : ImageUtils::to1D(row - 1, column, columns);
    case ContourTraceSide::West:
        return column == 0 ? InvalidPixel : ImageUtils::to1D(row, column - 1, columns);
    case ContourTraceSide::East:
        return column == columns - 1 ? InvalidPixel : ImageUtils::to1D(row, column + 1, columns);
    case ContourTraceSide::South:
        return row == rows - 1 ? InvalidPixel : ImageUtils::to1D(row + 1, column, columns);
    }
    return InvalidPixel;
}

std::vector<int> expectedEdgesForNode(const MorphologicalTree& tree, NodeId nodeId) {
    static constexpr std::array<ContourTraceSide, 4> sides{ContourTraceSide::North, ContourTraceSide::West, ContourTraceSide::East, ContourTraceSide::South};

    std::vector<uint8_t> mask(static_cast<std::size_t>(tree.numPixels()), 0);
    for (PixelId pixel : pixelsOfConnectedComponent(tree, nodeId)) {
        mask[static_cast<std::size_t>(pixel)] = 1;
    }

    std::vector<int> expected;
    for (PixelId pixel = 0; pixel < tree.numPixels(); ++pixel) {
        if (!mask[static_cast<std::size_t>(pixel)]) {
            continue;
        }
        for (ContourTraceSide side : sides) {
            const PixelId neighbor = neighborPixel(tree, pixel, side);
            if (neighbor == InvalidPixel || !mask[static_cast<std::size_t>(neighbor)]) {
                expected.push_back(ContourTraceComputation::packEdge(pixel, side));
            }
        }
    }
    std::sort(expected.begin(), expected.end());
    return expected;
}

std::vector<int> edgeVector(const ContourTraceComputation::IncrementalContourTraces& traces, NodeId nodeId) {
    std::vector<int> values;
    for (ContourTraceEdge edge : traces.getEdges(nodeId)) {
        values.push_back(packEdge(edge));
    }
    std::sort(values.begin(), values.end());
    return values;
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

std::vector<PixelId> contourVector(const ContoursComputedIncrementally::IncrementalContours& contours, NodeId nodeId) {
    std::vector<PixelId> values;
    for (PixelId pixel : contours.getContour(nodeId)) {
        values.push_back(pixel);
    }
    std::sort(values.begin(), values.end());
    return values;
}

std::vector<int> loopSignature(const ContourTraceComputation::IncrementalContourTraces& traces, NodeId nodeId) {
    std::vector<std::array<int, 4>> loopSummaries;
    for (const ContourTraceLoop& loop : traces.getLoops(nodeId)) {
        int edgeSum = 0;
        for (ContourTraceEdge edge : traces.getLoopEdges(loop)) {
            edgeSum += packEdge(edge) + 1;
        }
        loopSummaries.push_back(std::array<int, 4>{static_cast<int>(loop.kind), static_cast<int>(loop.edgeCount), loop.signedArea2, edgeSum});
    }
    std::sort(loopSummaries.begin(), loopSummaries.end());

    std::vector<int> signature;
    signature.reserve(loopSummaries.size() * 4);
    for (const auto& summary : loopSummaries) {
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
    auto traces = ContourTraceComputation::extract(tree);
    auto contours = ContoursComputedIncrementally::extractCompactContours(tree);
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
                     counts[static_cast<std::size_t>(sideIndex(ContourTraceSide::North))], label + " north side count node " + std::to_string(nodeId));
        requireEqual(static_cast<int>(values[static_cast<std::size_t>(names.linearIndex(nodeId, ContourSideWest))]),
                     counts[static_cast<std::size_t>(sideIndex(ContourTraceSide::West))], label + " west side count node " + std::to_string(nodeId));
        requireEqual(static_cast<int>(values[static_cast<std::size_t>(names.linearIndex(nodeId, ContourSideEast))]),
                     counts[static_cast<std::size_t>(sideIndex(ContourTraceSide::East))], label + " east side count node " + std::to_string(nodeId));
        requireEqual(static_cast<int>(values[static_cast<std::size_t>(names.linearIndex(nodeId, ContourSideSouth))]),
                     counts[static_cast<std::size_t>(sideIndex(ContourTraceSide::South))], label + " south side count node " + std::to_string(nodeId));
    }
}

void verifyOnePixelLoop() {
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

    auto traces = ContourTraceComputation::extract(*tree);
    auto loops = traces.getLoops(nodeId);
    requireEqual(loops.size(), std::size_t{1}, "one-pixel support loop count");
    requireEqual(static_cast<int>(loops[0].kind), static_cast<int>(ContourLoopKind::External), "one-pixel support loop kind");
    requireEqual(loops[0].edgeCount, uint32_t{4}, "one-pixel support loop edge count");
    requireEqual(loops[0].signedArea2, 2, "one-pixel support signed area");
    requireEqual(traces.getLoopEdges(loops[0]).size(), std::size_t{4}, "one-pixel support loop edge range");
}

void verifyRingLoopSeparation() {
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

    auto traces = ContourTraceComputation::extract(*tree);
    auto loops = traces.getLoops(ringNode);
    requireEqual(loops.size(), std::size_t{2}, "ring support loop count");

    int externalLoops = 0;
    int internalLoops = 0;
    int externalEdges = 0;
    int internalEdges = 0;
    int totalEdges = 0;
    for (const ContourTraceLoop& loop : loops) {
        totalEdges += static_cast<int>(loop.edgeCount);
        if (loop.kind == ContourLoopKind::External) {
            ++externalLoops;
            externalEdges += static_cast<int>(loop.edgeCount);
            require(loop.signedArea2 > 0, "ring external loop must have positive signed area");
        } else {
            ++internalLoops;
            internalEdges += static_cast<int>(loop.edgeCount);
            require(loop.signedArea2 < 0, "ring internal loop must have negative signed area");
        }
    }

    requireEqual(externalLoops, 1, "ring external loop count");
    requireEqual(internalLoops, 1, "ring internal loop count");
    requireEqual(externalEdges, 12, "ring external edge count");
    requireEqual(internalEdges, 4, "ring internal edge count");
    requireEqual(totalEdges, 16, "ring total edge count");
}

void requireLoopSummary(const MorphologicalTree& tree, NodeId nodeId, int expectedExternalLoops, int expectedInternalLoops, int expectedExternalEdges,
                        int expectedInternalEdges, const std::string& label) {
    auto traces = ContourTraceComputation::extract(tree);
    auto loops = traces.getLoops(nodeId);

    int externalLoops = 0;
    int internalLoops = 0;
    int externalEdges = 0;
    int internalEdges = 0;
    for (const ContourTraceLoop& loop : loops) {
        if (loop.kind == ContourLoopKind::External) {
            ++externalLoops;
            externalEdges += static_cast<int>(loop.edgeCount);
            require(loop.signedArea2 > 0, label + " external loop signed area");
        } else {
            ++internalLoops;
            internalEdges += static_cast<int>(loop.edgeCount);
            require(loop.signedArea2 < 0, label + " internal loop signed area");
        }
    }

    requireEqual(externalLoops, expectedExternalLoops, label + " external loop count");
    requireEqual(internalLoops, expectedInternalLoops, label + " internal loop count");
    requireEqual(externalEdges, expectedExternalEdges, label + " external edge count");
    requireEqual(internalEdges, expectedInternalEdges, label + " internal edge count");
    requireEqual(externalEdges + internalEdges, static_cast<int>(edgeVector(traces, nodeId).size()), label + " loop edges must cover materialized edges");
}

void verifyBorderTouchingLoop() {
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

    requireLoopSummary(*tree, nodeId, 1, 0, 8, 0, "border-touching support");
}

void verifyMultipleInternalLoops() {
    auto image = makeImage(5, 5,
                           {
                               2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2,
                           });
    auto tree = makeComponentTree(image, true);
    const NodeId nodeId = findNodeByArea(*tree, 23);
    require(nodeId != InvalidNode, "two-hole fixture must have an area-23 node");

    requireLoopSummary(*tree, nodeId, 1, 2, 20, 8, "two-hole support");
}

void verifyDiagonalTouchingLoopsAreDeterministic() {
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

    requireLoopSummary(*tree, nodeId, 2, 0, 8, 0, "diagonal-touching support");
}

void verifyLoopAccessOrderIndependence() {
    auto image = makeComponentTreeFixture();

    for (bool isMaxtree : {true, false}) {
        auto tree = makeComponentTree(image, isMaxtree);
        auto baseline = ContourTraceComputation::extract(*tree);
        baseline.materializeAll();

        auto lazy = ContourTraceComputation::extract(*tree);
        std::vector<NodeId> reverseOrder = collectNodeIds(tree->aliveNodeIds());
        std::reverse(reverseOrder.begin(), reverseOrder.end());
        for (NodeId nodeId : reverseOrder) {
            static_cast<void>(lazy.getLoops(nodeId));
        }

        const std::string label = isMaxtree ? "max-tree" : "min-tree";
        for (NodeId nodeId : tree->aliveNodeIds()) {
            requireVectorEqual(loopSignature(lazy, nodeId), loopSignature(baseline, nodeId),
                               label + " trace loop access-order independence node " + std::to_string(nodeId));
        }
    }
}

void verifyLoopResultsOwnTheirStorage() {
    auto image = makeComponentTreeFixture();
    auto tree = makeComponentTree(image, true);
    auto traces = ContourTraceComputation::extract(*tree);
    std::vector<NodeId> nodes = collectNodeIds(tree->aliveNodeIds());
    require(nodes.size() > 1, "owned-loop fixture must have multiple nodes");

    auto heldLoops = traces.getLoops(nodes.front());
    require(!heldLoops.empty(), "owned-loop fixture first node must expose a loop");
    const ContourTraceLoop expected = heldLoops.front();

    for (std::size_t i = 1; i < nodes.size(); ++i) {
        static_cast<void>(traces.getLoops(nodes[i]));
    }

    requireEqual(static_cast<int>(heldLoops.front().kind), static_cast<int>(expected.kind), "owned loop kind after lazy materialization");
    requireEqual(heldLoops.front().edgeOffset, expected.edgeOffset, "owned loop edge offset after lazy materialization");
    requireEqual(heldLoops.front().edgeCount, expected.edgeCount, "owned loop edge count after lazy materialization");
    requireEqual(heldLoops.front().signedArea2, expected.signedArea2, "owned loop signed area after lazy materialization");
    requireEqual(traces.getLoopEdges(heldLoops.front()).size(), static_cast<std::size_t>(expected.edgeCount),
                 "owned loop must remain usable after lazy materialization");
}

void verifyGetLoopsMaterializesOnlyRequestedNode() {
    auto image = makeComponentTreeFixture();

    for (bool isMaxtree : {true, false}) {
        auto tree = makeComponentTree(image, isMaxtree);
        auto traces = ContourTraceComputation::extract(*tree);
        const NodeId root = tree->root();

        static_cast<void>(traces.getLoops(root));
        require(traces.isNodeTraced(root), "getLoops(root) must trace root");

        int tracedNodes = 0;
        int liveNodes = 0;
        for (NodeId nodeId : tree->aliveNodeIds()) {
            ++liveNodes;
            if (traces.isNodeTraced(nodeId)) {
                ++tracedNodes;
            }
        }
        require(liveNodes > 1, "lazy loop materialization fixture must have more than one node");
        requireEqual(tracedNodes, 1, isMaxtree ? "max-tree getLoops(root) must not trace descendants" : "min-tree getLoops(root) must not trace descendants");

        traces.materializeAll();
        require(traces.isMaterialized(), "materializeAll must still trace every live node");
    }
}

void verifyScratchReleaseAfterGlobalEdgeMaterialization() {
    auto image = makeComponentTreeFixture();
    auto tree = makeComponentTree(image, true);
    auto traces = ContourTraceComputation::extract(*tree);

    auto extractedStats = traces.storageStats();
    require(extractedStats.addDeltaValues > 0, "trace fixture must have addition deltas before materialization");
    require(extractedStats.removeDeltaValues > 0, "trace fixture must have removal deltas before materialization");

    static_cast<void>(traces.getEdges(tree->root()));
    auto edgeStats = traces.storageStats();
    requireEqual(edgeStats.cachedEdgeReadyNodes, static_cast<std::size_t>(tree->numNodes()), "root edge materialization must prepare every live edge cache");
    requireEqual(edgeStats.addDeltaValues, std::size_t{0}, "edge deltas must be released after all edges are ready");
    requireEqual(edgeStats.removeDeltaValues, std::size_t{0}, "removal deltas must be released after all edges are ready");

    traces.materializeAll();
    require(traces.isMaterialized(), "scratch release must preserve full loop materialization");
    for (NodeId nodeId : tree->aliveNodeIds()) {
        static_cast<void>(traces.getLoops(nodeId));
    }
}

void verifyNodeLocalLoopTracingUsesSparseAdjacency() {
    auto image = ImageUInt8::create(64, 64, std::uint8_t{0});
    (*image)[ImageUtils::to1D(32, 32, 64)] = std::uint8_t{1};
    auto tree = makeComponentTree(image, true);
    require(tree->numNodes() > 1, "sparse trace fixture must have more than one node");

    auto traces = ContourTraceComputation::extract(*tree);
    static_cast<void>(traces.getLoops(tree->root()));
    const auto rootLoopStats = traces.storageStats();
    requireEqual(rootLoopStats.traceDenseOutgoingSlots, std::size_t{0}, "node-local root loop tracing must not allocate dense outgoing heads");
    require(rootLoopStats.traceSparseOutgoingSlots > 0, "node-local root loop tracing must allocate sparse outgoing heads");
    require(traces.isNodeTraced(tree->root()), "sparse node-local trace must mark root traced");
    require(!traces.isMaterialized(), "sparse node-local trace must not trace every node");

    traces.materializeAll();
    const auto finalStats = traces.storageStats();
    requireEqual(finalStats.traceDenseOutgoingSlots, std::size_t{0}, "completed global loop tracing must release dense outgoing heads");
    requireEqual(finalStats.traceSparseOutgoingSlots, std::size_t{0}, "completed global loop tracing must release sparse outgoing heads");
}

void verifyRepeatedNodeLocalLoopTracingSwitchesToDenseAdjacency() {
    auto image = ImageUInt8::create(64, 64, std::uint8_t{0});
    for (int i = 0; i < 12; ++i) {
        const int row = 4 + (i / 4) * 16;
        const int column = 4 + (i % 4) * 16;
        (*image)[ImageUtils::to1D(row, column, 64)] = static_cast<std::uint8_t>(i + 1);
    }
    auto tree = makeComponentTree(image, true);
    std::vector<NodeId> nodes = collectNodeIds(tree->aliveNodeIds());
    require(nodes.size() > 10, "sparse-to-dense fixture must expose enough nodes");

    auto traces = ContourTraceComputation::extract(*tree);
    for (std::size_t i = 0; i < 8; ++i) {
        static_cast<void>(traces.getLoops(nodes[i]));
    }
    const auto sparseStats = traces.storageStats();
    requireEqual(sparseStats.traceDenseOutgoingSlots, std::size_t{0}, "first node-local loop queries should keep dense outgoing heads unallocated");
    require(sparseStats.traceSparseOutgoingSlots > 0, "first node-local loop queries should use sparse outgoing heads");

    static_cast<void>(traces.getLoops(nodes[8]));
    const auto denseStats = traces.storageStats();
    require(denseStats.traceDenseOutgoingSlots > 0, "repeated node-local loop queries should switch to dense outgoing heads");
    requireEqual(denseStats.traceSparseOutgoingSlots, std::size_t{0}, "switching to dense outgoing heads should release sparse scratch");
}

} // namespace

int main() {
    auto image = makeComponentTreeFixture();

    for (bool isMaxtree : {true, false}) {
        auto tree = makeComponentTree(image, isMaxtree);
        verifyTraceEdgesAgainstSupportMasks(*tree, isMaxtree ? "max-tree" : "min-tree");

        if (isMaxtree && contract::validationsEnabled) {
            auto staleTraces = ContourTraceComputation::extract(*tree);
            tree->mergeNodeIntoParent(4);
            requireThrows<std::logic_error>([&]() { static_cast<void>(staleTraces.isMaterialized()); },
                                            "trace materialization status must reject topology mutation");
            requireThrows<std::logic_error>([&]() { static_cast<void>(staleTraces.getEdges(tree->root())); },
                                            "trace edge access must reject topology mutation");
            requireThrows<std::logic_error>([&]() { staleTraces.materializeAll(); }, "trace materializeAll must reject topology mutation");
        }

        auto valuedTree = makeValuedComponentTree(image, isMaxtree);
        auto topologyTraces = ContourTraceComputation::extract(valuedTree->topology());
        auto viewTraces = ContourTraceComputation::extract(valuedTree->asView());
        for (NodeId nodeId : valuedTree->topology().aliveNodeIds()) {
            requireVectorEqual(edgeVector(viewTraces, nodeId), edgeVector(topologyTraces, nodeId),
                               isMaxtree ? "max-tree traces via view" : "min-tree traces via view");
        }

        if (isMaxtree && contract::validationsEnabled) {
            auto staleValuedTree = makeValuedComponentTree(image, true);
            const auto staleView = staleValuedTree->asView();
            staleValuedTree->mergeNodeIntoParent(4);
            requireThrows<std::logic_error>([&]() { static_cast<void>(ContourTraceComputation::extract(staleView)); },
                                            "trace extraction must reject stale ValuedMorphologicalTreeView");
        }
    }

    {
        auto tree = makeTreeOfShapes(image, TestTopographicImmersion::SelfDualSpan);
        verifyTraceEdgesAgainstSupportMasks(*tree, "tree of shapes");
    }

    verifyOnePixelLoop();
    verifyRingLoopSeparation();
    verifyBorderTouchingLoop();
    verifyMultipleInternalLoops();
    verifyDiagonalTouchingLoopsAreDeterministic();
    verifyLoopAccessOrderIndependence();
    verifyLoopResultsOwnTheirStorage();
    verifyGetLoopsMaterializesOnlyRequestedNode();
    verifyScratchReleaseAfterGlobalEdgeMaterialization();
    verifyNodeLocalLoopTracingUsesSparseAdjacency();
    verifyRepeatedNodeLocalLoopTracingSwitchesToDenseAdjacency();

    return 0;
}
