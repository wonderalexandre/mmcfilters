#include "support/TestSupport.hpp"

#include "mmcfilters/trees/saliency/HierarchySaliencyMapValidation.hpp"
#include "mmcfilters/trees/saliency/HierarchySaliencyMapProjection.hpp"
#include "mmcfilters/trees/saliency/HierarchySaliencyMap.hpp"

#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <utility>
#include <vector>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

namespace {

ValuedMorphologicalTree<std::uint8_t> makeThreePixelValuedTree(MorphologicalTreeKind kind, std::optional<RegularGridAdjacency2D> adjacency) {
    const std::vector<NodeId> parent = {3, 4, 4, 5, 5, 5};
    const std::vector<std::uint8_t> altitude = {5, 4, 4, 5, 4, 0};
    return MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(parent), std::span<const std::uint8_t>(altitude), 1, 3, kind,
                                                           std::move(adjacency));
}

ValuedMorphologicalTree<std::uint8_t> makeThreePixelMinTree() {
    const std::vector<NodeId> parent = {3, 4, 4, 5, 5, 5};
    const std::vector<std::uint8_t> altitude = {0, 1, 1, 0, 1, 5};
    return MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(parent), std::span<const std::uint8_t>(altitude), 1, 3,
                                                           MorphologicalTreeKind::MinTree, RegularGridAdjacency2D(1, 3, 1.0));
}

ValuedMorphologicalTree<std::uint8_t> makeTwoSingletonHierarchy() {
    const std::vector<NodeId> nodeParent{2, 2, 2};
    const std::vector<NodeId> smallestNodeMap{0, 1};
    const std::vector<std::uint8_t> altitude{2, 1, 0};
    return MorphologicalTreeFactory::createFromNativeTopology(std::span<const NodeId>(nodeParent), std::span<const NodeId>(smallestNodeMap),
                                                              std::span<const std::uint8_t>(altitude), 2, 1, 2,
                                                              makeMorphologicalTreeSemantics(
                                                                  MorphologicalTreeKind::MaxTree,
                                                                  SharedAdjacencyContext{RegularGridAdjacency2D(1, 2, 1.0)}));
}

ValuedMorphologicalTree<std::uint8_t> makeSpatiallyDisconnectedHierarchy() {
    const std::vector<NodeId> nodeParent{1, 1};
    const std::vector<NodeId> smallestNodeMap{0, 1, 0};
    const std::vector<std::uint8_t> altitude{1, 0};
    return MorphologicalTreeFactory::createFromNativeTopology(std::span<const NodeId>(nodeParent), std::span<const NodeId>(smallestNodeMap),
                                                              std::span<const std::uint8_t>(altitude), 1, 1, 3,
                                                              makeMorphologicalTreeSemantics(
                                                                  MorphologicalTreeKind::MaxTree,
                                                                  SharedAdjacencyContext{RegularGridAdjacency2D(1, 3, 1.0)}));
}

ValuedMorphologicalTree<std::uint8_t> makeDirectionalThreePixelTree(bool equivalentDirections) {
    const std::vector<NodeId> nodeParent{2, 2, 2};
    const std::vector<NodeId> smallestNodeMap{0, 1, 1};
    const std::vector<std::uint8_t> altitude{0, 20, 10};
    return MorphologicalTreeFactory::createFromNativeTopology(
        std::span<const NodeId>(nodeParent), std::span<const NodeId>(smallestNodeMap), std::span<const std::uint8_t>(altitude), 2, 1, 3,
        MorphologicalTreeSemantics{
            MorphologicalTreeKind::TreeOfShapes, NodeAltitudeOrder::Unconstrained,
            TopographicConvention{ComplementaryGridImmersion{ComplementaryAdjacencies{
                                      RegularGridAdjacency2D(1, 3, 1.0), RegularGridAdjacency2D(1, 3, equivalentDirections ? 1.0 : 1.5)}},
                                  TopographicDomainExtension::ExteriorRing, PixelId{0}}});
}

template <class Value> std::vector<int> computeQfzLabels(const EdgeSaliencyMap<Value>& edgeMap, Value threshold) {
    const int numVertices = edgeMap.numRows * edgeMap.numColumns;
    std::vector<int> parent(static_cast<std::size_t>(numVertices));
    for (int i = 0; i < numVertices; ++i) {
        parent[static_cast<std::size_t>(i)] = i;
    }

    auto findRoot = [&parent](int node) {
        int root = node;
        while (parent[static_cast<std::size_t>(root)] != root) {
            root = parent[static_cast<std::size_t>(root)];
        }
        while (parent[static_cast<std::size_t>(node)] != node) {
            const int next = parent[static_cast<std::size_t>(node)];
            parent[static_cast<std::size_t>(node)] = root;
            node = next;
        }
        return root;
    };

    for (std::size_t i = 0; i < edgeMap.size(); ++i) {
        if (edgeMap.values[i] < threshold) {
            const int sourceRoot = findRoot(edgeMap.sources[i]);
            const int targetRoot = findRoot(edgeMap.targets[i]);
            if (sourceRoot != targetRoot) {
                parent[static_cast<std::size_t>(targetRoot)] = sourceRoot;
            }
        }
    }

    std::vector<int> labels(static_cast<std::size_t>(numVertices));
    for (int i = 0; i < numVertices; ++i) {
        labels[static_cast<std::size_t>(i)] = findRoot(i);
    }
    return labels;
}

} // namespace

int main() {
    auto valuedTree = makeThreePixelValuedTree(MorphologicalTreeKind::MaxTree, RegularGridAdjacency2D(1, 3, 1.0));

    std::vector<float> nodeValuation(static_cast<std::size_t>(valuedTree.topology().numInternalNodeSlots()), 0.0f);
    nodeValuation[0] = 10.0f;
    nodeValuation[1] = 20.0f;
    nodeValuation[2] = 30.0f;
    auto scoreMap = HierarchySaliencyMap::computeSaliencyEdgeMap(valuedTree.topology(), std::span<const float>(nodeValuation));
    requireVectorEqual(scoreMap.sources, {0, 1}, "valuation saliency sources");
    requireVectorEqual(scoreMap.targets, {1, 2}, "valuation saliency targets");
    requireVectorEqual(scoreMap.values, {30.0f, 0.0f}, "valuation saliency values");

    auto topologicalMap = HierarchySaliencyMap::computeTopologicalLevelEdgeMap(valuedTree.topology());
    requireVectorEqual(topologicalMap.values, {1, 0}, "topological-level saliency values");

    std::vector<int> topologicalValuation = HierarchySaliencyMap::computeTopologicalLevels(valuedTree.topology());
    HierarchySaliencyMapValidation::validateHierarchyValuation(valuedTree.topology(), std::span<const int>(topologicalValuation),
                                                               HierarchyValuationPolicy::RequireStrictHierarchy);
    auto strictTopologicalMap = HierarchySaliencyMap::computeSaliencyEdgeMap(valuedTree.topology(), std::span<const int>(topologicalValuation),
                                                                             HierarchyValuationPolicy::RequireStrictHierarchy);
    requireVectorEqual(strictTopologicalMap.values, {1, 0}, "strict topological saliency values");
    requireVectorEqual(computeQfzLabels(strictTopologicalMap, 0), {0, 1, 2}, "qfz lambda 0 recovers singleton partition");
    requireVectorEqual(computeQfzLabels(strictTopologicalMap, 1), {0, 1, 1}, "qfz lambda 1 recovers intermediate partition");
    requireVectorEqual(computeQfzLabels(strictTopologicalMap, 2), {0, 0, 0}, "qfz lambda 2 recovers root partition");

    std::vector<int> rankedScores = HierarchySaliencyMapValidation::rankHierarchyValuation(valuedTree.topology(), std::span<const float>(nodeValuation),
                                                                                           HierarchyValuationPolicy::RequireStrictHierarchy);
    requireVectorEqual(rankedScores, {0, 1, 2}, "ranked strict hierarchy valuation");
    std::vector<double> normalizedScores = HierarchySaliencyMapValidation::computeNormalizedScores(valuedTree.topology(), std::span<const float>(nodeValuation),
                                                                                                   HierarchyValuationPolicy::RequireStrictHierarchy);
    requireNear(normalizedScores[0], 0.0, 1e-12, "normalized strict valuation leaf level");
    requireNear(normalizedScores[1], 0.5, 1e-12, "normalized strict valuation intermediate level");
    requireNear(normalizedScores[2], 1.0, 1e-12, "normalized strict valuation root level");
    auto rankedScoreMap =
        HierarchySaliencyMap::computeSaliencyEdgeMap(valuedTree.topology(), std::span<const int>(rankedScores), HierarchyValuationPolicy::RequireStrictHierarchy);
    requireVectorEqual(rankedScoreMap.values, {2, 0}, "ranked valuation saliency values");

    auto canonicalRankedScoreMap = HierarchySaliencyMap::computeCanonicalRankedSaliencyEdgeMap(valuedTree.topology(), std::span<const float>(nodeValuation),
                                                                                               HierarchyValuationPolicy::RequireStrictHierarchy);
    requireVectorEqual(canonicalRankedScoreMap.values, {1, 0}, "canonical ranks use only effective edge levels");
    requireVectorEqual(computeQfzLabels(canonicalRankedScoreMap, 0), {0, 1, 2}, "canonical qfz lambda 0 is singleton partition");
    requireVectorEqual(computeQfzLabels(canonicalRankedScoreMap, 1), {0, 1, 1}, "canonical qfz lambda 1 is intermediate partition");
    requireVectorEqual(computeQfzLabels(canonicalRankedScoreMap, 2), {0, 0, 0}, "canonical qfz lambda 2 is root partition");

    std::vector<int> appearanceLevels = HierarchySaliencyMap::computePartitionAppearanceLevels(valuedTree.topology());
    auto appearanceMap =
        HierarchySaliencyMap::computeSaliencyEdgeMap(valuedTree.topology(), std::span<const int>(appearanceLevels),
                                                     HierarchyValuationPolicy::RequireStrictHierarchy, HierarchyLevelConvention::PartitionAppearanceLevel);
    requireVectorEqual(appearanceMap.values, {1, 0}, "partition appearance convention applies Cousty level(LCA)-1");

    auto singletonHierarchy = makeTwoSingletonHierarchy();
    auto singletonTopologicalMap = HierarchySaliencyMap::computeTopologicalLevelEdgeMap(singletonHierarchy.topology());
    requireVectorEqual(singletonTopologicalMap.values, {0}, "canonical topological scale removes a redundant singleton-leaf level");

    auto disconnectedHierarchy = makeSpatiallyDisconnectedHierarchy();
    requireThrowsContaining<std::invalid_argument>(
        [&]() { HierarchySaliencyMapValidation::validateHierarchyConnectivity(disconnectedHierarchy.topology(), RegularGridAdjacency2D(1, 3, 1.0)); },
        "disconnected support", "formal hierarchy connectivity rejects a spatially disconnected node support");
    std::vector<int> disconnectedLevels = HierarchySaliencyMap::computeTopologicalLevels(disconnectedHierarchy.topology());
    requireThrowsContaining<std::invalid_argument>(
        [&]() { static_cast<void>(HierarchySaliencyMap::computeSaliencyEdgeMap(disconnectedHierarchy.topology(), std::span<const int>(disconnectedLevels))); },
        "disconnected support", "formal saliency projection validates graph connectivity by default");

    std::vector<float> collapsingScores = {30.0f, 30.0f, 30.0f};
    HierarchySaliencyMapValidation::validateHierarchyValuation(valuedTree.topology(), std::span<const float>(collapsingScores),
                                                               HierarchyValuationPolicy::AllowLevelCollapse);
    requireThrows<std::invalid_argument>(
        [&]() {
            HierarchySaliencyMapValidation::validateHierarchyValuation(valuedTree.topology(), std::span<const float>(collapsingScores),
                                                                       HierarchyValuationPolicy::RequireStrictHierarchy);
        },
        "strict hierarchy valuation must reject equal parent-child levels");

    std::vector<float> negativeScores = {-10.0f, -5.0f, 0.0f};
    HierarchySaliencyMapValidation::validateHierarchyValuation(valuedTree.topology(), std::span<const float>(negativeScores));
    requireThrows<std::invalid_argument>(
        [&]() { static_cast<void>(HierarchySaliencyMap::computeSaliencyEdgeMap(valuedTree.topology(), std::span<const float>(negativeScores))); },
        "formal saliency projection must reject negative valuations");
    requireThrows<std::invalid_argument>(
        [&]() {
            HierarchySaliencyMapValidation::validateHierarchyValuation(valuedTree.topology(), std::span<const float>(negativeScores),
                                                                       HierarchyValuationPolicy::AllowLevelCollapse,
                                                                       HierarchyValuationRangePolicy::RequireNonNegative);
        },
        "non-negative hierarchy valuation policy must reject negative values");
    std::vector<int> rankedNegativeScores = HierarchySaliencyMapValidation::rankHierarchyValuation(valuedTree.topology(), std::span<const float>(negativeScores));
    requireVectorEqual(rankedNegativeScores, {0, 1, 2}, "ranked negative valuation becomes non-negative");
    std::vector<double> normalizedNegativeScores =
        HierarchySaliencyMapValidation::computeNormalizedScores(valuedTree.topology(), std::span<const float>(negativeScores));
    requireNear(normalizedNegativeScores[0], 0.0, 1e-12, "normalized negative valuation leaf level");
    requireNear(normalizedNegativeScores[1], 0.5, 1e-12, "normalized negative valuation intermediate level");
    requireNear(normalizedNegativeScores[2], 1.0, 1e-12, "normalized negative valuation root level");
    requireThrows<std::invalid_argument>(
        [&]() {
            static_cast<void>(HierarchySaliencyMapValidation::computeNormalizedScores(valuedTree.topology(), std::span<const float>(negativeScores),
                                                                                      HierarchyValuationPolicy::AllowLevelCollapse,
                                                                                      HierarchyValuationRangePolicy::RequireNonNegative));
        },
        "non-negative normalized valuation policy must reject negative values");

    std::vector<std::int64_t> extremeIntegerScores = {std::numeric_limits<std::int64_t>::min(), 0, std::numeric_limits<std::int64_t>::max()};
    std::vector<double> normalizedExtremeIntegerScores = HierarchySaliencyMapValidation::computeNormalizedScores(
        valuedTree.topology(), std::span<const std::int64_t>(extremeIntegerScores), HierarchyValuationPolicy::RequireStrictHierarchy);
    requireNear(normalizedExtremeIntegerScores[0], 0.0, 1e-12, "normalized integer minimum");
    requireNear(normalizedExtremeIntegerScores[1], 0.5, 1e-12, "normalized integer midpoint");
    requireNear(normalizedExtremeIntegerScores[2], 1.0, 1e-12, "normalized integer maximum");

    const double maxFiniteDouble = std::numeric_limits<double>::max();
    std::vector<double> extremeDoubleScores = {-maxFiniteDouble, 0.0, maxFiniteDouble};
    std::vector<double> normalizedExtremeDoubleScores = HierarchySaliencyMapValidation::computeNormalizedScores(
        valuedTree.topology(), std::span<const double>(extremeDoubleScores), HierarchyValuationPolicy::RequireStrictHierarchy);
    requireNear(normalizedExtremeDoubleScores[0], 0.0, 1e-12, "normalized finite double minimum");
    requireNear(normalizedExtremeDoubleScores[1], 0.5, 1e-12, "normalized finite double midpoint");
    requireNear(normalizedExtremeDoubleScores[2], 1.0, 1e-12, "normalized finite double maximum");
    for (double score : normalizedExtremeDoubleScores) {
        require(std::isfinite(score) && score >= 0.0 && score <= 1.0, "normalized finite double extremes must remain finite and inside [0, 1]");
    }

    const long double maxFiniteLongDouble = std::numeric_limits<long double>::max();
    std::vector<long double> extremeLongDoubleScores = {-maxFiniteLongDouble, 0.0L, maxFiniteLongDouble};
    std::vector<double> normalizedExtremeLongDoubleScores = HierarchySaliencyMapValidation::computeNormalizedScores(
        valuedTree.topology(), std::span<const long double>(extremeLongDoubleScores), HierarchyValuationPolicy::RequireStrictHierarchy);
    requireNear(normalizedExtremeLongDoubleScores[0], 0.0, 1e-12, "normalized finite long double minimum");
    requireNear(normalizedExtremeLongDoubleScores[1], 0.5, 1e-12, "normalized finite long double midpoint");
    requireNear(normalizedExtremeLongDoubleScores[2], 1.0, 1e-12, "normalized finite long double maximum");
    for (double score : normalizedExtremeLongDoubleScores) {
        require(std::isfinite(score) && score >= 0.0 && score <= 1.0, "normalized finite long double extremes must remain finite and inside [0, 1]");
    }

    std::vector<float> invalidScores = {10.0f, 40.0f, 30.0f};
    requireThrows<std::invalid_argument>(
        [&]() { static_cast<void>(HierarchySaliencyMap::computeSaliencyEdgeMap(valuedTree.topology(), std::span<const float>(invalidScores))); },
        "formal saliency must reject valuations that decrease toward the root");
    requireThrows<std::invalid_argument>(
        [&]() { static_cast<void>(HierarchySaliencyMapValidation::computeNormalizedScores(valuedTree.topology(), std::span<const float>(invalidScores))); },
        "normalized saliency valuation must reject valuations that decrease toward the root");

    std::vector<double> nonFiniteScores = {10.0, 20.0, std::numeric_limits<double>::quiet_NaN()};
    requireThrows<std::invalid_argument>(
        [&]() { HierarchySaliencyMapValidation::validateHierarchyValuation(valuedTree.topology(), std::span<const double>(nonFiniteScores)); },
        "formal saliency valuation must reject non-finite values");

    std::vector<double> normalizedMaxScores = HierarchySaliencyMapValidation::computeNormalizedScores(valuedTree);
    requireNear(normalizedMaxScores[0], 0.0, 1e-12, "normalized max-tree altitude leaf score");
    requireNear(normalizedMaxScores[1], 0.2, 1e-12, "normalized max-tree altitude intermediate score");
    requireNear(normalizedMaxScores[2], 1.0, 1e-12, "normalized max-tree altitude root score");
    auto normalizedMaxMap = HierarchySaliencyMap::computeNormalizedAltitudeEdgeMap(valuedTree);
    requireNear(normalizedMaxMap.values[0], 1.0, 1e-12, "normalized max-tree saliency root edge");
    requireNear(normalizedMaxMap.values[1], 0.0, 1e-12, "normalized max-tree saliency same-smallest-node edge");

    auto maxPixelImage = HierarchySaliencyMapProjection::edgeMapToPixelImage(normalizedMaxMap, EdgeToPixelReducer::Max);
    requireNear((*maxPixelImage)[0], 1.0, 1e-12, "edge-to-pixel max source");
    requireNear((*maxPixelImage)[1], 1.0, 1e-12, "edge-to-pixel max shared endpoint");
    requireNear((*maxPixelImage)[2], 0.0, 1e-12, "edge-to-pixel max target");
    auto meanPixelImage = HierarchySaliencyMapProjection::edgeMapToPixelImage(normalizedMaxMap, EdgeToPixelReducer::Mean);
    requireNear((*meanPixelImage)[0], 1.0, 1e-12, "edge-to-pixel mean source");
    requireNear((*meanPixelImage)[1], 0.5, 1e-12, "edge-to-pixel mean shared endpoint");
    requireNear((*meanPixelImage)[2], 0.0, 1e-12, "edge-to-pixel mean target");

    auto highCut = HierarchySaliencyMapProjection::thresholdCut(normalizedMaxMap, 0.5);
    requireEqual(highCut.numRows, 1, "high-threshold contour rows");
    requireEqual(highCut.numColumns, 3, "high-threshold contour columns");
    requireEqual(highCut.size(), static_cast<std::size_t>(1), "high-threshold contour edge count");
    requireVectorEqual(highCut.sources, {0}, "high-threshold contour sources");
    requireVectorEqual(highCut.targets, {1}, "high-threshold contour targets");

    auto lowCut = HierarchySaliencyMapProjection::thresholdCut(normalizedMaxMap, 0.199);
    requireVectorEqual(lowCut.sources, {0}, "low-threshold contour sources");
    requireVectorEqual(lowCut.targets, {1}, "low-threshold contour targets");

    auto nodeContours = HierarchySaliencyMapProjection::nodeContourEdges(valuedTree);
    requireEqual(nodeContours.size(), static_cast<std::size_t>(1), "node contour edge count");
    requireVectorEqual(nodeContours.sources, {0}, "node contour sources");
    requireVectorEqual(nodeContours.targets, {1}, "node contour targets");
    requireVectorEqual(nodeContours.nodes, {2}, "node contour LCA nodes");

    auto incrementalContours = HierarchySaliencyMapProjection::computeIncrementalNodeContours(valuedTree);
    requireEqual(incrementalContours.numRows, 1, "incremental contour rows");
    requireEqual(incrementalContours.numColumns, 3, "incremental contour columns");
    requireEqual(incrementalContours.numNodeSlots, valuedTree.topology().numInternalNodeSlots(), "incremental contour node slots");
    requireVectorEqual(incrementalContours.offsets, {std::size_t{0}, std::size_t{0}, std::size_t{0}, std::size_t{1}}, "incremental contour offsets");
    requireVectorEqual(incrementalContours.sources, {0}, "incremental contour grouped sources");
    requireVectorEqual(incrementalContours.targets, {1}, "incremental contour grouped targets");

    auto projectedScoreMap = HierarchySaliencyMapProjection::projectNodeValuation(incrementalContours, std::span<const float>(nodeValuation));
    requireVectorEqual(projectedScoreMap.sources, {0}, "incremental projected score sources");
    requireVectorEqual(projectedScoreMap.targets, {1}, "incremental projected score targets");
    requireVectorEqual(projectedScoreMap.values, {30.0f}, "incremental projected score values");

    auto scoreCut = HierarchySaliencyMapProjection::thresholdByNodeValuation(incrementalContours, std::span<const float>(nodeValuation), 25.0f);
    requireVectorEqual(scoreCut.sources, {0}, "incremental score-cut sources");
    requireVectorEqual(scoreCut.targets, {1}, "incremental score-cut targets");

    auto minValuedTree = makeThreePixelMinTree();
    std::vector<double> normalizedMinScores = HierarchySaliencyMapValidation::computeNormalizedScores(minValuedTree);
    requireNear(normalizedMinScores[0], 0.0, 1e-12, "normalized min-tree altitude leaf score");
    requireNear(normalizedMinScores[1], 0.2, 1e-12, "normalized min-tree altitude intermediate score");
    requireNear(normalizedMinScores[2], 1.0, 1e-12, "normalized min-tree altitude root score");
    auto normalizedMinMap = HierarchySaliencyMap::computeNormalizedAltitudeEdgeMap(minValuedTree);
    requireNear(normalizedMinMap.values[0], 1.0, 1e-12, "normalized min-tree saliency root edge");
    requireNear(normalizedMinMap.values[1], 0.0, 1e-12, "normalized min-tree saliency same-smallest-node edge");

    auto explicitAdjacencyMap = HierarchySaliencyMap::computeTopologicalLevelEdgeMap(valuedTree.topology(), RegularGridAdjacency2D(1, 3, 1.5));
    requireVectorEqual(explicitAdjacencyMap.values, {1, 0}, "explicit adjacency topological saliency values");

    requireThrows<std::invalid_argument>(
        [&]() { static_cast<void>(HierarchySaliencyMap::computeTopologicalLevelEdgeMap(valuedTree.topology(), RegularGridAdjacency2D(2, 2, 1.0))); },
        "hierarchy saliency must reject adjacency domains that differ from the tree domain");

    auto noAdjacencyTree = makeThreePixelValuedTree(MorphologicalTreeKind::TreeOfShapes, std::nullopt);
    requireThrows<std::invalid_argument>([&]() { static_cast<void>(HierarchySaliencyMap::computeTopologicalLevelEdgeMap(noAdjacencyTree.topology())); },
                                         "hierarchy saliency must require explicit adjacency when the tree has none");
    std::vector<int> tosTopologicalValuation = HierarchySaliencyMap::computeTopologicalLevels(noAdjacencyTree.topology());
    auto tosFormalMap =
        HierarchySaliencyMap::computeSaliencyEdgeMap(noAdjacencyTree.topology(), RegularGridAdjacency2D(1, 3, 1.0),
                                                     std::span<const int>(tosTopologicalValuation), HierarchyValuationPolicy::RequireStrictHierarchy);
    requireVectorEqual(tosFormalMap.values, {1, 0}, "tree-of-shapes topological valuation saliency values");
    requireThrows<std::invalid_argument>(
        [&]() { static_cast<void>(HierarchySaliencyMap::computeNormalizedAltitudeEdgeMap(noAdjacencyTree, RegularGridAdjacency2D(1, 3, 1.0))); },
        "normalized altitude saliency must reject trees without max/min polarity");

    {
        MorphologicalTree stagedTree = valuedTree.topology().clone();
        auto editor = stagedTree.edit();
        editor.detach(0);

        requireThrows<std::invalid_argument>([&] { static_cast<void>(HierarchySaliencyMap::computeTopologicalLevels(stagedTree)); },
                                             "topological saliency levels must reject an active staged forest");
        requireThrows<std::invalid_argument>(
            [&] { HierarchySaliencyMapValidation::validateHierarchyValuation(stagedTree, std::span<const float>(nodeValuation)); },
            "saliency valuation validation must reject an active staged forest");
        requireThrows<std::invalid_argument>(
            [&] { static_cast<void>(HierarchySaliencyMap::computeEdgeMap(stagedTree, RegularGridAdjacency2D(1, 3, 1.0), [](NodeId node) { return node; })); },
            "LCA edge projection must reject an active staged forest");
        editor.rollback();
    }

    {
        auto equivalentDirectional = makeDirectionalThreePixelTree(true);
        const auto storedMap = HierarchySaliencyMap::computeTopologicalLevelEdgeMap(equivalentDirectional.topology());
        requireVectorEqual(storedMap.values, {1, 0}, "equivalent directional adjacency may define the stored saliency graph");

        auto distinctDirectional = makeDirectionalThreePixelTree(false);
        requireThrowsContaining<std::invalid_argument>(
            [&] { static_cast<void>(HierarchySaliencyMap::computeTopologicalLevelEdgeMap(distinctDirectional.topology())); }, "distinct minimum/maximum",
            "distinct directional adjacencies require an explicit saliency graph");
        const auto explicitMap = HierarchySaliencyMap::computeTopologicalLevelEdgeMap(distinctDirectional.topology(), RegularGridAdjacency2D(1, 3, 1.0));
        requireVectorEqual(explicitMap.values, {1, 0}, "explicit graph resolves directional adjacency ambiguity");
    }

    {
        const std::vector<NodeId> nodeParent{0, 0};
        const std::vector<NodeId> smallestNodeMap{1};
        const long double rootAltitude = 1.0L;
        const long double childAltitude = std::nextafter(rootAltitude, 2.0L);
        const std::vector<long double> altitude{rootAltitude, childAltitude};
        auto preciseValuedTree = MorphologicalTreeFactory::createFromNativeTopology(std::span<const NodeId>(nodeParent), std::span<const NodeId>(smallestNodeMap),
                                                                                  std::span<const long double>(altitude), 0,
                                                                                  makeMorphologicalTreeSemantics(MorphologicalTreeKind::MaxTree));

        const std::vector<double> preciseNormalized = HierarchySaliencyMapValidation::computeNormalizedScores(preciseValuedTree);
        requireNear(preciseNormalized[0], 1.0, 0.0, "valuedTree normalization must preserve a long-double root level");
        requireNear(preciseNormalized[1], 0.0, 0.0, "valuedTree normalization must preserve a distinct long-double child level");
    }

    return 0;
}
