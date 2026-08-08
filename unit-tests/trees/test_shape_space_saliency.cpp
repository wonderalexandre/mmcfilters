#include "support/TestSupport.hpp"

#include "mmcfilters/trees/saliency/HierarchySaliencyMap.hpp"
#include "mmcfilters/trees/saliency/ShapeSpaceSaliency.hpp"

#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

namespace {

WeightedMorphologicalTree<std::uint8_t> makeNestedChain(MorphologicalTreeKind kind = MorphologicalTreeKind::MAX_TREE,
                                                        std::optional<RegularGridAdjacency2D> adjacency = RegularGridAdjacency2D(1, 4, 1.0)) {
    // Proper-part supports form B subset A subset Omega:
    // B={0}, A={0,1}, Omega={0,1,2,3}.
    // The last two neighbouring proper parts have the same direct owner.
    const std::vector<NodeId> parent = {4, 5, 6, 6, 5, 6, 6};
    const std::vector<std::uint8_t> altitude = {3, 2, 1, 1, 3, 2, 1};
    return MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(parent), std::span<const std::uint8_t>(altitude), 1, 4, kind,
                                                           std::move(adjacency));
}

WeightedMorphologicalTree<std::uint8_t> makeTiedBranchingTree() {
    // Two sibling singleton regions separated by a root-owned proper part.
    const std::vector<NodeId> parent = {3, 5, 4, 5, 5, 5};
    const std::vector<std::uint8_t> altitude = {2, 1, 2, 2, 2, 1};
    return MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(parent), std::span<const std::uint8_t>(altitude), 1, 3,
                                                           MorphologicalTreeKind::MAX_TREE, RegularGridAdjacency2D(1, 3, 1.0));
}

WeightedMorphologicalTree<std::uint8_t> makeMultiStageTieTree() {
    // Shape-space topology:
    //
    //          4
    //        /   \
    //       2     3
    //     /   \
    //    0     1
    //
    // Proper parts are ordered so every represented region is connected in
    // the 1-D image domain: owner sequence [0, 2, 1, 4, 3].
    const std::vector<NodeId> parent = {5, 7, 6, 9, 8, 7, 7, 9, 9, 9};
    const std::vector<std::uint8_t> altitude = {3, 2, 3, 1, 2, 3, 3, 2, 2, 1};
    return MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(parent), std::span<const std::uint8_t>(altitude), 1, 5,
                                                           MorphologicalTreeKind::MAX_TREE, RegularGridAdjacency2D(1, 5, 1.0));
}

void verifyNestedChainPipeline() {
    auto weighted = makeNestedChain();
    const MorphologicalTree& tree = weighted.topology();
    const std::vector<float> attribute = {1.0f, 4.0f, 10.0f};

    const auto extinction = ShapeSpaceSaliency::computeExtinctionValues(tree, std::span<const float>(attribute), ShapeSpaceExtremaPolarity::Minima);
    requireEqual(extinction.extrema.size(), std::size_t{1}, "nested-chain minimum count");
    requireVectorEqual(extinction.nodeScores, {9.0f, 0.0f, 0.0f}, "nested-chain extinction scores");
    requireEqual(extinction.extrema[0].representative, NodeId{0}, "nested-chain minimum representative");
    requireEqual(extinction.extrema[0].birthLevel, 1.0f, "nested-chain minimum birth level");
    requireEqual(extinction.extrema[0].deathLevel, 10.0f, "nested-chain minimum death level");
    requireEqual(extinction.extrema[0].extinction, 9.0f, "nested-chain minimum extinction");

    const auto explicitProjection =
        ShapeSpaceSaliency::projectContourScores(tree, std::span<const float>(extinction.nodeScores), RegularGridAdjacency2D(1, 4, 1.0));
    requireVectorEqual(explicitProjection.sources, {0, 1, 2}, "nested-chain contour sources");
    requireVectorEqual(explicitProjection.targets, {1, 2, 3}, "nested-chain contour targets");
    requireVectorEqual(explicitProjection.values, {9.0f, 0.0f, 0.0f}, "nested-chain contour scores");

    const auto storedProjection = ShapeSpaceSaliency::projectContourScores(tree, std::span<const float>(extinction.nodeScores));
    requireVectorEqual(storedProjection.sources, explicitProjection.sources, "stored-adjacency contour sources");
    requireVectorEqual(storedProjection.targets, explicitProjection.targets, "stored-adjacency contour targets");
    requireVectorEqual(storedProjection.values, explicitProjection.values, "stored-adjacency contour scores");

    const auto explicitResult =
        ShapeSpaceSaliency::compute(tree, std::span<const float>(attribute), ShapeSpaceExtremaPolarity::Minima, RegularGridAdjacency2D(1, 4, 1.0));
    requireVectorEqual(explicitResult.nodeScores, extinction.nodeScores, "one-shot explicit extinction scores");
    requireVectorEqual(explicitResult.edgeMap.values, explicitProjection.values, "one-shot explicit contour scores");
    requireEqual(explicitResult.extrema.size(), extinction.extrema.size(), "one-shot explicit extrema count");

    const auto storedResult = ShapeSpaceSaliency::compute(tree, std::span<const float>(attribute), ShapeSpaceExtremaPolarity::Minima);
    requireVectorEqual(storedResult.nodeScores, extinction.nodeScores, "one-shot stored extinction scores");
    requireVectorEqual(storedResult.edgeMap.values, explicitProjection.values, "one-shot stored contour scores");

    // LCA projection answers a different question. Extending the one minimum's
    // score monotonically to its ancestors marks both hierarchy transitions,
    // while full-contour accumulation marks only the boundary of B.
    const std::vector<float> lcaValuation = {9.0f, 9.0f, 9.0f};
    const auto lcaMap = HierarchySaliencyMap::computeSaliencyEdgeMap(tree, std::span<const float>(lcaValuation));
    requireVectorEqual(lcaMap.values, {9.0f, 9.0f, 0.0f}, "LCA projection differs from shaping contour accumulation");
    require(lcaMap.values != storedResult.edgeMap.values, "shaping saliency must not collapse to LCA projection");

    require(std::isfinite(extinction.nodeScores[0]), "dominant extinction must be finite");
    requireEqual(explicitProjection.values[2], 0.0f, "same-owner proper-part edge must have zero saliency");
}

void verifyMinimaAndMaximaPolarity() {
    auto weighted = makeNestedChain();
    const MorphologicalTree& tree = weighted.topology();
    const std::vector<double> attribute = {1.0, 5.0, 2.0};

    const auto minima = ShapeSpaceSaliency::compute(tree, std::span<const double>(attribute), ShapeSpaceExtremaPolarity::Minima);
    requireVectorEqual(minima.nodeScores, {4.0, 0.0, 3.0}, "minimum-oriented extinction scores");
    requireVectorEqual(minima.edgeMap.values, {4.0, 0.0, 0.0}, "minimum-oriented contour scores");
    requireEqual(minima.extrema.size(), std::size_t{2}, "minimum-oriented extrema count");
    requireEqual(minima.extrema[0].representative, NodeId{0}, "first minimum representative");
    requireEqual(minima.extrema[0].birthLevel, 1.0, "first minimum birth level");
    requireEqual(minima.extrema[0].deathLevel, 5.0, "first minimum death level");
    requireEqual(minima.extrema[0].extinction, 4.0, "first minimum extinction");
    requireEqual(minima.extrema[1].representative, NodeId{2}, "second minimum representative");
    requireEqual(minima.extrema[1].birthLevel, 2.0, "second minimum birth level");
    requireEqual(minima.extrema[1].deathLevel, 5.0, "second minimum death level");
    requireEqual(minima.extrema[1].extinction, 3.0, "second minimum extinction");

    const auto maxima = ShapeSpaceSaliency::compute(tree, std::span<const double>(attribute), ShapeSpaceExtremaPolarity::Maxima);
    requireVectorEqual(maxima.nodeScores, {0.0, 4.0, 0.0}, "maximum-oriented extinction scores");
    requireVectorEqual(maxima.edgeMap.values, {0.0, 4.0, 0.0}, "maximum-oriented contour scores");
    requireEqual(maxima.extrema.size(), std::size_t{1}, "maximum-oriented extrema count");
    requireEqual(maxima.extrema[0].representative, NodeId{1}, "maximum representative");
    requireEqual(maxima.extrema[0].birthLevel, 5.0, "maximum birth level stays in original domain");
    requireEqual(maxima.extrema[0].deathLevel, 1.0, "maximum death level stays in original domain");
    requireEqual(maxima.extrema[0].extinction, 4.0, "maximum extinction is birth minus death");
    require(std::isfinite(maxima.nodeScores[1]), "dominant maximum extinction must be finite");
}

void verifyPlateauxAndTiesAreDeterministic() {
    auto chain = makeNestedChain();
    const std::vector<float> plateauAttribute = {1.0f, 1.0f, 5.0f};
    const auto plateauFirst =
        ShapeSpaceSaliency::computeExtinctionValues(chain.topology(), std::span<const float>(plateauAttribute), ShapeSpaceExtremaPolarity::Minima);
    const auto plateauSecond =
        ShapeSpaceSaliency::computeExtinctionValues(chain.topology(), std::span<const float>(plateauAttribute), ShapeSpaceExtremaPolarity::Minima);
    requireEqual(plateauFirst.extrema.size(), std::size_t{1}, "flat minimum must be one extremum");
    requireVectorEqual(plateauFirst.nodeScores, {0.0f, 4.0f, 0.0f}, "flat minimum canonical representative");
    requireEqual(plateauFirst.extrema[0].representative, NodeId{1}, "flat minimum uses the outer plateau representative");
    requireVectorEqual(plateauSecond.nodeScores, plateauFirst.nodeScores, "flat minimum result must be deterministic");

    auto branching = makeTiedBranchingTree();
    const std::vector<float> tiedAttribute = {1.0f, 1.0f, 5.0f};
    const auto tiedFirst =
        ShapeSpaceSaliency::computeExtinctionValues(branching.topology(), std::span<const float>(tiedAttribute), ShapeSpaceExtremaPolarity::Minima);
    const auto tiedSecond =
        ShapeSpaceSaliency::computeExtinctionValues(branching.topology(), std::span<const float>(tiedAttribute), ShapeSpaceExtremaPolarity::Minima);
    requireEqual(tiedFirst.extrema.size(), std::size_t{2}, "tied sibling minima count");
    requireVectorEqual(tiedFirst.nodeScores, {4.0f, 4.0f, 0.0f}, "tied sibling minima scores");
    requireEqual(tiedFirst.extrema[0].representative, NodeId{0}, "first tied minimum representative");
    requireEqual(tiedFirst.extrema[1].representative, NodeId{1}, "second tied minimum representative");
    requireVectorEqual(tiedSecond.nodeScores, tiedFirst.nodeScores, "tied sibling result must be deterministic");
}

void verifyMultiStageBranchMerges() {
    auto weighted = makeMultiStageTieTree();
    const MorphologicalTree& tree = weighted.topology();

    // Nodes 0 and 1 are tied minima born at level 1. They first meet through
    // node 2 at level 5, so node 1 dies there and node 0 wins by NodeId. The
    // surviving node 0 then meets the stronger minimum 3 at level 10.
    const std::vector<double> attribute = {1.0, 1.0, 5.0, 0.0, 10.0};
    const auto result = ShapeSpaceSaliency::compute(tree, std::span<const double>(attribute), ShapeSpaceExtremaPolarity::Minima);

    requireEqual(result.extrema.size(), std::size_t{3}, "multi-stage extrema count");
    requireVectorEqual(result.nodeScores, {9.0, 4.0, 0.0, 10.0, 0.0}, "multi-stage extinction scores");

    requireEqual(result.extrema[0].representative, NodeId{0}, "multi-stage tie winner representative");
    requireEqual(result.extrema[0].birthLevel, 1.0, "multi-stage tie winner birth level");
    requireEqual(result.extrema[0].deathLevel, 10.0, "multi-stage tie winner survives to later merge");
    requireEqual(result.extrema[0].extinction, 9.0, "multi-stage tie winner extinction");

    requireEqual(result.extrema[1].representative, NodeId{1}, "multi-stage tied loser representative");
    requireEqual(result.extrema[1].birthLevel, 1.0, "multi-stage tied loser birth level");
    requireEqual(result.extrema[1].deathLevel, 5.0, "multi-stage tied loser dies at first merge");
    requireEqual(result.extrema[1].extinction, 4.0, "multi-stage tied loser extinction");

    requireEqual(result.extrema[2].representative, NodeId{3}, "multi-stage dominant representative");
    requireEqual(result.extrema[2].birthLevel, 0.0, "multi-stage dominant birth level");
    requireEqual(result.extrema[2].deathLevel, 10.0, "multi-stage dominant finite death level");
    requireEqual(result.extrema[2].extinction, 10.0, "multi-stage dominant extinction");

    requireVectorEqual(result.edgeMap.sources, {0, 1, 2, 3}, "multi-stage contour sources");
    requireVectorEqual(result.edgeMap.targets, {1, 2, 3, 4}, "multi-stage contour targets");
    requireVectorEqual(result.edgeMap.values, {9.0, 4.0, 4.0, 10.0}, "multi-stage contour scores");
}

void verifyValidation() {
    auto weighted = makeNestedChain();
    const MorphologicalTree& tree = weighted.topology();
    const std::vector<float> attribute = {1.0f, 4.0f, 10.0f};
    const std::vector<float> shortAttribute = {1.0f, 4.0f};
    requireThrows<std::invalid_argument>(
        [&]() {
            static_cast<void>(ShapeSpaceSaliency::computeExtinctionValues(tree, std::span<const float>(shortAttribute), ShapeSpaceExtremaPolarity::Minima));
        },
        "shape-space extinction must reject a short attribute buffer");

    std::vector<float> nanAttribute = attribute;
    nanAttribute[1] = std::numeric_limits<float>::quiet_NaN();
    requireThrows<std::invalid_argument>(
        [&]() {
            static_cast<void>(ShapeSpaceSaliency::computeExtinctionValues(tree, std::span<const float>(nanAttribute), ShapeSpaceExtremaPolarity::Minima));
        },
        "shape-space extinction must reject NaN attributes");

    const std::vector<float> scores = {9.0f, 0.0f, 0.0f};
    requireThrows<std::invalid_argument>(
        [&]() { static_cast<void>(ShapeSpaceSaliency::projectContourScores(tree, std::span<const float>(scores), RegularGridAdjacency2D(1, 3, 1.0))); },
        "shape-space contour projection must reject an adjacency-domain mismatch");

    std::vector<float> negativeScores = scores;
    negativeScores[1] = -1.0f;
    requireThrows<std::invalid_argument>([&]() { static_cast<void>(ShapeSpaceSaliency::projectContourScores(tree, std::span<const float>(negativeScores))); },
                                         "shape-space contour projection must reject negative node scores");

    std::vector<float> nanScores = scores;
    nanScores[1] = std::numeric_limits<float>::quiet_NaN();
    requireThrows<std::invalid_argument>([&]() { static_cast<void>(ShapeSpaceSaliency::projectContourScores(tree, std::span<const float>(nanScores))); },
                                         "shape-space contour projection must reject NaN node scores");

    auto noAdjacency = makeNestedChain(MorphologicalTreeKind::TREE_OF_SHAPES, std::nullopt);
    const auto explicitTreeOfShapes = ShapeSpaceSaliency::compute(noAdjacency.topology(), std::span<const float>(attribute), ShapeSpaceExtremaPolarity::Minima,
                                                                  RegularGridAdjacency2D(1, 4, 1.0));
    requireVectorEqual(explicitTreeOfShapes.edgeMap.values, {9.0f, 0.0f, 0.0f}, "shape-space saliency accepts a tree of shapes with explicit adjacency");

    requireThrows<std::invalid_argument>(
        [&]() { static_cast<void>(ShapeSpaceSaliency::compute(noAdjacency.topology(), std::span<const float>(attribute), ShapeSpaceExtremaPolarity::Minima)); },
        "shape-space one-shot stored overload must require stored adjacency");

    const double maxFinite = std::numeric_limits<double>::max();
    const std::vector<double> overflowingAttribute = {-maxFinite, 0.0, maxFinite};
    requireThrows<std::overflow_error>(
        [&]() {
            static_cast<void>(
                ShapeSpaceSaliency::computeExtinctionValues(tree, std::span<const double>(overflowingAttribute), ShapeSpaceExtremaPolarity::Minima));
        },
        "shape-space extinction must reject a non-finite level difference");

    {
        auto editor = weighted.edit();
        requireThrows<std::invalid_argument>(
            [&]() {
                static_cast<void>(ShapeSpaceSaliency::computeExtinctionValues(tree, std::span<const float>(attribute), ShapeSpaceExtremaPolarity::Minima));
            },
            "shape-space extinction must reject an active edit session");
        editor.commit();
    }
}

void verifyDeadNodeSlots() {
    auto weighted = makeNestedChain();
    weighted.mergeNodeIntoParent(1);
    const MorphologicalTree& tree = weighted.topology();
    require(!tree.isAlive(1), "dead-slot fixture must remove the middle node");

    const std::vector<float> attribute = {1.0f, 0.0f, 10.0f};
    const auto result = ShapeSpaceSaliency::compute(tree, std::span<const float>(attribute), ShapeSpaceExtremaPolarity::Minima);
    requireVectorEqual(result.nodeScores, {9.0f, 0.0f, 0.0f}, "dead slots remain zero in node scores");
    requireVectorEqual(result.edgeMap.values, {9.0f, 0.0f, 0.0f}, "dead slots are skipped by contour projection");
}

} // namespace

int main() {
    verifyNestedChainPipeline();
    verifyMinimaAndMaximaPolarity();
    verifyPlateauxAndTiesAreDeterministic();
    verifyMultiStageBranchMerges();
    verifyValidation();
    verifyDeadNodeSlots();
    return 0;
}
