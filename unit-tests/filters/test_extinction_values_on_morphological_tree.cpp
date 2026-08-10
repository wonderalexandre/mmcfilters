#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/AttributeComputation.hpp"
#include "mmcfilters/filters/ExtinctionValues.hpp"
#include "mmcfilters/trees/saliency/HierarchySaliencyMapValidation.hpp"
#include "mmcfilters/trees/saliency/HierarchySaliencyMap.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <stdexcept>
#include <vector>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

namespace {

WeightedMorphologicalTree<std::uint8_t> makeFlatThreeLeafMaxTreeFixture() {
    const std::vector<NodeId> parent{3, 4, 5, 6, 6, 6, 6};
    const std::vector<std::uint8_t> altitude{5, 5, 5, 5, 5, 5, 0};
    return MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(parent), std::span<const std::uint8_t>(altitude), 1, 3,
                                                           MorphologicalTreeKind::MAX_TREE, RegularGridAdjacency2D(1, 3, 1.5));
}

template <class Real> std::vector<Real> expectedExtinctionValueAttribute(const MorphologicalTree& tree, const std::vector<RegionalExtremaNode<Real>>& records) {
    std::vector<Real> expected(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), Real{0});
    for (const RegionalExtremaNode<Real>& record : records) {
        Real& value = expected[static_cast<std::size_t>(record.leaf)];
        if (value < record.extinction) {
            value = record.extinction;
        }
    }

    for (NodeId nodeId : tree.getPostOrderNodes()) {
        Real& nodeValue = expected[static_cast<std::size_t>(nodeId)];
        for (NodeId childId : tree.getChildren(nodeId)) {
            const Real childValue = expected[static_cast<std::size_t>(childId)];
            if (nodeValue < childValue) {
                nodeValue = childValue;
            }
        }
    }
    return expected;
}

void verifyDeterministicTieOrdering() {
    auto weighted = makeFlatThreeLeafMaxTreeFixture();
    const MorphologicalTree& tree = weighted.topology();
    std::vector<float> attr(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), 1.0f);

    ExtinctionValues<std::uint8_t> extinction(weighted, attr);
    const auto& records = extinction.getRegionalExtrema();

    requireEqual(static_cast<int>(records.size()), 3, "flat extinction record count");
    requireEqual(records[0].extinction, std::numeric_limits<float>::max(), "dominant extinction sentinel");
    requireEqual(records[1].extinction, 1.0f, "first tied loser extinction");
    requireEqual(records[2].extinction, 1.0f, "second tied loser extinction");
    require(records[1].leaf < records[2].leaf, "equal extinction records must be sorted deterministically by leaf NodeId");

    std::vector<float> valuation = extinction.getExtinctionValueAttribute();
    requireVectorEqual(valuation, expectedExtinctionValueAttribute(tree, records), "extinction value attribute must extend leaf extinctions to ancestors");
    for (const RegionalExtremaNode<float>& record : records) {
        requireEqual(valuation[static_cast<std::size_t>(record.leaf)], record.extinction, "each leaf extremum must receive its own extinction value");
    }
    HierarchySaliencyMapValidation::validateHierarchyValuation(tree, std::span<const float>(valuation), HierarchyValuationPolicy::AllowLevelCollapse,
                                                               HierarchyValuationRangePolicy::RequireNonNegative);
    requireEqual(valuation[static_cast<std::size_t>(tree.getRoot())], std::numeric_limits<float>::max(),
                 "dominant extinction sentinel becomes the root hierarchy valuation");

    std::vector<int> rankedValuation = extinction.computeRankedExtinctionValueAttribute();
    HierarchySaliencyMapValidation::validateHierarchyValuation(tree, std::span<const int>(rankedValuation), HierarchyValuationPolicy::AllowLevelCollapse,
                                                               HierarchyValuationRangePolicy::RequireNonNegative);
    requireEqual(rankedValuation[static_cast<std::size_t>(tree.getRoot())], 1, "ranked extinction valuation preserves the dominant root level");

    auto hierarchicalWatershed = extinction.computeFormalSaliencyEdgeMap();
    requireVectorEqual(hierarchicalWatershed.values, {1.0f, 1.0f}, "Cousty persistence uses the losing child extinction at each MST merge");
    auto monotoneProjection = extinction.computeMonotoneExtinctionProjection();
    requireVectorEqual(monotoneProjection.values, {std::numeric_limits<float>::max(), std::numeric_limits<float>::max()},
                       "legacy monotone projection remains explicitly available");
}

void verifyHandComputedPersistencePath() {
    auto image = makeImage(1, 5, {0, 2, 1, 2, 0});
    auto weighted = MorphologicalTreeFactory::createMinTree(image, RegularGridAdjacency2D(1, 5, 1.0));
    const MorphologicalTree& tree = weighted.topology();
    std::vector<float> leafExtinction(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), 0.0f);
    requireEqual(static_cast<int>(tree.getLeaves().size()), 3, "three-minimum persistence fixture leaf count");
    for (NodeId leaf : tree.getLeaves()) {
        NodeId representative = InvalidNode;
        for (NodeId properPart : tree.getProperParts(leaf)) {
            representative = representative == InvalidNode ? properPart : std::min(representative, properPart);
        }
        require(representative != InvalidNode, "each regional minimum owns a pixel");
        if (representative == 0) {
            leafExtinction[static_cast<std::size_t>(leaf)] = 9.0f;
        } else if (representative == 2) {
            leafExtinction[static_cast<std::size_t>(leaf)] = 5.0f;
        } else if (representative == 4) {
            leafExtinction[static_cast<std::size_t>(leaf)] = 2.0f;
        } else {
            throw std::runtime_error("unexpected regional-minimum representative");
        }
    }

    auto saliency = HierarchicalWatershedSaliency::compute(weighted.asView(), std::span<const float>(leafExtinction), RegularGridAdjacency2D(1, 5, 1.0));
    requireVectorEqual(saliency.sources, {0, 1, 2, 3}, "persistence path sources");
    requireVectorEqual(saliency.targets, {1, 2, 3, 4}, "persistence path targets");
    requireVectorEqual(saliency.values, {0.0f, 5.0f, 0.0f, 2.0f}, "hand-computed min-child persistence map");

    auto ranked = HierarchicalWatershedSaliency::computeRanked(weighted.asView(), std::span<const float>(leafExtinction), RegularGridAdjacency2D(1, 5, 1.0));
    requireVectorEqual(ranked.values, {0, 2, 0, 1}, "canonical persistence edge ranks");
}

void verifyDominantExtremumDoesNotBecomeAnEdgeLevel() {
    auto image = makeImage(1, 3, {0, 1, 2});
    auto weighted = MorphologicalTreeFactory::createMinTree(image, RegularGridAdjacency2D(1, 3, 1.0));
    const MorphologicalTree& tree = weighted.topology();
    requireEqual(static_cast<int>(tree.getLeaves().size()), 1, "single-basin fixture leaf count");

    std::vector<float> attribute(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), 1.0f);
    ExtinctionValues<std::uint8_t> extinction(weighted, attribute);
    requireEqual(extinction.getRegionalExtrema().front().extinction, std::numeric_limits<float>::max(), "single basin has the dominant extinction sentinel");

    const auto saliency = extinction.computeFormalSaliencyEdgeMap();
    requireVectorEqual(saliency.values, {0.0f, 0.0f}, "dominant sentinel is not a hierarchical-watershed edge level");
}

void verifyFormalAttributeRejectsInvalidExtinctionValues() {
    auto weighted = makeFlatThreeLeafMaxTreeFixture();
    const MorphologicalTree& tree = weighted.topology();
    std::vector<float> attr(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), -1.0f);

    ExtinctionValues<std::uint8_t> extinction(weighted, attr);
    requireThrows<std::invalid_argument>([&]() { static_cast<void>(extinction.getExtinctionValueAttribute()); },
                                         "formal extinction value attribute must reject negative extinction records");
    requireThrows<std::invalid_argument>([&]() { static_cast<void>(extinction.computeFormalSaliencyEdgeMap()); },
                                         "formal extinction saliency must reject negative extinction records");
}

void verifyTreeOfShapesIsRejected() {
    auto image = makeComponentTreeFixture();
    auto tos = makeWeightedTreeOfShapes(image, ToSInterpolation::SelfDual);
    std::vector<float> tosAttr(static_cast<std::size_t>(tos->topology().getNumInternalNodeSlots()), 1.0f);
    requireThrows<std::invalid_argument>([&]() { ExtinctionValues<std::uint8_t> invalid(*tos, tosAttr); }, "ExtinctionValues must reject tree of shapes");
}

} // namespace

int main() {
    verifyHandComputedPersistencePath();
    verifyDominantExtremumDoesNotBecomeAnEdgeLevel();
    auto image = makeComponentTreeFixture();

    for (bool isMaxtree : {true, false}) {
        auto weighted = makeWeightedComponentTree(image, isMaxtree);
        auto [attrNames, attrBuffer] = AttributeComputation::computeSingleAttribute(*weighted, LEVEL);
        (void)attrNames;

        if (isMaxtree && contract::validationsEnabled) {
            requireThrows<std::invalid_argument>([&]() { ExtinctionValues<std::uint8_t> invalidExtinction(*weighted, std::vector<float>{1.0f}); },
                                                 "ExtinctionValues<std::uint8_t> must reject short vector attribute buffer");
        }

        ExtinctionValues<std::uint8_t> extinction(*weighted, attrBuffer);
        const auto keepAll = ExtinctionSelectionPolicy<float>::byTopK(1024);
        if constexpr (contract::validationsEnabled) {
            requireThrows<std::invalid_argument>([&]() { static_cast<void>(extinction.filtering(ExtinctionSelectionPolicy<float>::byTopK(-1))); },
                                                 isMaxtree ? "weighted max-tree extinction filtering must reject negative keep count"
                                                           : "weighted min-tree extinction filtering must reject negative keep count");
            requireThrows<std::invalid_argument>(
                [&]() { static_cast<void>(extinction.contourMap(ExtinctionSelectionPolicy<float>::byTopK(-1), ExtinctionContourScorePolicy::RankScore)); },
                isMaxtree ? "weighted max-tree extinction saliency must reject negative keep count"
                          : "weighted min-tree extinction saliency must reject negative keep count");
            requireThrows<std::invalid_argument>(
                [&]() { static_cast<void>(extinction.filtering(ExtinctionSelectionPolicy<float>::byThreshold(std::numeric_limits<float>::quiet_NaN()))); },
                isMaxtree ? "weighted max-tree extinction threshold filtering must reject NaN"
                          : "weighted min-tree extinction threshold filtering must reject NaN");
        }

        std::vector<float> hierarchyValuation = extinction.getExtinctionValueAttribute();
        requireVectorEqual(hierarchyValuation, expectedExtinctionValueAttribute(weighted->topology(), extinction.getRegionalExtrema()),
                           isMaxtree ? "weighted max-tree extinction attribute extends leaf records"
                                     : "weighted min-tree extinction attribute extends leaf records");
        HierarchySaliencyMapValidation::validateHierarchyValuation(weighted->topology(), std::span<const float>(hierarchyValuation),
                                                                   HierarchyValuationPolicy::AllowLevelCollapse,
                                                                   HierarchyValuationRangePolicy::RequireNonNegative);
        auto formalExtinctionMap = extinction.computeFormalSaliencyEdgeMap();
        auto directFormalExtinctionMap = HierarchySaliencyMap::computeSaliencyEdgeMap(weighted->topology(), std::span<const float>(hierarchyValuation),
                                                                                      HierarchyValuationPolicy::AllowLevelCollapse);
        auto monotoneProjection = extinction.computeMonotoneExtinctionProjection();
        requireVectorEqual(monotoneProjection.values, directFormalExtinctionMap.values,
                           isMaxtree ? "weighted max-tree preserves explicit monotone extinction projection"
                                     : "weighted min-tree preserves explicit monotone extinction projection");
        for (float value : formalExtinctionMap.values) {
            require(std::isfinite(value) && value >= 0.0f, isMaxtree ? "weighted max-tree persistence saliency is finite and non-negative"
                                                                     : "weighted min-tree persistence saliency is finite and non-negative");
        }

        std::vector<int> rankedHierarchyValuation = extinction.computeRankedExtinctionValueAttribute();
        HierarchySaliencyMapValidation::validateHierarchyValuation(weighted->topology(), std::span<const int>(rankedHierarchyValuation),
                                                                   HierarchyValuationPolicy::AllowLevelCollapse,
                                                                   HierarchyValuationRangePolicy::RequireNonNegative);
        auto rankedFormalExtinctionMap = extinction.computeRankedFormalSaliencyEdgeMap();
        auto expectedRankedFormalExtinctionMap = HierarchySaliencyMap::rankEdgeSaliencyMap(formalExtinctionMap);
        requireVectorEqual(rankedFormalExtinctionMap.values, expectedRankedFormalExtinctionMap.values,
                           isMaxtree ? "weighted max-tree ranked persistence saliency uses effective edge ranks"
                                     : "weighted min-tree ranked persistence saliency uses effective edge ranks");
        auto rankedMonotoneProjection = extinction.computeRankedMonotoneExtinctionProjection();
        auto directRankedMonotoneProjection = HierarchySaliencyMap::computeCanonicalRankedSaliencyEdgeMap(
            weighted->topology(), std::span<const float>(hierarchyValuation), HierarchyValuationPolicy::AllowLevelCollapse);
        requireVectorEqual(rankedMonotoneProjection.values, directRankedMonotoneProjection.values,
                           isMaxtree ? "weighted max-tree ranked monotone projection uses canonical edge ranks"
                                     : "weighted min-tree ranked monotone projection uses canonical edge ranks");

        auto filtered = extinction.filtering(keepAll);
        auto reconstructed = weighted->reconstructionImage();
        const auto& records = extinction.getRegionalExtrema();
        if (!records.empty()) {
            auto thresholdKeepAll = extinction.filtering(ExtinctionSelectionPolicy<float>::byThreshold(-1.0f));
            requireVectorEqual(collectImageValues(thresholdKeepAll), collectImageValues(filtered),
                               isMaxtree ? "weighted max-tree extinction threshold below all values keeps all extrema"
                                         : "weighted min-tree extinction threshold below all values keeps all extrema");

            auto thresholdDominant = extinction.filtering(ExtinctionSelectionPolicy<float>::byThreshold(records.front().extinction));
            auto strongestOnly = extinction.filtering(ExtinctionSelectionPolicy<float>::byTopK(1));
            requireVectorEqual(collectImageValues(thresholdDominant), collectImageValues(strongestOnly),
                               isMaxtree ? "weighted max-tree extinction threshold at dominant value keeps strongest extremum"
                                         : "weighted min-tree extinction threshold at dominant value keeps strongest extremum");
        }

        requireVectorEqual(collectImageValues(filtered), collectImageValues(reconstructed),
                           isMaxtree ? "weighted max-tree extinction keep-all reconstruction" : "weighted min-tree extinction keep-all reconstruction");

        if (isMaxtree && contract::validationsEnabled) {
            auto staleWeighted = makeWeightedComponentTree(image, true);
            auto [staleNames, staleAttr] = AttributeComputation::computeSingleAttribute(*staleWeighted, LEVEL);
            (void)staleNames;
            ExtinctionValues<std::uint8_t> staleExtinction(*staleWeighted, staleAttr);
            staleWeighted->mergeNodeIntoParent(4);
            requireThrows<std::logic_error>([&]() { static_cast<void>(staleExtinction.filtering(keepAll)); },
                                            "ExtinctionValues<std::uint8_t> filtering must reject use after topology mutation");
            requireThrows<std::logic_error>([&]() { static_cast<void>(staleExtinction.contourMap(keepAll, ExtinctionContourScorePolicy::RankScore)); },
                                            "ExtinctionValues<std::uint8_t> saliency must reject use after topology mutation");
            requireThrows<std::logic_error>([&]() { static_cast<void>(staleExtinction.getExtinctionValueAttribute()); },
                                            "ExtinctionValues<std::uint8_t> extinction value attribute must reject use after topology mutation");
            requireThrows<std::logic_error>([&]() { static_cast<void>(staleExtinction.computeFormalSaliencyEdgeMap()); },
                                            "ExtinctionValues<std::uint8_t> formal saliency must reject use after topology mutation");
            requireThrows<std::logic_error>([&]() { static_cast<void>(staleExtinction.getRegionalExtrema()); },
                                            "ExtinctionValues<std::uint8_t> tuples must reject use after topology mutation");
        }

        const auto [invalidParent, exportedAltitude] = weighted->exportHigraHierarchy();
        std::vector<int> validHigraAltitude(exportedAltitude.begin(), exportedAltitude.end());
        auto invalidWeighted = MorphologicalTreeFactory::createFromHigraParent<int>(
            std::span<const NodeId>(invalidParent), std::span<const int>(validHigraAltitude), image->getNumRows(), image->getNumCols(),
            isMaxtree ? MorphologicalTreeKind::MAX_TREE : MorphologicalTreeKind::MIN_TREE,
            RegularGridAdjacency2D(image->getNumRows(), image->getNumCols(), 1.5));
        auto [invalidLevelNames, invalidLevelBuffer] = AttributeComputation::computeSingleAttribute(invalidWeighted.asView(), LEVEL);
        (void)invalidLevelNames;
        ExtinctionValues<int> invalidExtinction(invalidWeighted, invalidLevelBuffer);
        auto invalidSaliencyBefore = collectImageValues(invalidExtinction.contourMap(keepAll, ExtinctionContourScorePolicy::ExtinctionValue));

        std::vector<int> invalidAltitude = invalidWeighted.getAltitudeBuffer();
        const int invalidOffset = isMaxtree ? 300 : -300;
        for (int& level : invalidAltitude) {
            level += invalidOffset;
        }
        invalidWeighted.setAltitudeBuffer(std::move(invalidAltitude));
        requireVectorEqual(collectImageValues(invalidExtinction.contourMap(keepAll, ExtinctionContourScorePolicy::ExtinctionValue)), invalidSaliencyBefore,
                           isMaxtree ? "weighted max-tree extinction saliency ignores replaced altitude"
                                     : "weighted min-tree extinction saliency ignores replaced altitude");
        requireVectorEqual(collectImageValues(invalidExtinction.filtering(keepAll)), collectImageValues(invalidWeighted.reconstructionImage()),
                           isMaxtree ? "ExtinctionValues<std::uint8_t> object must preserve altitude above 255"
                                     : "ExtinctionValues<std::uint8_t> object must preserve negative altitude");
    }

    {
        auto weighted = makeWeightedComponentTree(image, true);
        weighted->mergeNodeIntoParent(4);

        auto expectedAfterMerge = ImageUInt8::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        for (int pixel = 0; pixel < expectedAfterMerge->getSize(); ++pixel) {
            const NodeId nodeId = weighted->topology().getProperPartOwner(pixel);
            (*expectedAfterMerge)[pixel] = static_cast<uint8_t>(weighted->getAltitude(nodeId));
        }

        auto [attrNames, attrBuffer] = AttributeComputation::computeSingleAttribute(*weighted, LEVEL);
        (void)attrNames;

        ExtinctionValues<std::uint8_t> extinction(*weighted, attrBuffer);
        auto filtered = extinction.filtering(ExtinctionSelectionPolicy<float>::byTopK(1024));
        auto saliency = extinction.contourMap(ExtinctionSelectionPolicy<float>::byTopK(1024), ExtinctionContourScorePolicy::ExtinctionValue);

        requireVectorEqual(collectImageValues(filtered), collectImageValues(expectedAfterMerge),
                           "weighted extinction keep-all reconstruction after middle-slot merge");
        requireVectorEqual(collectImageValues(weighted->reconstructionImage()), collectImageValues(expectedAfterMerge),
                           "weighted reconstruction after middle-slot merge");
        requireImageShape(saliency, weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D());
    }

    verifyDeterministicTieOrdering();
    verifyFormalAttributeRejectsInvalidExtinctionValues();
    verifyTreeOfShapesIsRejected();

    return 0;
}
