#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/AttributeComputation.hpp"
#include "mmcfilters/filters/ExtinctionValues.hpp"

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
}

void verifyTreeOfShapesIsRejected() {
    auto image = makeComponentTreeFixture();
    auto tos = makeWeightedTreeOfShapes(image, ToSInterpolation::SelfDual);
    std::vector<float> tosAttr(static_cast<std::size_t>(tos->topology().getNumInternalNodeSlots()), 1.0f);
    requireThrows<std::invalid_argument>([&]() { ExtinctionValues<std::uint8_t> invalid(*tos, tosAttr); }, "ExtinctionValues must reject tree of shapes");
}

} // namespace

int main() {
    auto image = makeComponentTreeFixture();

    for (bool isMaxtree : {true, false}) {
        auto weighted = makeWeightedComponentTree(image, isMaxtree);
        auto [attrNames, attrBuffer] = AttributeComputation::computeSingleAttribute(*weighted, LEVEL);
        (void)attrNames;

        if (isMaxtree) {
            requireThrows<std::invalid_argument>([&]() { ExtinctionValues<std::uint8_t> invalidExtinction(*weighted, std::vector<float>{1.0f}); },
                                                 "ExtinctionValues<std::uint8_t> must reject short vector attribute buffer");
        }

        ExtinctionValues<std::uint8_t> extinction(*weighted, attrBuffer);
        const auto keepAll = ExtinctionSelectionPolicy<float>::byTopK(1024);
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(extinction.filtering(ExtinctionSelectionPolicy<float>::byTopK(-1))); },
                                             isMaxtree ? "weighted max-tree extinction filtering must reject negative keep count"
                                                       : "weighted min-tree extinction filtering must reject negative keep count");
        requireThrows<std::invalid_argument>(
            [&]() { static_cast<void>(extinction.contourMap(ExtinctionSelectionPolicy<float>::byTopK(-1), ExtinctionContourScorePolicy::RankScore)); },
            isMaxtree ? "weighted max-tree extinction contour map must reject negative keep count"
                      : "weighted min-tree extinction contour map must reject negative keep count");
        requireThrows<std::invalid_argument>(
            [&]() { static_cast<void>(extinction.filtering(ExtinctionSelectionPolicy<float>::byThreshold(std::numeric_limits<float>::quiet_NaN()))); },
            isMaxtree ? "weighted max-tree extinction threshold filtering must reject NaN"
                      : "weighted min-tree extinction threshold filtering must reject NaN");

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

        if (isMaxtree) {
            auto staleWeighted = makeWeightedComponentTree(image, true);
            auto [staleNames, staleAttr] = AttributeComputation::computeSingleAttribute(*staleWeighted, LEVEL);
            (void)staleNames;
            ExtinctionValues<std::uint8_t> staleExtinction(*staleWeighted, staleAttr);
            staleWeighted->mergeNodeIntoParent(4);
            requireThrows<std::logic_error>([&]() { static_cast<void>(staleExtinction.filtering(keepAll)); },
                                            "ExtinctionValues<std::uint8_t> filtering must reject use after topology mutation");
            requireThrows<std::logic_error>([&]() { static_cast<void>(staleExtinction.contourMap(keepAll, ExtinctionContourScorePolicy::RankScore)); },
                                            "ExtinctionValues<std::uint8_t> contour map must reject use after topology mutation");
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
        auto invalidContourBefore = collectImageValues(invalidExtinction.contourMap(keepAll, ExtinctionContourScorePolicy::ExtinctionValue));

        std::vector<int> invalidAltitude = invalidWeighted.getAltitudeBuffer();
        const int invalidOffset = isMaxtree ? 300 : -300;
        for (int& level : invalidAltitude) {
            level += invalidOffset;
        }
        invalidWeighted.setAltitudeBuffer(std::move(invalidAltitude));
        requireVectorEqual(collectImageValues(invalidExtinction.contourMap(keepAll, ExtinctionContourScorePolicy::ExtinctionValue)), invalidContourBefore,
                           isMaxtree ? "weighted max-tree extinction contour map ignores replaced altitude"
                                     : "weighted min-tree extinction contour map ignores replaced altitude");
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
        auto contourMap = extinction.contourMap(ExtinctionSelectionPolicy<float>::byTopK(1024), ExtinctionContourScorePolicy::ExtinctionValue);

        requireVectorEqual(collectImageValues(filtered), collectImageValues(expectedAfterMerge),
                           "weighted extinction keep-all reconstruction after middle-slot merge");
        requireVectorEqual(collectImageValues(weighted->reconstructionImage()), collectImageValues(expectedAfterMerge),
                           "weighted reconstruction after middle-slot merge");
        requireImageShape(contourMap, weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D());
    }

    verifyDeterministicTieOrdering();
    verifyTreeOfShapesIsRejected();

    return 0;
}
