#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/AttributeComputation.hpp"
#include "mmcfilters/filters/ExtinctionValues.hpp"

#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <vector>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

namespace {

WeightedMorphologicalTree<std::uint8_t> makeFlatThreeLeafMaxTreeFixture() {
    const std::vector<NodeId> parent{
        3, 4, 5,
        6, 6, 6, 6};
    const std::vector<std::uint8_t> altitude{
        5, 5, 5,
        5, 5, 5, 0};
    return MorphologicalTreeFactory::createFromHigraParent(
        std::span<const NodeId>(parent),
        std::span<const std::uint8_t>(altitude),
        1,
        3,
        MorphologicalTreeKind::MAX_TREE,
        AdjacencyRelation(1, 3, 1.5));
}

void verifyDeterministicTieOrdering() {
    auto weighted = makeFlatThreeLeafMaxTreeFixture();
    const MorphologicalTree& tree = weighted.topology();
    std::vector<float> attr(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), 1.0f);

    ExtinctionValues<std::uint8_t> extinction(weighted, attr);
    auto& records = extinction.getExtinctionValues();

    requireEqual(static_cast<int>(records.size()), 3, "flat extinction record count");
    requireEqual(records[0].extinction, std::numeric_limits<float>::max(), "dominant extinction sentinel");
    requireEqual(records[1].extinction, 1.0f, "first tied loser extinction");
    requireEqual(records[2].extinction, 1.0f, "second tied loser extinction");
    require(records[1].leaf < records[2].leaf, "equal extinction records must be sorted deterministically by leaf NodeId");
}

void verifySelfDualTreesAreRejected() {
    auto image = makeComponentTreeFixture();
    auto tos = makeWeightedTreeOfShapes(image, ToSInterpolation::SelfDual);
    std::vector<float> tosAttr(static_cast<std::size_t>(tos->topology().getNumInternalNodeSlots()), 1.0f);
    requireThrows<std::invalid_argument>(
        [&]() { ExtinctionValues<std::uint8_t> invalid(*tos, tosAttr); },
        "ExtinctionValues must reject tree of shapes");

    auto residualImage = makeImage(1, 2, {0, 1});
    auto residual = MorphologicalTreeFactory::createSelfDualResidualTree(residualImage, 1.0);
    std::vector<float> residualAttr(static_cast<std::size_t>(residual.topology().getNumInternalNodeSlots()), 1.0f);
    requireThrows<std::invalid_argument>(
        [&]() { ExtinctionValues<std::uint8_t> invalid(residual, residualAttr); },
        "ExtinctionValues must reject self-dual residual trees");
}

} // namespace

int main() {
    auto image = makeComponentTreeFixture();

    for (bool isMaxtree : {true, false}) {
        auto weighted = makeWeightedComponentTree(image, isMaxtree);
        auto [attrNames, attrBuffer] = AttributeComputation::computeSingleAttribute(*weighted, LEVEL);
        (void)attrNames;

        if (isMaxtree) {
            requireThrows<std::invalid_argument>(
                [&]() { ExtinctionValues<std::uint8_t> invalidExtinction(*weighted, std::vector<float>{1.0f}); },
                "ExtinctionValues<std::uint8_t> must reject short vector attribute buffer");
        }

        ExtinctionValues<std::uint8_t> extinction(*weighted, attrBuffer);
        requireThrows<std::invalid_argument>(
            [&]() { static_cast<void>(extinction.filtering(-1)); },
            isMaxtree ? "weighted max-tree extinction filtering must reject negative keep count" : "weighted min-tree extinction filtering must reject negative keep count");
        requireThrows<std::invalid_argument>(
            [&]() { static_cast<void>(extinction.saliencyMap(-1)); },
            isMaxtree ? "weighted max-tree extinction saliency must reject negative keep count" : "weighted min-tree extinction saliency must reject negative keep count");
        auto filtered = extinction.filtering(1024);
        auto reconstructed = weighted->reconstructionImage();

        requireVectorEqual(
            collectImageValues(filtered),
            collectImageValues(reconstructed),
            isMaxtree ? "weighted max-tree extinction keep-all reconstruction" : "weighted min-tree extinction keep-all reconstruction"
        );

        if (isMaxtree) {
            auto staleWeighted = makeWeightedComponentTree(image, true);
            auto [staleNames, staleAttr] = AttributeComputation::computeSingleAttribute(*staleWeighted, LEVEL);
            (void)staleNames;
            ExtinctionValues<std::uint8_t> staleExtinction(*staleWeighted, staleAttr);
            staleWeighted->mergeNodeIntoParent(4);
            requireThrows<std::logic_error>(
                [&]() { static_cast<void>(staleExtinction.filtering(1024)); },
                "ExtinctionValues<std::uint8_t> filtering must reject use after topology mutation");
            requireThrows<std::logic_error>(
                [&]() { static_cast<void>(staleExtinction.saliencyMap(1024)); },
                "ExtinctionValues<std::uint8_t> saliency must reject use after topology mutation");
            requireThrows<std::logic_error>(
                [&]() { static_cast<void>(staleExtinction.getExtinctionValues()); },
                "ExtinctionValues<std::uint8_t> tuples must reject use after topology mutation");
        }

        const auto [invalidParent, exportedAltitude] = weighted->exportHigraHierarchy();
        std::vector<int> validHigraAltitude(exportedAltitude.begin(), exportedAltitude.end());
        auto invalidWeighted = MorphologicalTreeFactory::createFromHigraParent<int>(
            std::span<const NodeId>(invalidParent),
            std::span<const int>(validHigraAltitude),
            image->getNumRows(),
            image->getNumCols(),
            isMaxtree ? MorphologicalTreeKind::MAX_TREE : MorphologicalTreeKind::MIN_TREE,
            AdjacencyRelation(image->getNumRows(), image->getNumCols(), 1.5));
        auto [invalidLevelNames, invalidLevelBuffer] =
            AttributeComputation::computeSingleAttribute(invalidWeighted.asView(), LEVEL);
        (void)invalidLevelNames;
        ExtinctionValues<int> invalidExtinction(invalidWeighted, invalidLevelBuffer);
        auto invalidSaliencyBefore = collectImageValues(invalidExtinction.saliencyMap(1024, false));

        std::vector<int> invalidAltitude = invalidWeighted.getAltitudeBuffer();
        const int invalidOffset = isMaxtree ? 300 : -300;
        for (int& level : invalidAltitude) {
            level += invalidOffset;
        }
        invalidWeighted.setAltitudeBuffer(std::move(invalidAltitude));
        requireVectorEqual(
            collectImageValues(invalidExtinction.saliencyMap(1024, false)),
            invalidSaliencyBefore,
            isMaxtree ? "weighted max-tree extinction saliency ignores replaced altitude" : "weighted min-tree extinction saliency ignores replaced altitude"
        );
        requireVectorEqual(
            collectImageValues(invalidExtinction.filtering(1024)),
            collectImageValues(invalidWeighted.reconstructionImage()),
            isMaxtree ? "ExtinctionValues<std::uint8_t> object must preserve altitude above 255" : "ExtinctionValues<std::uint8_t> object must preserve negative altitude");
    }

    {
        auto weighted = makeWeightedComponentTree(image, true);
        weighted->mergeNodeIntoParent(4);

        auto expectedAfterMerge = ImageUInt8::create(weighted->topology().getNumRowsOfImage(), weighted->topology().getNumColsOfImage(), 0);
        for (int pixel = 0; pixel < expectedAfterMerge->getSize(); ++pixel) {
            const NodeId nodeId = weighted->topology().getProperPartOwner(pixel);
            (*expectedAfterMerge)[pixel] = static_cast<uint8_t>(weighted->getAltitude(nodeId));
        }

        auto [attrNames, attrBuffer] = AttributeComputation::computeSingleAttribute(*weighted, LEVEL);
        (void)attrNames;

        ExtinctionValues<std::uint8_t> extinction(*weighted, attrBuffer);
        auto filtered = extinction.filtering(1024);
        auto saliency = extinction.saliencyMap(1024, false);

        requireVectorEqual(
            collectImageValues(filtered),
            collectImageValues(expectedAfterMerge),
            "weighted extinction keep-all reconstruction after middle-slot merge");
        requireVectorEqual(
            collectImageValues(weighted->reconstructionImage()),
            collectImageValues(expectedAfterMerge),
            "weighted reconstruction after middle-slot merge"
        );
        requireImageShape(saliency, weighted->topology().getNumRowsOfImage(), weighted->topology().getNumColsOfImage());
    }

    verifyDeterministicTieOrdering();
    verifySelfDualTreesAreRejected();

    return 0;
}
