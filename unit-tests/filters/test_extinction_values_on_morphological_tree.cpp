#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/AttributeComputation.hpp"
#include "mmcfilters/filters/ExtinctionValues.hpp"

#include <memory>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

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

    return 0;
}
