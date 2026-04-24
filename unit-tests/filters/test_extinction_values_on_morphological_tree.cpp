#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/AttributeComputedIncrementally.hpp"
#include "mmcfilters/filters/ExtinctionValues.hpp"

#include <memory>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

int main() {
    auto image = makeComponentTreeFixture();

    for (bool isMaxtree : {true, false}) {
        auto weighted = std::make_shared<WeightedMorphologicalTree>(image, isMaxtree);
        auto [attrNames, attrBuffer] = AttributeComputedIncrementally::computeSingleAttribute(*weighted, LEVEL);
        (void)attrNames;

        ExtinctionValues extinction(*weighted, attrBuffer);
        auto filtered = extinction.filtering(1024);
        auto reconstructed = weighted->reconstructionImage();

        requireVectorEqual(
            collectImageValues(filtered),
            collectImageValues(reconstructed),
            isMaxtree ? "weighted max-tree extinction keep-all reconstruction" : "weighted min-tree extinction keep-all reconstruction"
        );
    }

    {
        auto weighted = std::make_shared<WeightedMorphologicalTree>(image, true);
        weighted->tree.mergeNodeIntoParent(4);

        auto expectedAfterMerge = ImageUInt8::create(weighted->tree.getNumRowsOfImage(), weighted->tree.getNumColsOfImage(), 0);
        for (int pixel = 0; pixel < expectedAfterMerge->getSize(); ++pixel) {
            const NodeId nodeId = weighted->tree.getSmallestComponent(pixel);
            (*expectedAfterMerge)[pixel] = static_cast<uint8_t>(weighted->getAltitude(nodeId));
        }

        auto [attrNames, attrBuffer] = AttributeComputedIncrementally::computeSingleAttribute(*weighted, LEVEL);
        (void)attrNames;

        ExtinctionValues extinction(*weighted, attrBuffer);
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
        requireImageShape(saliency, weighted->tree.getNumRowsOfImage(), weighted->tree.getNumColsOfImage());
    }

    return 0;
}
