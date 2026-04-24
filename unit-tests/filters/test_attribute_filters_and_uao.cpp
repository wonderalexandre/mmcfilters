#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/AttributeComputedIncrementally.hpp"
#include "mmcfilters/filters/AttributeFilters.hpp"
#include "mmcfilters/filters/AttributeOpeningPrimitivesFamily.hpp"
#include "mmcfilters/filters/UltimateAttributeOpening.hpp"

#include <memory>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

int main() {
    auto image = makeComponentTreeFixture();

    for (bool isMaxtree : {true, false}) {
        auto weighted = makeWeightedComponentTree(image, isMaxtree);
        auto reconstruction = weighted->reconstructionImage();
        AttributeFilters weightedFilters(*weighted);

        std::vector<bool> keepAll(weighted->tree.getNumInternalNodeSlots(), true);

        auto directViaObject = weightedFilters.filteringByDirectRule(keepAll);
        requireVectorEqual(
            collectImageValues(directViaObject),
            collectImageValues(reconstruction),
            isMaxtree ? "weighted max-tree direct-rule keep-all via object" : "weighted min-tree direct-rule keep-all via object"
        );

        auto direct = ImageUInt8::create(weighted->tree.getNumRowsOfImage(), weighted->tree.getNumColsOfImage(), 0);
        AttributeFilters::filteringByDirectRule(*weighted, keepAll, direct);
        requireVectorEqual(
            collectImageValues(direct),
            collectImageValues(reconstruction),
            isMaxtree ? "weighted max-tree direct-rule keep-all" : "weighted min-tree direct-rule keep-all"
        );

        auto subtractive = ImageUInt8::create(weighted->tree.getNumRowsOfImage(), weighted->tree.getNumColsOfImage(), 0);
        AttributeFilters::filteringBySubtractiveRule(*weighted, keepAll, subtractive);
        requireVectorEqual(
            collectImageValues(subtractive),
            collectImageValues(reconstruction),
            isMaxtree ? "weighted max-tree subtractive-rule keep-all" : "weighted min-tree subtractive-rule keep-all"
        );

        auto pruningMin = ImageUInt8::create(weighted->tree.getNumRowsOfImage(), weighted->tree.getNumColsOfImage(), 0);
        AttributeFilters::filteringByPruningMin(*weighted, keepAll, pruningMin);
        requireVectorEqual(
            collectImageValues(pruningMin),
            collectImageValues(reconstruction),
            isMaxtree ? "weighted max-tree pruning-min keep-all" : "weighted min-tree pruning-min keep-all"
        );

        auto pruningMax = ImageUInt8::create(weighted->tree.getNumRowsOfImage(), weighted->tree.getNumColsOfImage(), 0);
        AttributeFilters::filteringByPruningMax(*weighted, keepAll, pruningMax);
        requireVectorEqual(
            collectImageValues(pruningMax),
            collectImageValues(reconstruction),
            isMaxtree ? "weighted max-tree pruning-max keep-all" : "weighted min-tree pruning-max keep-all"
        );

        auto [attrNames, attr] = AttributeComputedIncrementally::computeSingleAttribute(weighted->tree, BOX_HEIGHT);
        (void)attrNames;
        float maxCriterion = static_cast<float>(weighted->tree.getNumRowsOfImage());

        auto weightedPrimitives = std::make_shared<AttributeOpeningPrimitivesFamily>(*weighted, attr, maxCriterion);
        AttributeOpeningPrimitivesFamily weightedPrimitivesRaw(*weighted, attr.data(), maxCriterion);
        require(weightedPrimitives->getNumPrimitives() >= 1, "weighted attribute-opening primitives must expose at least one primitive");
        requireEqual(
            weightedPrimitivesRaw.getNumPrimitives(),
            weightedPrimitives->getNumPrimitives(),
            "weighted raw attribute-opening primitives must match shared-buffer construction");

        UltimateAttributeOpening weightedUao(*weighted, attr);
        UltimateAttributeOpening weightedUaoRaw(*weighted, attr.data());
        weightedUao.execute(static_cast<int>(maxCriterion));
        weightedUaoRaw.execute(static_cast<int>(maxCriterion));
        requireImageShape(weightedUao.getMaxContrastImage(), weighted->tree.getNumRowsOfImage(), weighted->tree.getNumColsOfImage());
        requireImageShape(weightedUao.getAssociatedImage(), weighted->tree.getNumRowsOfImage(), weighted->tree.getNumColsOfImage());
        requireImageShape(weightedUaoRaw.getMaxContrastImage(), weighted->tree.getNumRowsOfImage(), weighted->tree.getNumColsOfImage());
        requireImageShape(weightedUaoRaw.getAssociatedImage(), weighted->tree.getNumRowsOfImage(), weighted->tree.getNumColsOfImage());

        weighted->tree.mergeNodeIntoParent(4);
        weighted->altitude[4] = AltitudeType{};
        std::vector<bool> keepAllAfterMerge(weighted->tree.getNumInternalNodeSlots(), true);
        auto expectedAfterMerge = ImageUInt8::create(weighted->tree.getNumRowsOfImage(), weighted->tree.getNumColsOfImage(), 0);
        for (int pixel = 0; pixel < expectedAfterMerge->getSize(); ++pixel) {
            const NodeId nodeId = weighted->tree.getSmallestComponent(pixel);
            (*expectedAfterMerge)[pixel] = static_cast<uint8_t>(weighted->getAltitude(nodeId));
        }

        auto directAfterMerge = ImageUInt8::create(weighted->tree.getNumRowsOfImage(), weighted->tree.getNumColsOfImage(), 0);
        AttributeFilters::filteringByDirectRule(*weighted, keepAllAfterMerge, directAfterMerge);
        requireVectorEqual(
            collectImageValues(directAfterMerge),
            collectImageValues(expectedAfterMerge),
            isMaxtree ? "weighted max-tree direct-rule keep-all after merge" : "weighted min-tree direct-rule keep-all after merge"
        );

        auto subtractiveAfterMerge = ImageUInt8::create(weighted->tree.getNumRowsOfImage(), weighted->tree.getNumColsOfImage(), 0);
        AttributeFilters::filteringBySubtractiveRule(*weighted, keepAllAfterMerge, subtractiveAfterMerge);
        requireVectorEqual(
            collectImageValues(subtractiveAfterMerge),
            collectImageValues(expectedAfterMerge),
            isMaxtree ? "weighted max-tree subtractive-rule keep-all after merge" : "weighted min-tree subtractive-rule keep-all after merge"
        );

        auto [mergedAttrNames, mergedAttr] = AttributeComputedIncrementally::computeSingleAttribute(weighted->tree, BOX_HEIGHT);
        (void)mergedAttrNames;
        UltimateAttributeOpening mergedWeightedUao(*weighted, mergedAttr);
        UltimateAttributeOpening mergedWeightedUaoRaw(*weighted, mergedAttr.data());
        mergedWeightedUao.execute(static_cast<int>(maxCriterion));
        mergedWeightedUaoRaw.execute(static_cast<int>(maxCriterion));
        requireImageShape(mergedWeightedUao.getMaxContrastImage(), weighted->tree.getNumRowsOfImage(), weighted->tree.getNumColsOfImage());
        requireImageShape(mergedWeightedUao.getAssociatedImage(), weighted->tree.getNumRowsOfImage(), weighted->tree.getNumColsOfImage());
        requireImageShape(mergedWeightedUaoRaw.getMaxContrastImage(), weighted->tree.getNumRowsOfImage(), weighted->tree.getNumColsOfImage());
        requireImageShape(mergedWeightedUaoRaw.getAssociatedImage(), weighted->tree.getNumRowsOfImage(), weighted->tree.getNumColsOfImage());
    }

    return 0;
}
