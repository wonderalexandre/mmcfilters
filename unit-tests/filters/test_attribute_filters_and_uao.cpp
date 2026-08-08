#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/AttributeComputation.hpp"
#include "mmcfilters/filters/AttributeFilters.hpp"
#include "mmcfilters/filters/detail/ViterbiDecision.hpp"
#include "mmcfilters/filters/ExtinctionValues.hpp"
#include "mmcfilters/filters/UltimateAttributeOpening.hpp"
#include "mmcfilters/trees/WeightedTreeView.hpp"

#include <cstdint>
#include <cmath>
#include <limits>
#include <memory>
#include <span>
#include <type_traits>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

template <class OutputValue, class ImagePtrT> std::vector<OutputValue> collectImageValuesAs(const ImagePtrT& image) {
    std::vector<OutputValue> output;
    output.reserve(static_cast<std::size_t>(image->getSize()));
    for (auto value : collectImageValues(image)) {
        output.push_back(static_cast<OutputValue>(value));
    }
    return output;
}

template <class T> void writeReferenceProperParts(const MorphologicalTree& tree, NodeId nodeId, std::shared_ptr<Image<T>>& image, T value) {
    for (int pixel : tree.getProperParts(nodeId)) {
        (*image)[pixel] = value;
    }
}

template <class T>
void propagateDirectReference(const WeightedMorphologicalTree<T>& weighted, NodeId nodeId, const std::vector<bool>& criterion, std::vector<T>& mapLevel) {
    const MorphologicalTree& tree = weighted.topology();
    for (NodeId childNodeId : tree.getChildren(nodeId)) {
        mapLevel[childNodeId] = criterion[childNodeId] ? weighted.getAltitude(childNodeId) : mapLevel[nodeId];
        propagateDirectReference(weighted, childNodeId, criterion, mapLevel);
    }
}

template <class T> std::vector<T> directReferenceImage(const WeightedMorphologicalTree<T>& weighted, const std::vector<bool>& criterion) {
    const MorphologicalTree& tree = weighted.topology();
    std::vector<T> mapLevel(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), T{});
    const NodeId rootNodeId = tree.getRoot();
    mapLevel[rootNodeId] = weighted.getAltitude(rootNodeId);
    propagateDirectReference(weighted, rootNodeId, criterion, mapLevel);

    auto image = Image<T>::create(tree.getNumRowsOfGridDomain2D(), tree.getNumColsOfGridDomain2D(), T{});
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        writeReferenceProperParts(tree, nodeId, image, mapLevel[nodeId]);
    }
    return collectImageValues(image);
}

template <class T>
void propagateSubtractiveReference(const WeightedMorphologicalTree<T>& weighted, NodeId nodeId, const std::vector<bool>& criterion,
                                   std::vector<AltitudeDiff<T>>& mapLevel) {
    const MorphologicalTree& tree = weighted.topology();
    for (NodeId childNodeId : tree.getChildren(nodeId)) {
        mapLevel[childNodeId] =
            criterion[childNodeId] ? static_cast<AltitudeDiff<T>>(mapLevel[nodeId] + weighted.getNodeResidue(childNodeId)) : mapLevel[nodeId];
        propagateSubtractiveReference(weighted, childNodeId, criterion, mapLevel);
    }
}

template <class T> std::vector<T> subtractiveReferenceImage(const WeightedMorphologicalTree<T>& weighted, const std::vector<bool>& criterion) {
    const MorphologicalTree& tree = weighted.topology();
    std::vector<AltitudeDiff<T>> mapLevel(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), AltitudeDiff<T>{});
    const NodeId rootNodeId = tree.getRoot();
    mapLevel[rootNodeId] = static_cast<AltitudeDiff<T>>(weighted.getAltitude(rootNodeId));
    propagateSubtractiveReference(weighted, rootNodeId, criterion, mapLevel);

    auto image = Image<T>::create(tree.getNumRowsOfGridDomain2D(), tree.getNumColsOfGridDomain2D(), T{});
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        writeReferenceProperParts(tree, nodeId, image, static_cast<T>(mapLevel[nodeId]));
    }
    return collectImageValues(image);
}

WeightedMorphologicalTree<std::uint8_t> makeViterbiChainFixture() {
    const std::vector<NodeId> parent = {1, 2, 3, 3};
    const std::vector<std::uint8_t> altitude = {5, 5, 3, 0};
    return MorphologicalTreeFactory::createFromHigraParent<std::uint8_t>(std::span<const NodeId>(parent), std::span<const std::uint8_t>(altitude), 1, 1,
                                                                         MorphologicalTreeKind::MAX_TREE, RegularGridAdjacency2D(1, 1, 1.5));
}

WeightedMorphologicalTree<std::uint8_t> makeViterbiBranchFixture() {
    const std::vector<NodeId> parent = {3, 4, 4, 5, 5, 5};
    const std::vector<std::uint8_t> altitude = {5, 4, 4, 5, 4, 0};
    return MorphologicalTreeFactory::createFromHigraParent<std::uint8_t>(std::span<const NodeId>(parent), std::span<const std::uint8_t>(altitude), 1, 3,
                                                                         MorphologicalTreeKind::MAX_TREE, RegularGridAdjacency2D(1, 3, 1.5));
}

int main() {
    {
        auto chain = makeViterbiChainFixture();
        AttributeFilters<std::uint8_t> filters(chain);
        std::vector<float> attr = {1.1f, 0.0f, 2.0f};
        auto costs = detail::makeThresholdViterbiCosts(chain.topology(), attr.data(), 1.0f);
        auto keep = detail::computeViterbiKeepCriterion(chain.topology(), costs);
        require(!keep[0] && !keep[1] && keep[2], "Viterbi chain must remove descendants below a removed ancestor");
        requireVectorEqual(collectImageValues(filters.filteringByViterbiRule(attr.data(), 1.0f)), std::vector<std::uint8_t>{0},
                           "Viterbi chain reconstruction with removed middle node");

        attr = {10.0f, 0.9f, 2.0f};
        keep = detail::computeViterbiKeepCriterion(chain.topology(), detail::makeThresholdViterbiCosts(chain.topology(), attr.data(), 1.0f));
        require(keep[0] && keep[1] && keep[2], "Viterbi chain must preserve a costly-to-remove descendant branch");
        requireVectorEqual(collectImageValues(filters.filteringByViterbiRule(attr.data(), 1.0f)), std::vector<std::uint8_t>{5},
                           "Viterbi chain reconstruction with preserved branch");

        attr = {1.0f, 1.0f, 1.0f};
        costs = detail::makeThresholdViterbiCosts(chain.topology(), attr.data(), 1.0f);
        keep = detail::computeViterbiKeepCriterion(chain.topology(), costs);
        require(!keep[0] && !keep[1] && keep[2], "Viterbi default tie-break must prefer removal except at the forced root");

        detail::ViterbiDecisionOptions preserveTies;
        preserveTies.tieBreak = detail::ViterbiTieBreak::PreferPreserve;
        keep = detail::computeViterbiKeepCriterion(chain.topology(), costs, preserveTies);
        require(keep[0] && keep[1] && keep[2], "Viterbi preserve tie-break must keep the full chain");

        std::vector<float> nanAttr = {1.0f, std::numeric_limits<float>::quiet_NaN(), 1.0f};
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(filters.filteringByViterbiRule(nanAttr.data(), 1.0f)); },
                                             "Viterbi rule must reject NaN attributes");
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(filters.filteringByViterbiRule(attr.data(), std::numeric_limits<float>::quiet_NaN())); },
                                             "Viterbi rule must reject NaN thresholds");
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(filters.filteringByViterbiRule(static_cast<const float*>(nullptr), 1.0f)); },
                                             "Viterbi rule must reject null attribute buffers");
    }

    {
        auto branch = makeViterbiBranchFixture();
        AttributeFilters<std::uint8_t> filters(branch);
        std::vector<double> attr = {2.0, 0.0, 2.0};
        auto keep = detail::computeViterbiKeepCriterion(branch.topology(), detail::makeThresholdViterbiCosts(branch.topology(), attr.data(), 1.0));
        require(keep[0] && !keep[1] && keep[2], "Viterbi branch must decide sibling subtrees independently");
        requireVectorEqual(collectImageValues(filters.filteringByViterbiRule(attr.data(), 1.0)), std::vector<std::uint8_t>{5, 0, 0},
                           "Viterbi branch reconstruction");
    }

    auto image = makeComponentTreeFixture();

    for (bool isMaxtree : {true, false}) {
        auto weighted = makeWeightedComponentTree(image, isMaxtree);
        auto reconstruction = weighted->reconstructionImage();
        AttributeFilters<std::uint8_t> weightedFilters(*weighted);
        const auto weightedView = weighted->asView();

        std::vector<bool> keepAll(weighted->topology().getNumInternalNodeSlots(), true);
        if (isMaxtree) {
            std::vector<bool> shortCriterion(1, true);
            std::vector<float> shortScores(1, 0.0f);
            auto wrongShape = ImageUInt8::create(1, 1, 0);
            requireThrows<std::invalid_argument>([&]() { static_cast<void>(weightedFilters.filteringByDirectRule(shortCriterion)); },
                                                 "AttributeFilters<std::uint8_t> object must reject short criterion");
            requireThrows<std::invalid_argument>([&]() { static_cast<void>(weightedFilters.filteringBySubtractiveScoreRule(shortScores)); },
                                                 "AttributeFilters<std::uint8_t> object must reject short score buffer");
            requireThrows<std::invalid_argument>([&]() { AttributeFilters<std::uint8_t>::filteringByDirectRule(*weighted, keepAll, wrongShape); },
                                                 "AttributeFilters<std::uint8_t> static direct rule must reject wrong output image shape");
            requireThrows<std::invalid_argument>([&]() { AttributeFilters<std::uint8_t>::filteringByDirectRule(weightedView, keepAll, wrongShape); },
                                                 "AttributeFilters<std::uint8_t> view direct rule must reject wrong output image shape");
            requireThrows<std::invalid_argument>(
                [&]() { AttributeFilters<std::uint8_t>::filteringByPruningMin(*weighted, static_cast<const float*>(nullptr), 1.0f, wrongShape); },
                "AttributeFilters<std::uint8_t> static pruning must reject null attribute pointer");

            auto staleWeighted = makeWeightedComponentTree(image, true);
            auto [staleNames, staleAttr] = AttributeComputation::computeSingleAttribute(*staleWeighted, LEVEL);
            (void)staleNames;
            AttributeFilters<std::uint8_t> staleFilters(*staleWeighted);
            UltimateAttributeOpening<std::uint8_t> staleUao(*staleWeighted, staleAttr);
            std::vector<bool> staleKeepAll(staleWeighted->topology().getNumInternalNodeSlots(), true);
            staleWeighted->mergeNodeIntoParent(4);
            requireThrows<std::logic_error>([&]() { static_cast<void>(staleFilters.filteringByDirectRule(staleKeepAll)); },
                                            "AttributeFilters<std::uint8_t> object must reject use after topology mutation");
            requireThrows<std::logic_error>([&]() { staleUao.execute(4); },
                                            "UltimateAttributeOpening<std::uint8_t> execute must reject use after topology mutation");
            requireThrows<std::logic_error>([&]() { static_cast<void>(staleUao.getMaxContrastImage()); },
                                            "UltimateAttributeOpening<std::uint8_t> output must reject use after topology mutation");

            auto staleViewWeighted = makeWeightedComponentTree(image, true);
            const auto staleView = staleViewWeighted->asView();
            std::vector<bool> staleViewKeepAll(staleViewWeighted->topology().getNumInternalNodeSlots(), true);
            auto staleViewOutput =
                ImageUInt8::create(staleViewWeighted->topology().getNumRowsOfGridDomain2D(), staleViewWeighted->topology().getNumColsOfGridDomain2D(), 0);
            staleViewWeighted->mergeNodeIntoParent(4);
            requireThrows<std::logic_error>([&]() { AttributeFilters<std::uint8_t>::filteringByDirectRule(staleView, staleViewKeepAll, staleViewOutput); },
                                            "AttributeFilters<std::uint8_t> static view API must reject stale WeightedTreeView");
            requireThrows<std::logic_error>([&]() { static_cast<void>(AttributeComputation::computeSingleAttribute(staleView, LEVEL)); },
                                            "AttributeComputation must reject stale WeightedTreeView");
        }

        auto directViaObject = weightedFilters.filteringByDirectRule(keepAll);
        requireVectorEqual(collectImageValues(directViaObject), collectImageValues(reconstruction),
                           isMaxtree ? "weighted max-tree direct-rule keep-all via object" : "weighted min-tree direct-rule keep-all via object");

        auto direct = ImageUInt8::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        AttributeFilters<std::uint8_t>::filteringByDirectRule(*weighted, keepAll, direct);
        requireVectorEqual(collectImageValues(direct), collectImageValues(reconstruction),
                           isMaxtree ? "weighted max-tree direct-rule keep-all" : "weighted min-tree direct-rule keep-all");

        auto directViaView = ImageUInt8::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        AttributeFilters<std::uint8_t>::filteringByDirectRule(weightedView, keepAll, directViaView);
        requireVectorEqual(collectImageValues(directViaView), collectImageValues(reconstruction),
                           isMaxtree ? "weighted max-tree direct-rule keep-all via view" : "weighted min-tree direct-rule keep-all via view");

        AltitudeBuffer<std::uint8_t> externalAltitude = weighted->getAltitudeBuffer();
        const WeightedTreeView<std::uint8_t> externalView(weighted->topology(), std::span<const std::uint8_t>(externalAltitude));
        std::vector<std::int16_t> int16Altitude;
        int16Altitude.reserve(externalAltitude.size());
        for (std::uint8_t level : externalAltitude) {
            int16Altitude.push_back(static_cast<std::int16_t>(level));
        }
        const WeightedTreeView<std::int16_t> int16View(weighted->topology(), std::span<const std::int16_t>(int16Altitude.data(), int16Altitude.size()));
        auto directViaExternalView = ImageUInt8::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        AttributeFilters<std::uint8_t>::filteringByDirectRule(externalView, keepAll, directViaExternalView);
        requireVectorEqual(collectImageValues(directViaExternalView), collectImageValues(reconstruction),
                           isMaxtree ? "weighted max-tree direct-rule keep-all via external view" : "weighted min-tree direct-rule keep-all via external view");

        AttributeFilters<std::uint8_t> externalViewFilters(externalView);
        auto directViaExternalViewObject = externalViewFilters.filteringByDirectRule(keepAll);
        requireVectorEqual(collectImageValues(directViaExternalViewObject), collectImageValues(reconstruction),
                           isMaxtree ? "weighted max-tree direct-rule keep-all via external view object"
                                     : "weighted min-tree direct-rule keep-all via external view object");

        auto directViaInt16View =
            Image<std::int16_t>::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        AttributeFilters<std::int16_t>::filteringByDirectRule(int16View, keepAll, directViaInt16View);
        requireVectorEqual(collectImageValues(directViaInt16View), collectImageValuesAs<std::int16_t>(reconstruction),
                           isMaxtree ? "weighted max-tree direct-rule keep-all via int16 view" : "weighted min-tree direct-rule keep-all via int16 view");

        AttributeFilters<std::int16_t> int16ViewFilters(int16View);
        auto directViaInt16ViewObject = int16ViewFilters.filteringByDirectRule(keepAll);
        requireVectorEqual(collectImageValues(directViaInt16ViewObject), collectImageValuesAs<std::int16_t>(reconstruction),
                           isMaxtree ? "weighted max-tree direct-rule keep-all via int16 view object"
                                     : "weighted min-tree direct-rule keep-all via int16 view object");

        std::vector<float> unitScores(weighted->topology().getNumInternalNodeSlots(), 1.0f);
        auto scoreViaObject = weightedFilters.filteringBySubtractiveScoreRule(unitScores);
        auto scoreViaView = ImageFloat::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        auto scoreViaExternalView = ImageFloat::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        auto scoreViaInt16View = ImageFloat::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        AttributeFilters<std::uint8_t>::filteringBySubtractiveScoreRule(weightedView, unitScores, scoreViaView);
        AttributeFilters<std::uint8_t>::filteringBySubtractiveScoreRule(externalView, unitScores, scoreViaExternalView);
        AttributeFilters<std::int16_t>::filteringBySubtractiveScoreRule(int16View, unitScores, scoreViaInt16View);
        requireVectorEqual(collectImageValues(scoreViaView), collectImageValues(scoreViaObject),
                           isMaxtree ? "weighted max-tree subtractive-score via view" : "weighted min-tree subtractive-score via view");
        requireVectorEqual(collectImageValues(scoreViaExternalView), collectImageValues(scoreViaObject),
                           isMaxtree ? "weighted max-tree subtractive-score via external view" : "weighted min-tree subtractive-score via external view");
        requireVectorEqual(collectImageValues(scoreViaInt16View), collectImageValues(scoreViaObject),
                           isMaxtree ? "weighted max-tree subtractive-score via int16 view" : "weighted min-tree subtractive-score via int16 view");

        {
            const auto [higraParent, higraAltitude] = weighted->exportHigraHierarchy();
            auto importedWeighted = MorphologicalTreeFactory::createFromHigraParent<std::uint8_t>(
                std::span<const NodeId>(higraParent), std::span<const std::uint8_t>(higraAltitude), image->getNumRows(), image->getNumCols(),
                isMaxtree ? MorphologicalTreeKind::MAX_TREE : MorphologicalTreeKind::MIN_TREE,
                RegularGridAdjacency2D(image->getNumRows(), image->getNumCols(), 1.5));
            AttributeFilters<std::uint8_t> importedFilters(importedWeighted);
            std::vector<bool> importedKeepAll(importedWeighted.topology().getNumInternalNodeSlots(), true);
            auto importedReconstruction = importedWeighted.reconstructionImage();

            requireVectorEqual(collectImageValues(importedFilters.filteringBySubtractiveRule(importedKeepAll)), collectImageValues(importedReconstruction),
                               isMaxtree ? "imported max-tree subtractive-rule keep-all" : "imported min-tree subtractive-rule keep-all");

            std::vector<float> importedUnitScores(importedWeighted.topology().getNumInternalNodeSlots(), 1.0f);
            requireVectorEqual(collectImageValues(importedFilters.filteringBySubtractiveScoreRule(importedUnitScores)),
                               collectImageValuesAs<float>(importedReconstruction),
                               isMaxtree ? "imported max-tree subtractive-score unit scores" : "imported min-tree subtractive-score unit scores");

            std::vector<bool> mixedCriterion(importedWeighted.topology().getNumInternalNodeSlots(), true);
            if (!mixedCriterion.empty()) {
                mixedCriterion[0] = false;
            }
            if (mixedCriterion.size() > 2) {
                mixedCriterion[2] = false;
            }

            requireVectorEqual(collectImageValues(importedFilters.filteringByDirectRule(mixedCriterion)),
                               directReferenceImage(importedWeighted, mixedCriterion),
                               isMaxtree ? "imported max-tree direct-rule mixed criterion" : "imported min-tree direct-rule mixed criterion");
            requireVectorEqual(collectImageValues(importedFilters.filteringBySubtractiveRule(mixedCriterion)),
                               subtractiveReferenceImage(importedWeighted, mixedCriterion),
                               isMaxtree ? "imported max-tree subtractive-rule mixed criterion" : "imported min-tree subtractive-rule mixed criterion");
        }

        auto subtractive = ImageUInt8::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        AttributeFilters<std::uint8_t>::filteringBySubtractiveRule(*weighted, keepAll, subtractive);
        requireVectorEqual(collectImageValues(subtractive), collectImageValues(reconstruction),
                           isMaxtree ? "weighted max-tree subtractive-rule keep-all" : "weighted min-tree subtractive-rule keep-all");

        auto subtractiveViaView = ImageUInt8::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        AttributeFilters<std::uint8_t>::filteringBySubtractiveRule(weightedView, keepAll, subtractiveViaView);
        requireVectorEqual(collectImageValues(subtractiveViaView), collectImageValues(subtractive),
                           isMaxtree ? "weighted max-tree subtractive-rule via view" : "weighted min-tree subtractive-rule via view");

        auto subtractiveViaExternalView =
            ImageUInt8::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        auto subtractiveViaInt16View =
            Image<std::int16_t>::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        AttributeFilters<std::uint8_t>::filteringBySubtractiveRule(externalView, keepAll, subtractiveViaExternalView);
        AttributeFilters<std::int16_t>::filteringBySubtractiveRule(int16View, keepAll, subtractiveViaInt16View);
        requireVectorEqual(collectImageValues(subtractiveViaExternalView), collectImageValues(subtractive),
                           isMaxtree ? "weighted max-tree subtractive-rule via external view" : "weighted min-tree subtractive-rule via external view");
        requireVectorEqual(collectImageValues(subtractiveViaInt16View), collectImageValuesAs<std::int16_t>(subtractive),
                           isMaxtree ? "weighted max-tree subtractive-rule via int16 view" : "weighted min-tree subtractive-rule via int16 view");

        auto pruningMin = ImageUInt8::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        AttributeFilters<std::uint8_t>::filteringByPruningMin(*weighted, keepAll, pruningMin);
        requireVectorEqual(collectImageValues(pruningMin), collectImageValues(reconstruction),
                           isMaxtree ? "weighted max-tree pruning-min keep-all" : "weighted min-tree pruning-min keep-all");

        auto pruningMinViaView = ImageUInt8::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        auto pruningMinViaExternalView =
            ImageUInt8::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        auto pruningMinViaInt16View =
            Image<std::int16_t>::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        AttributeFilters<std::uint8_t>::filteringByPruningMin(weightedView, keepAll, pruningMinViaView);
        AttributeFilters<std::uint8_t>::filteringByPruningMin(externalView, keepAll, pruningMinViaExternalView);
        AttributeFilters<std::int16_t>::filteringByPruningMin(int16View, keepAll, pruningMinViaInt16View);
        requireVectorEqual(collectImageValues(pruningMinViaView), collectImageValues(pruningMin),
                           isMaxtree ? "weighted max-tree pruning-min via view" : "weighted min-tree pruning-min via view");
        requireVectorEqual(collectImageValues(pruningMinViaExternalView), collectImageValues(pruningMin),
                           isMaxtree ? "weighted max-tree pruning-min via external view" : "weighted min-tree pruning-min via external view");
        requireVectorEqual(collectImageValues(pruningMinViaInt16View), collectImageValuesAs<std::int16_t>(pruningMin),
                           isMaxtree ? "weighted max-tree pruning-min via int16 view" : "weighted min-tree pruning-min via int16 view");

        auto pruningMax = ImageUInt8::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        AttributeFilters<std::uint8_t>::filteringByPruningMax(*weighted, keepAll, pruningMax);
        requireVectorEqual(collectImageValues(pruningMax), collectImageValues(reconstruction),
                           isMaxtree ? "weighted max-tree pruning-max keep-all" : "weighted min-tree pruning-max keep-all");

        auto pruningMaxViaView = ImageUInt8::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        auto pruningMaxViaExternalView =
            ImageUInt8::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        auto pruningMaxViaInt16View =
            Image<std::int16_t>::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        AttributeFilters<std::uint8_t>::filteringByPruningMax(weightedView, keepAll, pruningMaxViaView);
        AttributeFilters<std::uint8_t>::filteringByPruningMax(externalView, keepAll, pruningMaxViaExternalView);
        AttributeFilters<std::int16_t>::filteringByPruningMax(int16View, keepAll, pruningMaxViaInt16View);
        requireVectorEqual(collectImageValues(pruningMaxViaView), collectImageValues(pruningMax),
                           isMaxtree ? "weighted max-tree pruning-max via view" : "weighted min-tree pruning-max via view");
        requireVectorEqual(collectImageValues(pruningMaxViaExternalView), collectImageValues(pruningMax),
                           isMaxtree ? "weighted max-tree pruning-max via external view" : "weighted min-tree pruning-max via external view");
        requireVectorEqual(collectImageValues(pruningMaxViaInt16View), collectImageValuesAs<std::int16_t>(pruningMax),
                           isMaxtree ? "weighted max-tree pruning-max via int16 view" : "weighted min-tree pruning-max via int16 view");

        auto [attrNames, attr] = AttributeComputation::computeSingleAttribute(*weighted, BOX_HEIGHT);
        (void)attrNames;
        auto [attrNames64, attr64] = AttributeComputation::computeSingleAttribute<double>(*weighted, BOX_HEIGHT);
        (void)attrNames64;
        float maxCriterion = static_cast<float>(weighted->topology().getNumRowsOfGridDomain2D());
        double maxCriterion64 = static_cast<double>(weighted->topology().getNumRowsOfGridDomain2D());

        for (float threshold : {1.0f, 2.0f, 3.0f, 4.0f}) {
            std::vector<bool> keepByAttribute(weighted->topology().getNumInternalNodeSlots(), false);
            for (NodeId nodeId : weighted->topology().getAliveNodeIds()) {
                keepByAttribute[nodeId] = attr[static_cast<std::size_t>(nodeId)] > threshold;
            }

            requireVectorEqual(collectImageValues(weightedFilters.filteringByPruningMin(keepByAttribute)),
                               collectImageValues(weightedFilters.filteringByPruningMin(attr.data(), threshold)),
                               isMaxtree ? "weighted max-tree pruning-min criterion/attribute equivalence"
                                         : "weighted min-tree pruning-min criterion/attribute equivalence");
            requireVectorEqual(collectImageValues(weightedFilters.filteringByPruningMax(keepByAttribute)),
                               collectImageValues(weightedFilters.filteringByPruningMax(attr.data(), threshold)),
                               isMaxtree ? "weighted max-tree pruning-max criterion/attribute equivalence"
                                         : "weighted min-tree pruning-max criterion/attribute equivalence");
        }

        auto pruningMinAttrWeighted = ImageUInt8::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        auto pruningMinAttrView = ImageUInt8::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        auto pruningMinAttrExternalView =
            ImageUInt8::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        auto pruningMinAttrInt16View =
            Image<std::int16_t>::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        AttributeFilters<std::uint8_t>::filteringByPruningMin(*weighted, attr.data(), maxCriterion, pruningMinAttrWeighted);
        AttributeFilters<std::uint8_t>::filteringByPruningMin(weightedView, attr.data(), maxCriterion, pruningMinAttrView);
        AttributeFilters<std::uint8_t>::filteringByPruningMin(externalView, attr.data(), maxCriterion, pruningMinAttrExternalView);
        AttributeFilters<std::int16_t>::filteringByPruningMin(int16View, attr.data(), maxCriterion, pruningMinAttrInt16View);
        requireVectorEqual(collectImageValues(pruningMinAttrView), collectImageValues(pruningMinAttrWeighted),
                           isMaxtree ? "weighted max-tree pruning-min attr via view" : "weighted min-tree pruning-min attr via view");
        requireVectorEqual(collectImageValues(pruningMinAttrExternalView), collectImageValues(pruningMinAttrWeighted),
                           isMaxtree ? "weighted max-tree pruning-min attr via external view" : "weighted min-tree pruning-min attr via external view");
        requireVectorEqual(collectImageValues(pruningMinAttrInt16View), collectImageValuesAs<std::int16_t>(pruningMinAttrWeighted),
                           isMaxtree ? "weighted max-tree pruning-min attr via int16 view" : "weighted min-tree pruning-min attr via int16 view");

        auto pruningMaxAttrWeighted = ImageUInt8::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        auto pruningMaxAttrView = ImageUInt8::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        auto pruningMaxAttrExternalView =
            ImageUInt8::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        auto pruningMaxAttrInt16View =
            Image<std::int16_t>::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        AttributeFilters<std::uint8_t>::filteringByPruningMax(*weighted, attr.data(), maxCriterion, pruningMaxAttrWeighted);
        AttributeFilters<std::uint8_t>::filteringByPruningMax(weightedView, attr.data(), maxCriterion, pruningMaxAttrView);
        AttributeFilters<std::uint8_t>::filteringByPruningMax(externalView, attr.data(), maxCriterion, pruningMaxAttrExternalView);
        AttributeFilters<std::int16_t>::filteringByPruningMax(int16View, attr.data(), maxCriterion, pruningMaxAttrInt16View);
        requireVectorEqual(collectImageValues(pruningMaxAttrView), collectImageValues(pruningMaxAttrWeighted),
                           isMaxtree ? "weighted max-tree pruning-max attr via view" : "weighted min-tree pruning-max attr via view");
        requireVectorEqual(collectImageValues(pruningMaxAttrExternalView), collectImageValues(pruningMaxAttrWeighted),
                           isMaxtree ? "weighted max-tree pruning-max attr via external view" : "weighted min-tree pruning-max attr via external view");
        requireVectorEqual(collectImageValues(pruningMaxAttrInt16View), collectImageValuesAs<std::int16_t>(pruningMaxAttrWeighted),
                           isMaxtree ? "weighted max-tree pruning-max attr via int16 view" : "weighted min-tree pruning-max attr via int16 view");

        auto pruningMinAttr64 = ImageUInt8::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        auto pruningMaxAttr64 = ImageUInt8::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        AttributeFilters<std::uint8_t>::filteringByPruningMin(*weighted, attr64.data(), maxCriterion64, pruningMinAttr64);
        AttributeFilters<std::uint8_t>::filteringByPruningMax(*weighted, attr64.data(), maxCriterion64, pruningMaxAttr64);
        requireVectorEqual(collectImageValues(pruningMinAttr64), collectImageValues(pruningMinAttrWeighted),
                           isMaxtree ? "weighted max-tree pruning-min double attr" : "weighted min-tree pruning-min double attr");
        requireVectorEqual(collectImageValues(pruningMaxAttr64), collectImageValues(pruningMaxAttrWeighted),
                           isMaxtree ? "weighted max-tree pruning-max double attr" : "weighted min-tree pruning-max double attr");

        if (isMaxtree) {
            requireThrows<std::invalid_argument>([&]() { UltimateAttributeOpening<std::uint8_t> invalidUao(*weighted, std::vector<float>{1.0f}); },
                                                 "UltimateAttributeOpening<std::uint8_t> must reject short vector attribute buffer");
        }

        {
            const auto [invalidParent, exportedAltitude] = weighted->exportHigraHierarchy();
            std::vector<int> invalidHigraAltitude;
            invalidHigraAltitude.reserve(exportedAltitude.size());
            const int invalidOffset = isMaxtree ? 300 : -300;
            for (std::uint8_t level : exportedAltitude) {
                invalidHigraAltitude.push_back(static_cast<int>(level) + invalidOffset);
            }
            auto invalidWeighted = MorphologicalTreeFactory::createFromHigraParent<int>(
                std::span<const NodeId>(invalidParent), std::span<const int>(invalidHigraAltitude), image->getNumRows(), image->getNumCols(),
                isMaxtree ? MorphologicalTreeKind::MAX_TREE : MorphologicalTreeKind::MIN_TREE,
                RegularGridAdjacency2D(image->getNumRows(), image->getNumCols(), 1.5));
            AttributeFilters<int> invalidObjectFilters(invalidWeighted);
            UltimateAttributeOpening<int> invalidObjectUao(invalidWeighted, attr);
            const WeightedTreeView<int> invalidExternalView = invalidWeighted.asView();
            std::vector<bool> invalidKeepAll(invalidWeighted.topology().getNumInternalNodeSlots(), true);
            auto invalidOutput =
                Image<int>::create(invalidWeighted.topology().getNumRowsOfGridDomain2D(), invalidWeighted.topology().getNumColsOfGridDomain2D(), 0);

            auto invalidReconstruction = invalidWeighted.reconstructionImage();
            requireImageShape(invalidReconstruction, image->getNumRows(), image->getNumCols());
            for (int value : collectImageValues(invalidReconstruction)) {
                require(isMaxtree ? value > 255 : value < 0,
                        isMaxtree ? "typed reconstruction must preserve altitude above 255" : "typed reconstruction must preserve negative altitude");
            }
            AttributeFilters<int>::filteringByDirectRule(invalidWeighted, invalidKeepAll, invalidOutput);
            requireVectorEqual(collectImageValues(invalidOutput), collectImageValues(invalidReconstruction),
                               isMaxtree ? "AttributeFilters<std::uint8_t> must preserve altitude above 255"
                                         : "AttributeFilters<std::uint8_t> must preserve negative altitude");
            auto invalidViewOutput =
                Image<int>::create(invalidWeighted.topology().getNumRowsOfGridDomain2D(), invalidWeighted.topology().getNumColsOfGridDomain2D(), 0);
            AttributeFilters<int>::filteringByDirectRule(invalidExternalView, invalidKeepAll, invalidViewOutput);
            requireVectorEqual(collectImageValues(invalidViewOutput), collectImageValues(invalidReconstruction),
                               isMaxtree ? "AttributeFilters<std::uint8_t> view must preserve altitude above 255"
                                         : "AttributeFilters<std::uint8_t> view must preserve negative altitude");
            requireVectorEqual(collectImageValues(invalidObjectFilters.filteringByDirectRule(invalidKeepAll)), collectImageValues(invalidReconstruction),
                               isMaxtree ? "AttributeFilters<std::uint8_t> object must observe replaced altitude above 255"
                                         : "AttributeFilters<std::uint8_t> object must observe replaced negative altitude");

            ExtinctionValues<int> invalidExtinction(invalidWeighted, attr);
            requireVectorEqual(collectImageValues(invalidExtinction.filtering(ExtinctionSelectionPolicy<float>::byTopK(1024))),
                               collectImageValues(invalidReconstruction),
                               isMaxtree ? "ExtinctionValues<std::uint8_t> must preserve altitude above 255"
                                         : "ExtinctionValues<std::uint8_t> must preserve negative altitude");

            UltimateAttributeOpening<int> invalidUao(invalidWeighted, attr);
            invalidUao.execute(static_cast<int>(maxCriterion));
            invalidObjectUao.execute(static_cast<int>(maxCriterion));
            auto invalidContrast = invalidUao.getMaxContrastImage();
            static_assert(std::is_same_v<decltype(invalidContrast), ImagePtr<int>>);
            requireImageShape(invalidContrast, image->getNumRows(), image->getNumCols());
            requireVectorEqual(collectImageValues(invalidObjectUao.getMaxContrastImage()), collectImageValues(invalidContrast),
                               isMaxtree ? "UltimateAttributeOpening<std::uint8_t> object must observe typed altitude above 255"
                                         : "UltimateAttributeOpening<std::uint8_t> object must observe typed negative altitude");
        }

        UltimateAttributeOpening<std::uint8_t> weightedUao(*weighted, attr);
        UltimateAttributeOpening<std::uint8_t> weightedUaoRaw(*weighted, attr.data());
        UltimateAttributeOpening<std::uint8_t> weightedUaoSelected(*weighted, attr);
        UltimateAttributeOpening<std::uint8_t> weightedUaoMser(*weighted, attr);
        UltimateAttributeOpening<std::uint8_t, double> weightedUao64(*weighted, attr64);
        weightedUao.execute(static_cast<int>(maxCriterion));
        weightedUaoRaw.execute(static_cast<int>(maxCriterion));
        weightedUao64.execute(maxCriterion64);
        std::vector<uint8_t> selectedAll(weighted->topology().getNumInternalNodeSlots(), true);
        weightedUaoSelected.execute(static_cast<int>(maxCriterion), selectedAll);
        weightedUaoMser.executeWithMSER(static_cast<int>(maxCriterion), 1);
        requireImageShape(weightedUao.getMaxContrastImage(), weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D());
        requireImageShape(weightedUao.getAssociatedImage(), weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D());
        requireImageShape(weightedUaoRaw.getMaxContrastImage(), weighted->topology().getNumRowsOfGridDomain2D(),
                          weighted->topology().getNumColsOfGridDomain2D());
        requireImageShape(weightedUaoRaw.getAssociatedImage(), weighted->topology().getNumRowsOfGridDomain2D(),
                          weighted->topology().getNumColsOfGridDomain2D());
        requireVectorEqual(collectImageValues(weightedUao64.getMaxContrastImage()), collectImageValues(weightedUao.getMaxContrastImage()),
                           isMaxtree ? "weighted max-tree UAO double attr max contrast" : "weighted min-tree UAO double attr max contrast");
        requireVectorEqual(collectImageValues(weightedUao64.getAssociatedImage()), collectImageValues(weightedUao.getAssociatedImage()),
                           isMaxtree ? "weighted max-tree UAO double attr associated image" : "weighted min-tree UAO double attr associated image");
        requireVectorEqual(collectImageValues(weightedUaoSelected.getMaxContrastImage()), collectImageValues(weightedUao.getMaxContrastImage()),
                           isMaxtree ? "weighted max-tree UAO selected-all max contrast" : "weighted min-tree UAO selected-all max contrast");
        requireVectorEqual(collectImageValues(weightedUaoSelected.getAssociatedImage()), collectImageValues(weightedUao.getAssociatedImage()),
                           isMaxtree ? "weighted max-tree UAO selected-all associated image" : "weighted min-tree UAO selected-all associated image");
        requireImageShape(weightedUaoMser.getMaxContrastImage(), weighted->topology().getNumRowsOfGridDomain2D(),
                          weighted->topology().getNumColsOfGridDomain2D());
        requireImageShape(weightedUaoMser.getAssociatedImage(), weighted->topology().getNumRowsOfGridDomain2D(),
                          weighted->topology().getNumColsOfGridDomain2D());

        ExtinctionValues<std::uint8_t> weightedExtinction(*weighted, attr);
        ExtinctionValues<std::uint8_t, double> weightedExtinction64(*weighted, attr64);
        ExtinctionValues<std::uint8_t> viewExtinction(externalView, attr);
        ExtinctionValues<std::int16_t> int16Extinction(int16View, attr);
        const auto keepAllFloat = ExtinctionSelectionPolicy<float>::byTopK(1024);
        const auto keepAllDouble = ExtinctionSelectionPolicy<double>::byTopK(1024);
        static_assert(std::is_same_v<decltype(weightedExtinction64.contourMap(keepAllDouble, ExtinctionContourScorePolicy::RankScore)), ImagePtr<double>>);
        requireVectorEqual(collectImageValues(weightedExtinction64.filtering(keepAllDouble)), collectImageValues(weightedExtinction.filtering(keepAllFloat)),
                           isMaxtree ? "weighted max-tree ExtinctionValues double attr filtering" : "weighted min-tree ExtinctionValues double attr filtering");
        requireVectorEqual(collectImageValuesAs<float>(weightedExtinction64.contourMap(keepAllDouble, ExtinctionContourScorePolicy::RankScore)),
                           collectImageValues(weightedExtinction.contourMap(keepAllFloat, ExtinctionContourScorePolicy::RankScore)),
                           isMaxtree ? "weighted max-tree ExtinctionValues double attr contour map" : "weighted min-tree ExtinctionValues double attr contour map");
        requireVectorEqual(collectImageValues(viewExtinction.filtering(keepAllFloat)), collectImageValues(weightedExtinction.filtering(keepAllFloat)),
                           isMaxtree ? "weighted max-tree ExtinctionValues<std::uint8_t> filtering via external view object"
                                     : "weighted min-tree ExtinctionValues<std::uint8_t> filtering via external view object");
        requireVectorEqual(collectImageValues(viewExtinction.contourMap(keepAllFloat, ExtinctionContourScorePolicy::RankScore)),
                           collectImageValues(weightedExtinction.contourMap(keepAllFloat, ExtinctionContourScorePolicy::RankScore)),
                           isMaxtree ? "weighted max-tree ExtinctionValues<std::uint8_t> contour map via external view object"
                                     : "weighted min-tree ExtinctionValues<std::uint8_t> contour map via external view object");
        requireVectorEqual(collectImageValues(int16Extinction.filtering(keepAllFloat)),
                           collectImageValuesAs<std::int16_t>(weightedExtinction.filtering(keepAllFloat)),
                           isMaxtree ? "weighted max-tree ExtinctionValues<std::uint8_t> filtering via int16 view object"
                                     : "weighted min-tree ExtinctionValues<std::uint8_t> filtering via int16 view object");
        requireVectorEqual(collectImageValues(int16Extinction.contourMap(keepAllFloat, ExtinctionContourScorePolicy::RankScore)),
                           collectImageValues(weightedExtinction.contourMap(keepAllFloat, ExtinctionContourScorePolicy::RankScore)),
                           isMaxtree ? "weighted max-tree ExtinctionValues<std::uint8_t> contour map via int16 view object"
                                     : "weighted min-tree ExtinctionValues<std::uint8_t> contour map via int16 view object");

        UltimateAttributeOpening<std::uint8_t> viewUao(externalView, attr);
        viewUao.execute(static_cast<int>(maxCriterion));
        requireVectorEqual(collectImageValues(viewUao.getMaxContrastImage()), collectImageValues(weightedUao.getMaxContrastImage()),
                           isMaxtree ? "weighted max-tree UAO max contrast via external view object"
                                     : "weighted min-tree UAO max contrast via external view object");
        requireVectorEqual(collectImageValues(viewUao.getAssociatedImage()), collectImageValues(weightedUao.getAssociatedImage()),
                           isMaxtree ? "weighted max-tree UAO associated via external view object"
                                     : "weighted min-tree UAO associated via external view object");
        UltimateAttributeOpening<std::int16_t> int16Uao(int16View, attr);
        int16Uao.execute(static_cast<int>(maxCriterion));
        requireVectorEqual(collectImageValues(int16Uao.getMaxContrastImage()), collectImageValuesAs<std::int16_t>(weightedUao.getMaxContrastImage()),
                           isMaxtree ? "weighted max-tree UAO max contrast via int16 view object" : "weighted min-tree UAO max contrast via int16 view object");
        requireVectorEqual(collectImageValues(int16Uao.getAssociatedImage()), collectImageValues(weightedUao.getAssociatedImage()),
                           isMaxtree ? "weighted max-tree UAO associated via int16 view object" : "weighted min-tree UAO associated via int16 view object");
        std::vector<float> floatAltitude;
        floatAltitude.reserve(externalAltitude.size());
        for (std::uint8_t level : externalAltitude) {
            floatAltitude.push_back(static_cast<float>(level) * 1.3f);
        }
        const WeightedTreeView<float> floatView(weighted->topology(), std::span<const float>(floatAltitude.data(), floatAltitude.size()));
        UltimateAttributeOpening<float> floatUao(floatView, attr);
        floatUao.execute(static_cast<int>(maxCriterion));
        auto floatContrast = floatUao.getMaxContrastImage();
        static_assert(std::is_same_v<decltype(floatContrast), ImageFloatPtr>);
        const auto canonicalContrast = collectImageValues(weightedUao.getMaxContrastImage());
        const auto floatContrastValues = collectImageValues(floatContrast);
        requireEqual(floatContrastValues.size(), canonicalContrast.size(), "weighted UAO float contrast size");
        bool hasNonZeroContrast = false;
        bool hasFractionalExpected = false;
        for (std::size_t i = 0; i < floatContrastValues.size(); ++i) {
            const float expected = static_cast<float>(canonicalContrast[i]) * 1.3f;
            requireNear(floatContrastValues[i], expected, 1.0e-5f, isMaxtree ? "weighted max-tree UAO float contrast" : "weighted min-tree UAO float contrast");
            hasNonZeroContrast = hasNonZeroContrast || canonicalContrast[i] != 0;
            hasFractionalExpected = hasFractionalExpected || std::abs(expected - std::round(expected)) > 1.0e-5f;
        }
        require(hasNonZeroContrast, "weighted UAO fixture must produce non-zero contrast");
        require(hasFractionalExpected, "weighted UAO float contrast test must exercise fractional output");
        requireThrows<std::logic_error>([&]() { viewUao.executeWithMSER(static_cast<int>(maxCriterion), 1); },
                                        isMaxtree ? "weighted max-tree UAO view object must reject MSER"
                                                  : "weighted min-tree UAO view object must reject MSER");
        if (isMaxtree) {
            requireThrows<std::invalid_argument>([&]() { weightedUao.execute(static_cast<int>(maxCriterion), std::vector<uint8_t>{true}); },
                                                 "UltimateAttributeOpening<std::uint8_t> must reject short selectedForFiltering buffer");
        }

        weighted->mergeNodeIntoParent(4);
        std::vector<bool> keepAllAfterMerge(weighted->topology().getNumInternalNodeSlots(), true);
        auto expectedAfterMerge = ImageUInt8::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        for (int pixel = 0; pixel < expectedAfterMerge->getSize(); ++pixel) {
            const NodeId nodeId = weighted->topology().getProperPartOwner(pixel);
            (*expectedAfterMerge)[pixel] = static_cast<uint8_t>(weighted->getAltitude(nodeId));
        }

        auto directAfterMerge = ImageUInt8::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        AttributeFilters<std::uint8_t>::filteringByDirectRule(*weighted, keepAllAfterMerge, directAfterMerge);
        requireVectorEqual(collectImageValues(directAfterMerge), collectImageValues(expectedAfterMerge),
                           isMaxtree ? "weighted max-tree direct-rule keep-all after merge" : "weighted min-tree direct-rule keep-all after merge");

        auto subtractiveAfterMerge = ImageUInt8::create(weighted->topology().getNumRowsOfGridDomain2D(), weighted->topology().getNumColsOfGridDomain2D(), 0);
        AttributeFilters<std::uint8_t>::filteringBySubtractiveRule(*weighted, keepAllAfterMerge, subtractiveAfterMerge);
        requireVectorEqual(collectImageValues(subtractiveAfterMerge), collectImageValues(expectedAfterMerge),
                           isMaxtree ? "weighted max-tree subtractive-rule keep-all after merge" : "weighted min-tree subtractive-rule keep-all after merge");

        auto [mergedAttrNames, mergedAttr] = AttributeComputation::computeSingleAttribute(*weighted, BOX_HEIGHT);
        (void)mergedAttrNames;
        UltimateAttributeOpening<std::uint8_t> mergedWeightedUao(*weighted, mergedAttr);
        UltimateAttributeOpening<std::uint8_t> mergedWeightedUaoRaw(*weighted, mergedAttr.data());
        mergedWeightedUao.execute(static_cast<int>(maxCriterion));
        mergedWeightedUaoRaw.execute(static_cast<int>(maxCriterion));
        requireImageShape(mergedWeightedUao.getMaxContrastImage(), weighted->topology().getNumRowsOfGridDomain2D(),
                          weighted->topology().getNumColsOfGridDomain2D());
        requireImageShape(mergedWeightedUao.getAssociatedImage(), weighted->topology().getNumRowsOfGridDomain2D(),
                          weighted->topology().getNumColsOfGridDomain2D());
        requireImageShape(mergedWeightedUaoRaw.getMaxContrastImage(), weighted->topology().getNumRowsOfGridDomain2D(),
                          weighted->topology().getNumColsOfGridDomain2D());
        requireImageShape(mergedWeightedUaoRaw.getAssociatedImage(), weighted->topology().getNumRowsOfGridDomain2D(),
                          weighted->topology().getNumColsOfGridDomain2D());
    }

    auto nonSquare = makeImage(2, 3, {3, 3, 2, 1, 4, 5});
    auto nonSquareWeighted = makeWeightedComponentTree(nonSquare, true);
    auto [nonSquareAttrNames, nonSquareAttr] = AttributeComputation::computeSingleAttribute(*nonSquareWeighted, BOX_HEIGHT);
    (void)nonSquareAttrNames;
    UltimateAttributeOpening<std::uint8_t> nonSquareUao(*nonSquareWeighted, nonSquareAttr);
    nonSquareUao.execute(nonSquareWeighted->topology().getNumRowsOfGridDomain2D());
    requireImageShape(nonSquareUao.getMaxContrastImage(), 2, 3);
    requireImageShape(nonSquareUao.getAssociatedImage(), 2, 3);
    requireImageShape(nonSquareUao.getAssociatedColorImage(), 2, 9);

    return 0;
}
