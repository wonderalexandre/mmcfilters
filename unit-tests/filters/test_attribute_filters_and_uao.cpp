#include "support/TestSupport.hpp"

#include "mmcfilters/attributes/AttributeComputation.hpp"
#include "mmcfilters/filters/AttributeFilters.hpp"
#include "mmcfilters/filters/detail/ViterbiDecision.hpp"
#include "mmcfilters/filters/ExtinctionValues.hpp"
#include "mmcfilters/filters/UltimateAttributeOpening.hpp"
#include "mmcfilters/trees/ValuedMorphologicalTreeView.hpp"

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
    for (PixelId pixel : tree.properPart(nodeId)) {
        (*image)[pixel] = value;
    }
}

template <class T>
void propagateDirectReference(const ValuedMorphologicalTree<T>& valuedTree, NodeId nodeId, const std::vector<bool>& preservationDecisions, std::vector<T>& mapLevel) {
    const MorphologicalTree& tree = valuedTree.topology();
    for (NodeId childNodeId : tree.children(nodeId)) {
        mapLevel[childNodeId] = preservationDecisions[childNodeId] ? valuedTree.nodeAltitude(childNodeId) : mapLevel[nodeId];
        propagateDirectReference(valuedTree, childNodeId, preservationDecisions, mapLevel);
    }
}

template <class T> std::vector<T> directReferenceImage(const ValuedMorphologicalTree<T>& valuedTree, const std::vector<bool>& preservationDecisions) {
    const MorphologicalTree& tree = valuedTree.topology();
    std::vector<T> mapLevel(static_cast<std::size_t>(tree.numInternalNodeSlots()), T{});
    const NodeId rootNodeId = tree.root();
    mapLevel[rootNodeId] = valuedTree.nodeAltitude(rootNodeId);
    propagateDirectReference(valuedTree, rootNodeId, preservationDecisions, mapLevel);

    auto image = Image<T>::create(tree.numRows(), tree.numColumns(), T{});
    for (NodeId nodeId : tree.aliveNodeIds()) {
        writeReferenceProperParts(tree, nodeId, image, mapLevel[nodeId]);
    }
    return collectImageValues(image);
}

template <class T>
void propagateSubtractiveReference(const ValuedMorphologicalTree<T>& valuedTree, NodeId nodeId, const std::vector<bool>& preservationDecisions,
                                   std::vector<AltitudeDifference<T>>& mapLevel) {
    const MorphologicalTree& tree = valuedTree.topology();
    for (NodeId childNodeId : tree.children(nodeId)) {
        mapLevel[childNodeId] =
            preservationDecisions[childNodeId] ? static_cast<AltitudeDifference<T>>(mapLevel[nodeId] + valuedTree.nodeResidue(childNodeId)) : mapLevel[nodeId];
        propagateSubtractiveReference(valuedTree, childNodeId, preservationDecisions, mapLevel);
    }
}

template <class T>
std::vector<AltitudeDifference<T>> subtractiveReferenceImage(const ValuedMorphologicalTree<T>& valuedTree, const std::vector<bool>& preservationDecisions) {
    const MorphologicalTree& tree = valuedTree.topology();
    std::vector<AltitudeDifference<T>> mapLevel(static_cast<std::size_t>(tree.numInternalNodeSlots()), AltitudeDifference<T>{});
    const NodeId rootNodeId = tree.root();
    mapLevel[rootNodeId] = preservationDecisions[static_cast<std::size_t>(rootNodeId)] ? valuedTree.nodeResidue(rootNodeId) : AltitudeDifference<T>{};
    propagateSubtractiveReference(valuedTree, rootNodeId, preservationDecisions, mapLevel);

    auto image = Image<AltitudeDifference<T>>::create(tree.numRows(), tree.numColumns(), AltitudeDifference<T>{});
    for (NodeId nodeId : tree.aliveNodeIds()) {
        writeReferenceProperParts(tree, nodeId, image, mapLevel[static_cast<std::size_t>(nodeId)]);
    }
    return collectImageValues(image);
}

ValuedMorphologicalTree<std::uint8_t> makeViterbiChainFixture() {
    const std::vector<NodeId> parent = {1, 2, 3, 3};
    const std::vector<std::uint8_t> altitude = {5, 5, 3, 0};
    return MorphologicalTreeFactory::createFromHigraParent<std::uint8_t>(std::span<const NodeId>(parent), std::span<const std::uint8_t>(altitude), 1, 1,
                                                                         MorphologicalTreeKind::MaxTree, RegularGridAdjacency2D(1, 1, 1.5));
}

ValuedMorphologicalTree<std::uint8_t> makeViterbiBranchFixture() {
    const std::vector<NodeId> parent = {3, 4, 4, 5, 5, 5};
    const std::vector<std::uint8_t> altitude = {5, 4, 4, 5, 4, 0};
    return MorphologicalTreeFactory::createFromHigraParent<std::uint8_t>(std::span<const NodeId>(parent), std::span<const std::uint8_t>(altitude), 1, 3,
                                                                         MorphologicalTreeKind::MaxTree, RegularGridAdjacency2D(1, 3, 1.5));
}

int main() {
    static_assert(!std::is_convertible_v<NodePreservationMask, NodePruningMask>);
    static_assert(!std::is_convertible_v<NodePruningMask, NodePreservationMask>);

    {
        auto chain = makeViterbiChainFixture();
        AttributeFilters<std::uint8_t> filters(chain);
        std::vector<float> attr = {1.1f, 0.0f, 2.0f};
        auto costs = detail::makeThresholdViterbiCosts(chain.topology(), attr.data(), 1.0f);
        auto keep = detail::computeViterbiPreservationDecisions(chain.topology(), costs);
        require(!keep[0] && !keep[1] && keep[2], "Viterbi chain must remove descendants below a removed ancestor");
        requireVectorEqual(collectImageValues(filters.filteringByViterbiRule(attr.data(), 1.0f)), std::vector<std::uint8_t>{0},
                           "Viterbi chain reconstruction with removed middle node");

        attr = {10.0f, 0.9f, 2.0f};
        keep = detail::computeViterbiPreservationDecisions(chain.topology(), detail::makeThresholdViterbiCosts(chain.topology(), attr.data(), 1.0f));
        require(keep[0] && keep[1] && keep[2], "Viterbi chain must preserve a costly-to-remove descendant branch");
        requireVectorEqual(collectImageValues(filters.filteringByViterbiRule(attr.data(), 1.0f)), std::vector<std::uint8_t>{5},
                           "Viterbi chain reconstruction with preserved branch");

        attr = {1.0f, 1.0f, 1.0f};
        costs = detail::makeThresholdViterbiCosts(chain.topology(), attr.data(), 1.0f);
        keep = detail::computeViterbiPreservationDecisions(chain.topology(), costs);
        require(!keep[0] && !keep[1] && keep[2], "Viterbi default tie-break must prefer removal except at the forced root");

        detail::ViterbiDecisionOptions preserveTies;
        preserveTies.tieBreak = detail::ViterbiTieBreak::PreferPreserve;
        keep = detail::computeViterbiPreservationDecisions(chain.topology(), costs, preserveTies);
        require(keep[0] && keep[1] && keep[2], "Viterbi preserve tie-break must keep the full chain");

        if constexpr (contract::validationsEnabled) {
            std::vector<float> nanAttr = {1.0f, std::numeric_limits<float>::quiet_NaN(), 1.0f};
            requireThrows<std::invalid_argument>([&]() { static_cast<void>(filters.filteringByViterbiRule(nanAttr.data(), 1.0f)); },
                                                 "Viterbi rule must reject NaN attributes");
            requireThrows<std::invalid_argument>(
                [&]() { static_cast<void>(filters.filteringByViterbiRule(attr.data(), std::numeric_limits<float>::quiet_NaN())); },
                "Viterbi rule must reject NaN thresholds");
            requireThrows<std::invalid_argument>(
                [&]() { static_cast<void>(filters.filteringByViterbiRule(static_cast<const float*>(nullptr), 1.0f)); },
                "Viterbi rule must reject null attribute buffers");
        }
    }

    {
        auto branch = makeViterbiBranchFixture();
        AttributeFilters<std::uint8_t> filters(branch);
        std::vector<double> attr = {2.0, 0.0, 2.0};
        auto keep = detail::computeViterbiPreservationDecisions(branch.topology(), detail::makeThresholdViterbiCosts(branch.topology(), attr.data(), 1.0));
        require(keep[0] && !keep[1] && keep[2], "Viterbi branch must decide sibling subtrees independently");
        requireVectorEqual(collectImageValues(filters.filteringByViterbiRule(attr.data(), 1.0)), std::vector<std::uint8_t>{5, 0, 0},
                           "Viterbi branch reconstruction");
    }

    auto image = makeComponentTreeFixture();

    for (bool isMaxtree : {true, false}) {
        auto valuedTree = makeValuedComponentTree(image, isMaxtree);
        auto reconstruction = valuedTree->reconstructFromNodeAltitudes();
        AttributeFilters<std::uint8_t> valuedTreeFilters(*valuedTree);
        const auto valuedTreeView = valuedTree->asView();

        std::vector<bool> keepAll(valuedTree->topology().numInternalNodeSlots(), true);
        const NodePreservationMask keepAllMask(keepAll);
        DirectAttributeFilter<std::uint8_t> directFilter(*valuedTree);
        HardSubtractiveAttributeFilter<std::uint8_t> hardSubtractiveFilter(*valuedTree);
        SoftSubtractiveAttributeFilter<std::uint8_t, float> softSubtractiveFilter(*valuedTree);

        const NodePruningMask pruneNoneMask = toNodePruningMask(keepAllMask);
        const NodePreservationMask roundTripMask = toNodePreservationMask(pruneNoneMask);
        for (std::size_t index = 0; index < keepAll.size(); ++index) {
            require(!pruneNoneMask[index] && roundTripMask[index], "node decision mask conversions must explicitly complement decisions");
        }

        if (isMaxtree) {
            const NodePreservationMask shortMask(std::vector<bool>(1, true));
            std::vector<float> shortScores(1, 0.0f);
            auto wrongShape = ImageUInt8::create(1, 1, 0);
            if constexpr (contract::validationsEnabled) {
                requireThrows<std::invalid_argument>([&]() { static_cast<void>(directFilter.applyDirectAttributeFilter(shortMask)); },
                                                     "DirectAttributeFilter must reject a short node preservation mask");
                requireThrows<std::invalid_argument>([&]() { static_cast<void>(softSubtractiveFilter.applySoftSubtractiveAttributeFilter(shortScores)); },
                                                     "SoftSubtractiveAttributeFilter must reject a short score buffer");

                NodePreservationMask rejectRootMask(keepAll);
                std::vector<bool> rejectRoot = keepAll;
                rejectRoot[static_cast<std::size_t>(valuedTree->topology().root())] = false;
                rejectRootMask = NodePreservationMask(std::move(rejectRoot));
                requireThrows<std::invalid_argument>([&]() { static_cast<void>(directFilter.applyDirectAttributeFilter(rejectRootMask)); },
                                                     "DirectAttributeFilter must reject a mask that rejects the root");

                std::vector<float> invalidScores(keepAll.size(), 1.0f);
                invalidScores.front() = std::numeric_limits<float>::quiet_NaN();
                requireThrows<std::invalid_argument>(
                    [&]() { static_cast<void>(softSubtractiveFilter.applySoftSubtractiveAttributeFilter(invalidScores)); },
                    "SoftSubtractiveAttributeFilter must reject NaN scores");
                invalidScores.front() = -0.01f;
                requireThrows<std::invalid_argument>(
                    [&]() { static_cast<void>(softSubtractiveFilter.applySoftSubtractiveAttributeFilter(invalidScores)); },
                    "SoftSubtractiveAttributeFilter must reject scores below zero");
                invalidScores.front() = 1.01f;
                requireThrows<std::invalid_argument>(
                    [&]() { static_cast<void>(softSubtractiveFilter.applySoftSubtractiveAttributeFilter(invalidScores)); },
                    "SoftSubtractiveAttributeFilter must reject scores above one");
                requireThrows<std::invalid_argument>(
                    [&]() { AttributeFilters<std::uint8_t>::filteringByPruningMin(*valuedTree, static_cast<const float*>(nullptr), 1.0f, wrongShape); },
                    "AttributeFilters<std::uint8_t> static pruning must reject null attribute pointer");
            }

            if constexpr (contract::validationsEnabled) {
                auto staleValuedTree = makeValuedComponentTree(image, true);
                auto [staleNames, staleAttr] = AttributeComputation::computeSingleAttribute(*staleValuedTree, MeanGrayLevel);
                (void)staleNames;
                DirectAttributeFilter<std::uint8_t> staleFilter(*staleValuedTree);
                UltimateAttributeOpening<std::uint8_t> staleUao(*staleValuedTree, staleAttr);
                NodePreservationMask staleKeepAll(static_cast<std::size_t>(staleValuedTree->topology().numInternalNodeSlots()), true);
                staleValuedTree->mergeNodeIntoParent(4);
                requireThrows<std::logic_error>([&]() { static_cast<void>(staleFilter.applyDirectAttributeFilter(staleKeepAll)); },
                                                "DirectAttributeFilter must reject use after topology mutation");
                requireThrows<std::logic_error>([&]() { staleUao.execute(4); },
                                                "UltimateAttributeOpening<std::uint8_t> execute must reject use after topology mutation");
                requireThrows<std::logic_error>([&]() { static_cast<void>(staleUao.getMaxContrastImage()); },
                                                "UltimateAttributeOpening<std::uint8_t> output must reject use after topology mutation");

                auto staleValuedTreeView = makeValuedComponentTree(image, true);
                const auto staleView = staleValuedTreeView->asView();
                DirectAttributeFilter<std::uint8_t> staleViewFilter(staleView);
                NodePreservationMask staleViewKeepAll(static_cast<std::size_t>(staleValuedTreeView->topology().numInternalNodeSlots()), true);
                staleValuedTreeView->mergeNodeIntoParent(4);
                requireThrows<std::logic_error>(
                    [&]() { static_cast<void>(staleViewFilter.applyDirectAttributeFilter(staleViewKeepAll)); },
                    "DirectAttributeFilter view API must reject stale ValuedMorphologicalTreeView");
                requireThrows<std::logic_error>([&]() { static_cast<void>(AttributeComputation::computeSingleAttribute(staleView, MeanGrayLevel)); },
                                                "AttributeComputation must reject stale ValuedMorphologicalTreeView");
            }
        }

        auto directViaObject = directFilter.applyDirectAttributeFilter(keepAllMask);
        requireVectorEqual(collectImageValues(directViaObject), collectImageValues(reconstruction),
                           isMaxtree ? "valuedTree max-tree direct filter keep-all via object" : "valuedTree min-tree direct filter keep-all via object");

        auto direct = applyDirectAttributeFilter(*valuedTree, keepAllMask);
        requireVectorEqual(collectImageValues(direct), collectImageValues(reconstruction),
                           isMaxtree ? "valuedTree max-tree direct filter keep-all" : "valuedTree min-tree direct filter keep-all");

        auto directViaView = applyDirectAttributeFilter(valuedTreeView, keepAllMask);
        requireVectorEqual(collectImageValues(directViaView), collectImageValues(reconstruction),
                           isMaxtree ? "valuedTree max-tree direct filter keep-all via view" : "valuedTree min-tree direct filter keep-all via view");

        NodeAltitudeBuffer<std::uint8_t> externalAltitude = valuedTree->nodeAltitudes();
        const ValuedMorphologicalTreeView<std::uint8_t> externalView(valuedTree->topology(), std::span<const std::uint8_t>(externalAltitude));
        std::vector<std::int16_t> int16Altitude;
        int16Altitude.reserve(externalAltitude.size());
        for (std::uint8_t level : externalAltitude) {
            int16Altitude.push_back(static_cast<std::int16_t>(level));
        }
        const ValuedMorphologicalTreeView<std::int16_t> int16View(valuedTree->topology(), std::span<const std::int16_t>(int16Altitude.data(), int16Altitude.size()));
        auto directViaExternalView = applyDirectAttributeFilter(externalView, keepAllMask);
        requireVectorEqual(collectImageValues(directViaExternalView), collectImageValues(reconstruction),
                           isMaxtree ? "valuedTree max-tree direct filter keep-all via external view" : "valuedTree min-tree direct filter keep-all via external view");

        DirectAttributeFilter<std::uint8_t> externalViewFilter(externalView);
        auto directViaExternalViewObject = externalViewFilter.applyDirectAttributeFilter(keepAllMask);
        requireVectorEqual(collectImageValues(directViaExternalViewObject), collectImageValues(reconstruction),
                           isMaxtree ? "valuedTree max-tree direct filter keep-all via external view object"
                                     : "valuedTree min-tree direct filter keep-all via external view object");

        auto directViaInt16View = applyDirectAttributeFilter(int16View, keepAllMask);
        requireVectorEqual(collectImageValues(directViaInt16View), collectImageValuesAs<std::int16_t>(reconstruction),
                           isMaxtree ? "valuedTree max-tree direct filter keep-all via int16 view" : "valuedTree min-tree direct filter keep-all via int16 view");

        DirectAttributeFilter<std::int16_t> int16ViewFilter(int16View);
        auto directViaInt16ViewObject = int16ViewFilter.applyDirectAttributeFilter(keepAllMask);
        requireVectorEqual(collectImageValues(directViaInt16ViewObject), collectImageValuesAs<std::int16_t>(reconstruction),
                           isMaxtree ? "valuedTree max-tree direct filter keep-all via int16 view object"
                                     : "valuedTree min-tree direct filter keep-all via int16 view object");

        std::vector<float> unitScores(valuedTree->topology().numInternalNodeSlots(), 1.0f);
        auto scoreViaObject = softSubtractiveFilter.applySoftSubtractiveAttributeFilter(unitScores);
        auto scoreViaView = applySoftSubtractiveAttributeFilter(valuedTreeView, std::span<const float>(unitScores));
        auto scoreViaExternalView = applySoftSubtractiveAttributeFilter(externalView, std::span<const float>(unitScores));
        auto scoreViaInt16View = applySoftSubtractiveAttributeFilter(int16View, std::span<const float>(unitScores));
        requireVectorEqual(collectImageValues(scoreViaView), collectImageValues(scoreViaObject),
                           isMaxtree ? "valuedTree max-tree soft subtractive filter via view" : "valuedTree min-tree soft subtractive filter via view");
        requireVectorEqual(collectImageValues(scoreViaExternalView), collectImageValues(scoreViaObject),
                           isMaxtree ? "valuedTree max-tree soft subtractive filter via external view" : "valuedTree min-tree soft subtractive filter via external view");
        requireVectorEqual(collectImageValues(scoreViaInt16View), collectImageValues(scoreViaObject),
                           isMaxtree ? "valuedTree max-tree soft subtractive filter via int16 view" : "valuedTree min-tree soft subtractive filter via int16 view");

        const auto hardKeepAll = hardSubtractiveFilter.applyHardSubtractiveAttributeFilter(keepAllMask);
        requireVectorEqual(collectImageValues(hardKeepAll), collectImageValuesAs<AltitudeDifference<std::uint8_t>>(reconstruction),
                           isMaxtree ? "valuedTree max-tree hard subtractive filter keep-all"
                                     : "valuedTree min-tree hard subtractive filter keep-all");

        const NodePreservationMask rejectAllMask(keepAll.size(), false);
        const auto hardRejectAll = hardSubtractiveFilter.applyHardSubtractiveAttributeFilter(rejectAllMask);
        requireVectorEqual(collectImageValues(hardRejectAll),
                           std::vector<AltitudeDifference<std::uint8_t>>(static_cast<std::size_t>(hardRejectAll->getSize()), 0),
                           isMaxtree ? "valuedTree max-tree all-false hard mask must produce zero"
                                     : "valuedTree min-tree all-false hard mask must produce zero");

        const std::vector<float> zeroScores(keepAll.size(), 0.0f);
        const auto softRejectAll = softSubtractiveFilter.applySoftSubtractiveAttributeFilter(zeroScores);
        requireVectorEqual(collectImageValues(softRejectAll), std::vector<float>(static_cast<std::size_t>(softRejectAll->getSize()), 0.0f),
                           isMaxtree ? "valuedTree max-tree all-zero soft scores must produce zero"
                                     : "valuedTree min-tree all-zero soft scores must produce zero");

        std::vector<AltitudeDifference<std::uint8_t>> allResidues(keepAll.size(), 0);
        for (NodeId nodeId : valuedTree->topology().aliveNodeIds()) {
            allResidues[static_cast<std::size_t>(nodeId)] = valuedTree->nodeResidue(nodeId);
        }
        requireVectorEqual(collectImageValues(valuedTree->reconstructFromNodeContributions(std::span<const AltitudeDifference<std::uint8_t>>(allResidues))),
                           collectImageValues(hardKeepAll),
                           isMaxtree ? "valuedTree max-tree general contribution reconstruction"
                                     : "valuedTree min-tree general contribution reconstruction");

        {
            const auto [higraParent, higraAltitude] = valuedTree->exportHigraHierarchy();
            auto importedValuedTree = MorphologicalTreeFactory::createFromHigraParent<std::uint8_t>(
                std::span<const NodeId>(higraParent), std::span<const std::uint8_t>(higraAltitude), image->getNumRows(), image->getNumColumns(),
                isMaxtree ? MorphologicalTreeKind::MaxTree : MorphologicalTreeKind::MinTree,
                RegularGridAdjacency2D(image->getNumRows(), image->getNumColumns(), 1.5));
            std::vector<bool> importedKeepAll(importedValuedTree.topology().numInternalNodeSlots(), true);
            const NodePreservationMask importedKeepAllMask(importedKeepAll);
            auto importedReconstruction = importedValuedTree.reconstructFromNodeAltitudes();

            requireVectorEqual(collectImageValues(applyHardSubtractiveAttributeFilter(importedValuedTree, importedKeepAllMask)),
                               collectImageValuesAs<AltitudeDifference<std::uint8_t>>(importedReconstruction),
                               isMaxtree ? "imported max-tree hard subtractive filter keep-all"
                                         : "imported min-tree hard subtractive filter keep-all");

            std::vector<float> importedUnitScores(importedValuedTree.topology().numInternalNodeSlots(), 1.0f);
            requireVectorEqual(collectImageValues(applySoftSubtractiveAttributeFilter(importedValuedTree, std::span<const float>(importedUnitScores))),
                               collectImageValuesAs<float>(importedReconstruction),
                               isMaxtree ? "imported max-tree soft subtractive unit scores" : "imported min-tree soft subtractive unit scores");

            std::vector<bool> mixedPreservationDecisions(importedValuedTree.topology().numInternalNodeSlots(), true);
            int rejectedNonRootNodes = 0;
            for (NodeId nodeId : importedValuedTree.topology().aliveNodeIds()) {
                if (nodeId != importedValuedTree.topology().root() && rejectedNonRootNodes < 2) {
                    mixedPreservationDecisions[static_cast<std::size_t>(nodeId)] = false;
                    ++rejectedNonRootNodes;
                }
            }
            const NodePreservationMask mixedMask(mixedPreservationDecisions);

            requireVectorEqual(collectImageValues(applyDirectAttributeFilter(importedValuedTree, mixedMask)),
                               directReferenceImage(importedValuedTree, mixedPreservationDecisions),
                               isMaxtree ? "imported max-tree direct filter mixed mask" : "imported min-tree direct filter mixed mask");
            requireVectorEqual(collectImageValues(applyHardSubtractiveAttributeFilter(importedValuedTree, mixedMask)),
                               subtractiveReferenceImage(importedValuedTree, mixedPreservationDecisions),
                               isMaxtree ? "imported max-tree hard subtractive filter mixed mask"
                                         : "imported min-tree hard subtractive filter mixed mask");
        }

        requireVectorEqual(collectImageValues(applyHardSubtractiveAttributeFilter(valuedTreeView, keepAllMask)), collectImageValues(hardKeepAll),
                           isMaxtree ? "valuedTree max-tree hard subtractive filter via view"
                                     : "valuedTree min-tree hard subtractive filter via view");
        requireVectorEqual(collectImageValues(applyHardSubtractiveAttributeFilter(externalView, keepAllMask)), collectImageValues(hardKeepAll),
                           isMaxtree ? "valuedTree max-tree hard subtractive filter via external view"
                                     : "valuedTree min-tree hard subtractive filter via external view");
        requireVectorEqual(collectImageValues(applyHardSubtractiveAttributeFilter(int16View, keepAllMask)),
                           collectImageValuesAs<AltitudeDifference<std::int16_t>>(hardKeepAll),
                           isMaxtree ? "valuedTree max-tree hard subtractive filter via int16 view"
                                     : "valuedTree min-tree hard subtractive filter via int16 view");

        auto pruningMin = ImageUInt8::create(valuedTree->topology().numRows(), valuedTree->topology().numColumns(), 0);
        AttributeFilters<std::uint8_t>::filteringByPruningMin(*valuedTree, keepAllMask, pruningMin);
        requireVectorEqual(collectImageValues(pruningMin), collectImageValues(reconstruction),
                           isMaxtree ? "valuedTree max-tree pruning-min keep-all" : "valuedTree min-tree pruning-min keep-all");

        auto pruningMinViaView = ImageUInt8::create(valuedTree->topology().numRows(), valuedTree->topology().numColumns(), 0);
        auto pruningMinViaExternalView =
            ImageUInt8::create(valuedTree->topology().numRows(), valuedTree->topology().numColumns(), 0);
        auto pruningMinViaInt16View =
            Image<std::int16_t>::create(valuedTree->topology().numRows(), valuedTree->topology().numColumns(), 0);
        AttributeFilters<std::uint8_t>::filteringByPruningMin(valuedTreeView, keepAllMask, pruningMinViaView);
        AttributeFilters<std::uint8_t>::filteringByPruningMin(externalView, keepAllMask, pruningMinViaExternalView);
        AttributeFilters<std::int16_t>::filteringByPruningMin(int16View, keepAllMask, pruningMinViaInt16View);
        requireVectorEqual(collectImageValues(pruningMinViaView), collectImageValues(pruningMin),
                           isMaxtree ? "valuedTree max-tree pruning-min via view" : "valuedTree min-tree pruning-min via view");
        requireVectorEqual(collectImageValues(pruningMinViaExternalView), collectImageValues(pruningMin),
                           isMaxtree ? "valuedTree max-tree pruning-min via external view" : "valuedTree min-tree pruning-min via external view");
        requireVectorEqual(collectImageValues(pruningMinViaInt16View), collectImageValuesAs<std::int16_t>(pruningMin),
                           isMaxtree ? "valuedTree max-tree pruning-min via int16 view" : "valuedTree min-tree pruning-min via int16 view");

        auto pruningMax = ImageUInt8::create(valuedTree->topology().numRows(), valuedTree->topology().numColumns(), 0);
        AttributeFilters<std::uint8_t>::filteringByPruningMax(*valuedTree, keepAllMask, pruningMax);
        requireVectorEqual(collectImageValues(pruningMax), collectImageValues(reconstruction),
                           isMaxtree ? "valuedTree max-tree pruning-max keep-all" : "valuedTree min-tree pruning-max keep-all");

        auto pruningMaxViaView = ImageUInt8::create(valuedTree->topology().numRows(), valuedTree->topology().numColumns(), 0);
        auto pruningMaxViaExternalView =
            ImageUInt8::create(valuedTree->topology().numRows(), valuedTree->topology().numColumns(), 0);
        auto pruningMaxViaInt16View =
            Image<std::int16_t>::create(valuedTree->topology().numRows(), valuedTree->topology().numColumns(), 0);
        AttributeFilters<std::uint8_t>::filteringByPruningMax(valuedTreeView, keepAllMask, pruningMaxViaView);
        AttributeFilters<std::uint8_t>::filteringByPruningMax(externalView, keepAllMask, pruningMaxViaExternalView);
        AttributeFilters<std::int16_t>::filteringByPruningMax(int16View, keepAllMask, pruningMaxViaInt16View);
        requireVectorEqual(collectImageValues(pruningMaxViaView), collectImageValues(pruningMax),
                           isMaxtree ? "valuedTree max-tree pruning-max via view" : "valuedTree min-tree pruning-max via view");
        requireVectorEqual(collectImageValues(pruningMaxViaExternalView), collectImageValues(pruningMax),
                           isMaxtree ? "valuedTree max-tree pruning-max via external view" : "valuedTree min-tree pruning-max via external view");
        requireVectorEqual(collectImageValues(pruningMaxViaInt16View), collectImageValuesAs<std::int16_t>(pruningMax),
                           isMaxtree ? "valuedTree max-tree pruning-max via int16 view" : "valuedTree min-tree pruning-max via int16 view");

        auto [attrNames, attr] = AttributeComputation::computeSingleAttribute(*valuedTree, BoundingBoxHeight);
        (void)attrNames;
        auto [attrNames64, attr64] = AttributeComputation::computeSingleAttribute<double>(*valuedTree, BoundingBoxHeight);
        (void)attrNames64;
        float maximumAttributeThreshold = static_cast<float>(valuedTree->topology().numRows());
        double maximumAttributeThreshold64 = static_cast<double>(valuedTree->topology().numRows());

        for (float threshold : {1.0f, 2.0f, 3.0f, 4.0f}) {
            std::vector<bool> keepByAttribute(valuedTree->topology().numInternalNodeSlots(), false);
            for (NodeId nodeId : valuedTree->topology().aliveNodeIds()) {
                keepByAttribute[nodeId] = attr[static_cast<std::size_t>(nodeId)] > threshold;
            }
            const NodePreservationMask keepByAttributeMask(keepByAttribute);

            requireVectorEqual(collectImageValues(valuedTreeFilters.filteringByPruningMin(keepByAttributeMask)),
                               collectImageValues(valuedTreeFilters.filteringByPruningMin(attr.data(), threshold)),
                               isMaxtree ? "valuedTree max-tree pruning-min preservation-mask/attribute-threshold equivalence"
                                         : "valuedTree min-tree pruning-min preservation-mask/attribute-threshold equivalence");
            requireVectorEqual(collectImageValues(valuedTreeFilters.filteringByPruningMax(keepByAttributeMask)),
                               collectImageValues(valuedTreeFilters.filteringByPruningMax(attr.data(), threshold)),
                               isMaxtree ? "valuedTree max-tree pruning-max preservation-mask/attribute-threshold equivalence"
                                         : "valuedTree min-tree pruning-max preservation-mask/attribute-threshold equivalence");
        }

        auto pruningMinAttributeValuedTree = ImageUInt8::create(valuedTree->topology().numRows(), valuedTree->topology().numColumns(), 0);
        auto pruningMinAttrView = ImageUInt8::create(valuedTree->topology().numRows(), valuedTree->topology().numColumns(), 0);
        auto pruningMinAttrExternalView =
            ImageUInt8::create(valuedTree->topology().numRows(), valuedTree->topology().numColumns(), 0);
        auto pruningMinAttrInt16View =
            Image<std::int16_t>::create(valuedTree->topology().numRows(), valuedTree->topology().numColumns(), 0);
        AttributeFilters<std::uint8_t>::filteringByPruningMin(*valuedTree, attr.data(), maximumAttributeThreshold, pruningMinAttributeValuedTree);
        AttributeFilters<std::uint8_t>::filteringByPruningMin(valuedTreeView, attr.data(), maximumAttributeThreshold, pruningMinAttrView);
        AttributeFilters<std::uint8_t>::filteringByPruningMin(externalView, attr.data(), maximumAttributeThreshold, pruningMinAttrExternalView);
        AttributeFilters<std::int16_t>::filteringByPruningMin(int16View, attr.data(), maximumAttributeThreshold, pruningMinAttrInt16View);
        requireVectorEqual(collectImageValues(pruningMinAttrView), collectImageValues(pruningMinAttributeValuedTree),
                           isMaxtree ? "valuedTree max-tree pruning-min attr via view" : "valuedTree min-tree pruning-min attr via view");
        requireVectorEqual(collectImageValues(pruningMinAttrExternalView), collectImageValues(pruningMinAttributeValuedTree),
                           isMaxtree ? "valuedTree max-tree pruning-min attr via external view" : "valuedTree min-tree pruning-min attr via external view");
        requireVectorEqual(collectImageValues(pruningMinAttrInt16View), collectImageValuesAs<std::int16_t>(pruningMinAttributeValuedTree),
                           isMaxtree ? "valuedTree max-tree pruning-min attr via int16 view" : "valuedTree min-tree pruning-min attr via int16 view");

        auto pruningMaxAttributeValuedTree = ImageUInt8::create(valuedTree->topology().numRows(), valuedTree->topology().numColumns(), 0);
        auto pruningMaxAttrView = ImageUInt8::create(valuedTree->topology().numRows(), valuedTree->topology().numColumns(), 0);
        auto pruningMaxAttrExternalView =
            ImageUInt8::create(valuedTree->topology().numRows(), valuedTree->topology().numColumns(), 0);
        auto pruningMaxAttrInt16View =
            Image<std::int16_t>::create(valuedTree->topology().numRows(), valuedTree->topology().numColumns(), 0);
        AttributeFilters<std::uint8_t>::filteringByPruningMax(*valuedTree, attr.data(), maximumAttributeThreshold, pruningMaxAttributeValuedTree);
        AttributeFilters<std::uint8_t>::filteringByPruningMax(valuedTreeView, attr.data(), maximumAttributeThreshold, pruningMaxAttrView);
        AttributeFilters<std::uint8_t>::filteringByPruningMax(externalView, attr.data(), maximumAttributeThreshold, pruningMaxAttrExternalView);
        AttributeFilters<std::int16_t>::filteringByPruningMax(int16View, attr.data(), maximumAttributeThreshold, pruningMaxAttrInt16View);
        requireVectorEqual(collectImageValues(pruningMaxAttrView), collectImageValues(pruningMaxAttributeValuedTree),
                           isMaxtree ? "valuedTree max-tree pruning-max attr via view" : "valuedTree min-tree pruning-max attr via view");
        requireVectorEqual(collectImageValues(pruningMaxAttrExternalView), collectImageValues(pruningMaxAttributeValuedTree),
                           isMaxtree ? "valuedTree max-tree pruning-max attr via external view" : "valuedTree min-tree pruning-max attr via external view");
        requireVectorEqual(collectImageValues(pruningMaxAttrInt16View), collectImageValuesAs<std::int16_t>(pruningMaxAttributeValuedTree),
                           isMaxtree ? "valuedTree max-tree pruning-max attr via int16 view" : "valuedTree min-tree pruning-max attr via int16 view");

        auto pruningMinAttr64 = ImageUInt8::create(valuedTree->topology().numRows(), valuedTree->topology().numColumns(), 0);
        auto pruningMaxAttr64 = ImageUInt8::create(valuedTree->topology().numRows(), valuedTree->topology().numColumns(), 0);
        AttributeFilters<std::uint8_t>::filteringByPruningMin(*valuedTree, attr64.data(), maximumAttributeThreshold64, pruningMinAttr64);
        AttributeFilters<std::uint8_t>::filteringByPruningMax(*valuedTree, attr64.data(), maximumAttributeThreshold64, pruningMaxAttr64);
        requireVectorEqual(collectImageValues(pruningMinAttr64), collectImageValues(pruningMinAttributeValuedTree),
                           isMaxtree ? "valuedTree max-tree pruning-min double attr" : "valuedTree min-tree pruning-min double attr");
        requireVectorEqual(collectImageValues(pruningMaxAttr64), collectImageValues(pruningMaxAttributeValuedTree),
                           isMaxtree ? "valuedTree max-tree pruning-max double attr" : "valuedTree min-tree pruning-max double attr");

        if (isMaxtree && contract::validationsEnabled) {
            requireThrows<std::invalid_argument>([&]() { UltimateAttributeOpening<std::uint8_t> invalidUao(*valuedTree, std::vector<float>{1.0f}); },
                                                 "UltimateAttributeOpening<std::uint8_t> must reject short vector attribute buffer");
        }

        {
            const auto [invalidParent, exportedAltitude] = valuedTree->exportHigraHierarchy();
            std::vector<int> invalidHigraAltitude;
            invalidHigraAltitude.reserve(exportedAltitude.size());
            const int invalidOffset = isMaxtree ? 300 : -300;
            for (std::uint8_t level : exportedAltitude) {
                invalidHigraAltitude.push_back(static_cast<int>(level) + invalidOffset);
            }
            auto invalidValuedTree = MorphologicalTreeFactory::createFromHigraParent<int>(
                std::span<const NodeId>(invalidParent), std::span<const int>(invalidHigraAltitude), image->getNumRows(), image->getNumColumns(),
                isMaxtree ? MorphologicalTreeKind::MaxTree : MorphologicalTreeKind::MinTree,
                RegularGridAdjacency2D(image->getNumRows(), image->getNumColumns(), 1.5));
            UltimateAttributeOpening<int> invalidObjectUao(invalidValuedTree, attr);
            const ValuedMorphologicalTreeView<int> invalidExternalView = invalidValuedTree.asView();
            std::vector<bool> invalidKeepAll(invalidValuedTree.topology().numInternalNodeSlots(), true);
            const NodePreservationMask invalidKeepAllMask(invalidKeepAll);

            auto invalidReconstruction = invalidValuedTree.reconstructFromNodeAltitudes();
            requireImageShape(invalidReconstruction, image->getNumRows(), image->getNumColumns());
            for (int value : collectImageValues(invalidReconstruction)) {
                require(isMaxtree ? value > 255 : value < 0,
                        isMaxtree ? "typed reconstruction must preserve altitude above 255" : "typed reconstruction must preserve negative altitude");
            }
            const auto invalidOutput = applyDirectAttributeFilter(invalidValuedTree, invalidKeepAllMask);
            requireVectorEqual(collectImageValues(invalidOutput), collectImageValues(invalidReconstruction),
                               isMaxtree ? "AttributeFilters<std::uint8_t> must preserve altitude above 255"
                                         : "AttributeFilters<std::uint8_t> must preserve negative altitude");
            const auto invalidViewOutput = applyDirectAttributeFilter(invalidExternalView, invalidKeepAllMask);
            requireVectorEqual(collectImageValues(invalidViewOutput), collectImageValues(invalidReconstruction),
                               isMaxtree ? "AttributeFilters<std::uint8_t> view must preserve altitude above 255"
                                         : "AttributeFilters<std::uint8_t> view must preserve negative altitude");
            requireVectorEqual(collectImageValues(DirectAttributeFilter<int>(invalidValuedTree).applyDirectAttributeFilter(invalidKeepAllMask)),
                               collectImageValues(invalidReconstruction),
                               isMaxtree ? "AttributeFilters<std::uint8_t> object must observe replaced altitude above 255"
                                         : "AttributeFilters<std::uint8_t> object must observe replaced negative altitude");

            ExtinctionValues<int> invalidExtinction(invalidValuedTree, attr);
            requireVectorEqual(collectImageValues(invalidExtinction.filtering(ExtinctionSelectionPolicy<float>::byTopK(1024))),
                               collectImageValues(invalidReconstruction),
                               isMaxtree ? "ExtinctionValues<std::uint8_t> must preserve altitude above 255"
                                         : "ExtinctionValues<std::uint8_t> must preserve negative altitude");

            UltimateAttributeOpening<int> invalidUao(invalidValuedTree, attr);
            invalidUao.execute(static_cast<int>(maximumAttributeThreshold));
            invalidObjectUao.execute(static_cast<int>(maximumAttributeThreshold));
            auto invalidContrast = invalidUao.getMaxContrastImage();
            static_assert(std::is_same_v<decltype(invalidContrast), ImagePtr<int>>);
            requireImageShape(invalidContrast, image->getNumRows(), image->getNumColumns());
            requireVectorEqual(collectImageValues(invalidObjectUao.getMaxContrastImage()), collectImageValues(invalidContrast),
                               isMaxtree ? "UltimateAttributeOpening<std::uint8_t> object must observe typed altitude above 255"
                                         : "UltimateAttributeOpening<std::uint8_t> object must observe typed negative altitude");
        }

        UltimateAttributeOpening<std::uint8_t> valuedTreeUao(*valuedTree, attr);
        UltimateAttributeOpening<std::uint8_t> valuedTreeUaoRaw(*valuedTree, attr.data());
        UltimateAttributeOpening<std::uint8_t> valuedTreeUaoSelected(*valuedTree, attr);
        UltimateAttributeOpening<std::uint8_t> valuedTreeUaoMser(*valuedTree, attr);
        UltimateAttributeOpening<std::uint8_t, double> valuedTreeUao64(*valuedTree, attr64);
        valuedTreeUao.execute(static_cast<int>(maximumAttributeThreshold));
        valuedTreeUaoRaw.execute(static_cast<int>(maximumAttributeThreshold));
        valuedTreeUao64.execute(maximumAttributeThreshold64);
        std::vector<uint8_t> selectedAll(valuedTree->topology().numInternalNodeSlots(), true);
        valuedTreeUaoSelected.execute(static_cast<int>(maximumAttributeThreshold), selectedAll);
        valuedTreeUaoMser.executeWithMSER(static_cast<int>(maximumAttributeThreshold), 1);
        requireImageShape(valuedTreeUao.getMaxContrastImage(), valuedTree->topology().numRows(), valuedTree->topology().numColumns());
        requireImageShape(valuedTreeUao.getAssociatedImage(), valuedTree->topology().numRows(), valuedTree->topology().numColumns());
        requireImageShape(valuedTreeUaoRaw.getMaxContrastImage(), valuedTree->topology().numRows(),
                          valuedTree->topology().numColumns());
        requireImageShape(valuedTreeUaoRaw.getAssociatedImage(), valuedTree->topology().numRows(),
                          valuedTree->topology().numColumns());
        requireVectorEqual(collectImageValues(valuedTreeUao64.getMaxContrastImage()), collectImageValues(valuedTreeUao.getMaxContrastImage()),
                           isMaxtree ? "valuedTree max-tree UAO double attr max contrast" : "valuedTree min-tree UAO double attr max contrast");
        requireVectorEqual(collectImageValues(valuedTreeUao64.getAssociatedImage()), collectImageValues(valuedTreeUao.getAssociatedImage()),
                           isMaxtree ? "valuedTree max-tree UAO double attr associated image" : "valuedTree min-tree UAO double attr associated image");
        requireVectorEqual(collectImageValues(valuedTreeUaoSelected.getMaxContrastImage()), collectImageValues(valuedTreeUao.getMaxContrastImage()),
                           isMaxtree ? "valuedTree max-tree UAO selected-all max contrast" : "valuedTree min-tree UAO selected-all max contrast");
        requireVectorEqual(collectImageValues(valuedTreeUaoSelected.getAssociatedImage()), collectImageValues(valuedTreeUao.getAssociatedImage()),
                           isMaxtree ? "valuedTree max-tree UAO selected-all associated image" : "valuedTree min-tree UAO selected-all associated image");
        requireImageShape(valuedTreeUaoMser.getMaxContrastImage(), valuedTree->topology().numRows(),
                          valuedTree->topology().numColumns());
        requireImageShape(valuedTreeUaoMser.getAssociatedImage(), valuedTree->topology().numRows(),
                          valuedTree->topology().numColumns());

        ExtinctionValues<std::uint8_t> valuedTreeExtinction(*valuedTree, attr);
        ExtinctionValues<std::uint8_t, double> valuedTreeExtinction64(*valuedTree, attr64);
        ExtinctionValues<std::uint8_t> viewExtinction(externalView, attr);
        ExtinctionValues<std::int16_t> int16Extinction(int16View, attr);
        const auto keepAllFloat = ExtinctionSelectionPolicy<float>::byTopK(1024);
        const auto keepAllDouble = ExtinctionSelectionPolicy<double>::byTopK(1024);
        static_assert(std::is_same_v<decltype(valuedTreeExtinction64.contourMap(keepAllDouble, ExtinctionContourScorePolicy::RankScore)), ImagePtr<double>>);
        requireVectorEqual(collectImageValues(valuedTreeExtinction64.filtering(keepAllDouble)), collectImageValues(valuedTreeExtinction.filtering(keepAllFloat)),
                           isMaxtree ? "valuedTree max-tree ExtinctionValues double attr filtering" : "valuedTree min-tree ExtinctionValues double attr filtering");
        requireVectorEqual(collectImageValuesAs<float>(valuedTreeExtinction64.contourMap(keepAllDouble, ExtinctionContourScorePolicy::RankScore)),
                           collectImageValues(valuedTreeExtinction.contourMap(keepAllFloat, ExtinctionContourScorePolicy::RankScore)),
                           isMaxtree ? "valuedTree max-tree ExtinctionValues double attr saliency" : "valuedTree min-tree ExtinctionValues double attr saliency");
        requireVectorEqual(collectImageValues(viewExtinction.filtering(keepAllFloat)), collectImageValues(valuedTreeExtinction.filtering(keepAllFloat)),
                           isMaxtree ? "valuedTree max-tree ExtinctionValues<std::uint8_t> filtering via external view object"
                                     : "valuedTree min-tree ExtinctionValues<std::uint8_t> filtering via external view object");
        requireVectorEqual(collectImageValues(viewExtinction.contourMap(keepAllFloat, ExtinctionContourScorePolicy::RankScore)),
                           collectImageValues(valuedTreeExtinction.contourMap(keepAllFloat, ExtinctionContourScorePolicy::RankScore)),
                           isMaxtree ? "valuedTree max-tree ExtinctionValues<std::uint8_t> saliency via external view object"
                                     : "valuedTree min-tree ExtinctionValues<std::uint8_t> saliency via external view object");
        requireVectorEqual(collectImageValues(int16Extinction.filtering(keepAllFloat)),
                           collectImageValuesAs<std::int16_t>(valuedTreeExtinction.filtering(keepAllFloat)),
                           isMaxtree ? "valuedTree max-tree ExtinctionValues<std::uint8_t> filtering via int16 view object"
                                     : "valuedTree min-tree ExtinctionValues<std::uint8_t> filtering via int16 view object");
        requireVectorEqual(collectImageValues(int16Extinction.contourMap(keepAllFloat, ExtinctionContourScorePolicy::RankScore)),
                           collectImageValues(valuedTreeExtinction.contourMap(keepAllFloat, ExtinctionContourScorePolicy::RankScore)),
                           isMaxtree ? "valuedTree max-tree ExtinctionValues<std::uint8_t> saliency via int16 view object"
                                     : "valuedTree min-tree ExtinctionValues<std::uint8_t> saliency via int16 view object");

        UltimateAttributeOpening<std::uint8_t> viewUao(externalView, attr);
        viewUao.execute(static_cast<int>(maximumAttributeThreshold));
        requireVectorEqual(collectImageValues(viewUao.getMaxContrastImage()), collectImageValues(valuedTreeUao.getMaxContrastImage()),
                           isMaxtree ? "valuedTree max-tree UAO max contrast via external view object"
                                     : "valuedTree min-tree UAO max contrast via external view object");
        requireVectorEqual(collectImageValues(viewUao.getAssociatedImage()), collectImageValues(valuedTreeUao.getAssociatedImage()),
                           isMaxtree ? "valuedTree max-tree UAO associated via external view object"
                                     : "valuedTree min-tree UAO associated via external view object");
        UltimateAttributeOpening<std::int16_t> int16Uao(int16View, attr);
        int16Uao.execute(static_cast<int>(maximumAttributeThreshold));
        requireVectorEqual(collectImageValues(int16Uao.getMaxContrastImage()), collectImageValuesAs<std::int16_t>(valuedTreeUao.getMaxContrastImage()),
                           isMaxtree ? "valuedTree max-tree UAO max contrast via int16 view object" : "valuedTree min-tree UAO max contrast via int16 view object");
        requireVectorEqual(collectImageValues(int16Uao.getAssociatedImage()), collectImageValues(valuedTreeUao.getAssociatedImage()),
                           isMaxtree ? "valuedTree max-tree UAO associated via int16 view object" : "valuedTree min-tree UAO associated via int16 view object");
        std::vector<float> floatAltitude;
        floatAltitude.reserve(externalAltitude.size());
        for (std::uint8_t level : externalAltitude) {
            floatAltitude.push_back(static_cast<float>(level) * 1.3f);
        }
        const ValuedMorphologicalTreeView<float> floatView(valuedTree->topology(), std::span<const float>(floatAltitude.data(), floatAltitude.size()));
        UltimateAttributeOpening<float> floatUao(floatView, attr);
        floatUao.execute(static_cast<int>(maximumAttributeThreshold));
        auto floatContrast = floatUao.getMaxContrastImage();
        static_assert(std::is_same_v<decltype(floatContrast), ImageFloatPtr>);
        const auto canonicalContrast = collectImageValues(valuedTreeUao.getMaxContrastImage());
        const auto floatContrastValues = collectImageValues(floatContrast);
        requireEqual(floatContrastValues.size(), canonicalContrast.size(), "valuedTree UAO float contrast size");
        bool hasNonZeroContrast = false;
        bool hasFractionalExpected = false;
        for (std::size_t i = 0; i < floatContrastValues.size(); ++i) {
            const float expected = static_cast<float>(canonicalContrast[i]) * 1.3f;
            requireNear(floatContrastValues[i], expected, 1.0e-5f, isMaxtree ? "valuedTree max-tree UAO float contrast" : "valuedTree min-tree UAO float contrast");
            hasNonZeroContrast = hasNonZeroContrast || canonicalContrast[i] != 0;
            hasFractionalExpected = hasFractionalExpected || std::abs(expected - std::round(expected)) > 1.0e-5f;
        }
        require(hasNonZeroContrast, "valuedTree UAO fixture must produce non-zero contrast");
        require(hasFractionalExpected, "valuedTree UAO float contrast test must exercise fractional output");
        requireThrows<std::logic_error>([&]() { viewUao.executeWithMSER(static_cast<int>(maximumAttributeThreshold), 1); },
                                        isMaxtree ? "valuedTree max-tree UAO view object must reject MSER"
                                                  : "valuedTree min-tree UAO view object must reject MSER");
        if (isMaxtree && contract::validationsEnabled) {
            requireThrows<std::invalid_argument>([&]() { valuedTreeUao.execute(static_cast<int>(maximumAttributeThreshold), std::vector<uint8_t>{true}); },
                                                 "UltimateAttributeOpening<std::uint8_t> must reject short selectedForFiltering buffer");
        }

        valuedTree->mergeNodeIntoParent(4);
        std::vector<bool> keepAllAfterMerge(valuedTree->topology().numInternalNodeSlots(), true);
        const NodePreservationMask keepAllAfterMergeMask(keepAllAfterMerge);
        auto expectedAfterMerge = ImageUInt8::create(valuedTree->topology().numRows(), valuedTree->topology().numColumns(), 0);
        for (PixelId pixel = 0; pixel < expectedAfterMerge->getSize(); ++pixel) {
            const NodeId nodeId = valuedTree->topology().smallestNode(pixel);
            (*expectedAfterMerge)[pixel] = static_cast<uint8_t>(valuedTree->nodeAltitude(nodeId));
        }

        auto directAfterMerge = applyDirectAttributeFilter(*valuedTree, keepAllAfterMergeMask);
        requireVectorEqual(collectImageValues(directAfterMerge), collectImageValues(expectedAfterMerge),
                           isMaxtree ? "valuedTree max-tree direct filter keep-all after merge" : "valuedTree min-tree direct filter keep-all after merge");

        auto subtractiveAfterMerge = applyHardSubtractiveAttributeFilter(*valuedTree, keepAllAfterMergeMask);
        requireVectorEqual(collectImageValues(subtractiveAfterMerge), collectImageValuesAs<AltitudeDifference<std::uint8_t>>(expectedAfterMerge),
                           isMaxtree ? "valuedTree max-tree hard subtractive filter keep-all after merge"
                                     : "valuedTree min-tree hard subtractive filter keep-all after merge");

        auto [mergedAttrNames, mergedAttr] = AttributeComputation::computeSingleAttribute(*valuedTree, BoundingBoxHeight);
        (void)mergedAttrNames;
        UltimateAttributeOpening<std::uint8_t> mergedValuedTreeUao(*valuedTree, mergedAttr);
        UltimateAttributeOpening<std::uint8_t> mergedValuedTreeUaoRaw(*valuedTree, mergedAttr.data());
        mergedValuedTreeUao.execute(static_cast<int>(maximumAttributeThreshold));
        mergedValuedTreeUaoRaw.execute(static_cast<int>(maximumAttributeThreshold));
        requireImageShape(mergedValuedTreeUao.getMaxContrastImage(), valuedTree->topology().numRows(),
                          valuedTree->topology().numColumns());
        requireImageShape(mergedValuedTreeUao.getAssociatedImage(), valuedTree->topology().numRows(),
                          valuedTree->topology().numColumns());
        requireImageShape(mergedValuedTreeUaoRaw.getMaxContrastImage(), valuedTree->topology().numRows(),
                          valuedTree->topology().numColumns());
        requireImageShape(mergedValuedTreeUaoRaw.getAssociatedImage(), valuedTree->topology().numRows(),
                          valuedTree->topology().numColumns());
    }

    auto nonSquare = makeImage(2, 3, {3, 3, 2, 1, 4, 5});
    auto nonSquareValuedTree = makeValuedComponentTree(nonSquare, true);
    auto [nonSquareAttrNames, nonSquareAttr] = AttributeComputation::computeSingleAttribute(*nonSquareValuedTree, BoundingBoxHeight);
    (void)nonSquareAttrNames;
    UltimateAttributeOpening<std::uint8_t> nonSquareUao(*nonSquareValuedTree, nonSquareAttr);
    nonSquareUao.execute(nonSquareValuedTree->topology().numRows());
    requireImageShape(nonSquareUao.getMaxContrastImage(), 2, 3);
    requireImageShape(nonSquareUao.getAssociatedImage(), 2, 3);
    requireImageShape(nonSquareUao.getAssociatedColorImage(), 2, 9);

    return 0;
}
