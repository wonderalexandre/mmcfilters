#include "support/TestSupport.hpp"

#include "mmcfilters/filters/DepthStableRegionComputer.hpp"
#include "mmcfilters/filters/MSERComputer.hpp"
#include "mmcfilters/filters/NodePreservationStability.hpp"
#include "mmcfilters/trees/detail/TreeStabilityNeighborhood.hpp"

#include <cstdint>
#include <cmath>
#include <iostream>
#include <limits>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

namespace {

ValuedMorphologicalTree<std::uint8_t> makeLinearMaxTreeFixture() {
    const std::vector<NodeId> parent{5, 6, 7, 8, 9, 6, 7, 8, 9, 9};
    const std::vector<std::uint8_t> altitude{4, 3, 2, 1, 0, 4, 3, 2, 1, 0};
    return MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(parent), std::span<const std::uint8_t>(altitude), 1, 5,
                                                           MorphologicalTreeKind::MaxTree, RegularGridAdjacency2D(1, 5, 1.0));
}

ValuedMorphologicalTree<std::uint8_t> makeImportedTreeOfShapesFixture() {
    const std::vector<NodeId> parent{4, 5, 5, 6, 6, 6, 7, 7};
    const std::vector<std::uint8_t> altitude{11, 42, 42, 7, 11, 42, 7, 99};
    return MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(parent), std::span<const std::uint8_t>(altitude), 1, 4,
                                                           MorphologicalTreeKind::TreeOfShapes);
}

ValuedMorphologicalTree<std::uint8_t> makeSpatialTieFixture(bool smallestSpatialSupportUsesLowerNodeId) {
    const std::vector<NodeId> parent{0, 0, 0};
    const std::vector<NodeId> smallestNodeMap = smallestSpatialSupportUsesLowerNodeId ? std::vector<NodeId>{1, 1, 2, 2}
                                                                                     : std::vector<NodeId>{2, 2, 1, 1};
    const std::vector<std::uint8_t> altitude{0, 1, 1};
    return MorphologicalTreeFactory::createFromNativeTopology(
        std::span<const NodeId>(parent), std::span<const NodeId>(smallestNodeMap), std::span<const std::uint8_t>(altitude), 0, 1, 4,
        MorphologicalTreeSemantics{MorphologicalTreeKind::Generic, NodeAltitudeOrder::Increasing, NoConstructionContext{}});
}

void verifyStableTiePolicies() {
    for (const bool smallestSpatialSupportUsesLowerNodeId : {true, false}) {
        auto valuedTree = makeSpatialTieFixture(smallestSpatialSupportUsesLowerNodeId);
        const MorphologicalTree& tree = valuedTree.topology();
        const NodeId expectedRepresentative = tree.smallestNode(PixelId{0});

        const detail::StabilityNeighborhood depthNeighborhood = detail::computeDepthStabilityNeighborhood(tree, 1);
        requireEqual(depthNeighborhood.descendants[static_cast<std::size_t>(tree.root())], expectedRepresentative,
                     "depth-stability equal-support tie must use the smallest row-major support pixel");

        const detail::StabilityNeighborhood altitudeNeighborhood =
            detail::computeAltitudeStabilityNeighborhood(tree, std::span<const std::uint8_t>(valuedTree.nodeAltitudes()),
                                                          AltitudeDifference<std::uint8_t>{1});
        requireEqual(altitudeNeighborhood.descendants[static_cast<std::size_t>(tree.root())], expectedRepresentative,
                     "altitude-stability equal-support tie must use the smallest row-major support pixel");

        const NodeId otherRepresentative = expectedRepresentative == NodeId{1} ? NodeId{2} : NodeId{1};
        auto editor = valuedTree.edit();
        editor.movePixelToProperPart(otherRepresentative, expectedRepresentative, PixelId{1});
        const TreeValidationResult commitResult = editor.validateAndCommit();
        require(commitResult.ok, "support-metadata cache invalidation fixture must remain a valid committed tree");

        const detail::StabilityNeighborhood afterProperPartEdit = detail::computeDepthStabilityNeighborhood(tree, 1);
        requireEqual(afterProperPartEdit.descendants[static_cast<std::size_t>(tree.root())], otherRepresentative,
                     "stability descendant selection must refresh cached support metadata after a proper-part edit");
    }

    std::vector<float> variation(6, std::numeric_limits<float>::infinity());
    std::vector<NodeId> ancestors(6, InvalidNode);
    std::vector<NodeId> descendants(6, InvalidNode);
    constexpr NodeId center = 5;
    constexpr NodeId descendant = 1;
    constexpr NodeId ancestor = 2;
    descendants[center] = descendant;
    ancestors[center] = ancestor;
    variation[center] = 0.25f;
    variation[descendant] = 0.25f;
    variation[ancestor] = 0.25f;

    requireEqual(detail::nodeWithMinimumVariationInWindow(center, variation, ancestors, descendants), center,
                 "equal finite variations must retain the current node");

    variation[center] = std::numeric_limits<float>::quiet_NaN();
    requireEqual(detail::nodeWithMinimumVariationInWindow(center, variation, ancestors, descendants), descendant,
                 "an unavailable center must prefer the representative descendant on an equal-variation tie");

    variation[descendant] = std::numeric_limits<float>::quiet_NaN();
    requireEqual(detail::nodeWithMinimumVariationInWindow(center, variation, ancestors, descendants), ancestor,
                 "an unavailable center and descendant must select the finite ancestor");
}

void verifyClassicalMserContracts() {
    auto linear = makeLinearMaxTreeFixture();
    const MorphologicalTree& tree = linear.topology();

    MSERComputer<std::uint8_t> mser(linear);
    if constexpr (contract::validationsEnabled) {
        requireThrows<std::logic_error>([&]() { static_cast<void>(mser.getVariation(0)); }, "MSER variation requires computation");
        requireThrows<std::logic_error>([&]() { static_cast<void>(mser.getVariations()); }, "MSER variation buffer requires computation");
        requireThrows<std::logic_error>([&]() { static_cast<void>(mser.numNodes()); }, "MSER selected count requires computation");
        requireThrows<std::logic_error>([&]() { static_cast<void>(mser.nodeWithMinimumVariationInWindow(0)); }, "MSER minimum requires computation");
        requireThrows<std::logic_error>([&]() { static_cast<void>(mser.ancestorInStabilityWindow(0)); }, "MSER ancestor requires computation");
        requireThrows<std::logic_error>([&]() { static_cast<void>(mser.descendantInStabilityWindow(0)); }, "MSER descendant requires computation");
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(mser.computeMSER(0)); }, "MSERComputer must reject a zero altitude-window radius");
    }

    const std::vector<uint8_t> flags = mser.computeMSER(1);
    const std::vector<float>& variation = mser.getVariations();
    requireEqual(static_cast<int>(flags.size()), tree.numInternalNodeSlots(), "MSER flags size");
    requireNear(variation[1], 1.0f, 1e-6f, "linear MSER node 1 variation");
    requireNear(variation[2], 2.0f / 3.0f, 1e-6f, "linear MSER node 2 variation");
    requireNear(variation[3], 0.5f, 1e-6f, "linear MSER node 3 variation");

    const NodePreservationMask allPreserved(static_cast<std::size_t>(tree.numInternalNodeSlots()), true);
    const NodePreservationMask allPreservedAdjusted =
        adjustNodePreservationMaskByAltitudeStability(linear, allPreserved, AltitudeDifference<std::uint8_t>{1});
    requireVectorEqual(allPreservedAdjusted.decisions(), std::vector<bool>(allPreserved.size(), true),
                       "altitude stability must retain an all-preserved mask");

    std::vector<bool> rejectedAtMiddle(static_cast<std::size_t>(tree.numInternalNodeSlots()), true);
    rejectedAtMiddle[2] = false;
    const NodePreservationMask adjusted = adjustNodePreservationMaskByAltitudeStability(
        linear, NodePreservationMask(rejectedAtMiddle), AltitudeDifference<std::uint8_t>{1});
    require(!adjusted[3], "altitude stability must relocate rejection to the minimum-variation ancestor");
    require(adjusted[1], "altitude stability must preserve the higher-variation descendant");
    require(adjusted[2], "altitude stability must restore the input node after relocating its rejection");

    std::vector<bool> duplicateRelocationInput(static_cast<std::size_t>(tree.numInternalNodeSlots()), true);
    duplicateRelocationInput[2] = false;
    duplicateRelocationInput[3] = false;
    const NodePreservationMask duplicateRelocation = adjustNodePreservationMaskByAltitudeStability(
        linear, NodePreservationMask(duplicateRelocationInput), AltitudeDifference<std::uint8_t>{1});
    require(duplicateRelocation[2] && !duplicateRelocation[3],
            "duplicate altitude-stability relocations must collapse into one rejected target");

    const NodePreservationMask allRejected(static_cast<std::size_t>(tree.numInternalNodeSlots()), false);
    const NodePreservationMask allRejectedAdjusted = adjustNodePreservationMaskByAltitudeStability(
        linear, allRejected, AltitudeDifference<std::uint8_t>{1}, IncompleteStabilityWindowPolicy::PreserveInputDecision);
    requireVectorEqual(allRejectedAdjusted.decisions(), std::vector<bool>{false, true, false, false, false},
                       "all-rejected altitude input must use preservation polarity after relocation");

    const std::vector<float> nodeAttributes{1.0f, 2.0f, 2.0f, 3.0f, 4.0f};
    const NodePreservationMask thresholdMask = computeNodePreservationMask(std::span<const float>(nodeAttributes), 2.0f);
    requireVectorEqual(thresholdMask.decisions(), std::vector<bool>{false, true, true, true, true},
                       "node-preservation thresholding must use the inclusive greater-than-or-equal rule");

    MSERComputer<std::uint8_t> ownedCopy = [&]() {
        MSERComputer<std::uint8_t> source(linear, std::vector<float>{10.0f, 20.0f, 30.0f, 40.0f, 50.0f});
        return MSERComputer<std::uint8_t>(source);
    }();
    requireNear(ownedCopy.getAttrMSER(2), 30.0f, 1e-6f, "copied MSER owns an independent attribute buffer");

    std::vector<float> borrowedAttributes{1.0f, 2.0f, 3.0f, 4.0f, 5.0f};
    MSERComputer<std::uint8_t> borrowed(linear, borrowedAttributes.data());
    MSERComputer<std::uint8_t> borrowedCopy(borrowed);
    borrowedAttributes[2] = 7.0f;
    requireNear(borrowedCopy.getAttrMSER(2), 7.0f, 1e-6f, "copied MSER preserves an explicit borrowed view");
}

void verifyDepthStabilityContracts() {
    auto tos = makeImportedTreeOfShapesFixture();
    const MorphologicalTree& tree = tos.topology();

    requireThrows<std::invalid_argument>([&]() { MSERComputer<std::uint8_t> invalid(tos); }, "MSERComputer must reject unconstrained altitude order");
    if constexpr (contract::validationsEnabled) {
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(detail::computeDepthStabilityNeighborhood(tree, 0)); },
                                             "depth stability must reject a zero depth-window radius");
    }

    const detail::StabilityNeighborhood depth1 = detail::computeDepthStabilityNeighborhood(tree, 1);
    requireVectorEqual(depth1.ancestors, {2, 2, 3, InvalidNode}, "ToS depth-1 ancestors");
    requireVectorEqual(depth1.descendants, {InvalidNode, InvalidNode, 1, 2}, "ToS depth-1 descendants");

    const detail::StabilityNeighborhood depth2 = detail::computeDepthStabilityNeighborhood(tree, 2);
    requireVectorEqual(depth2.ancestors, {3, 3, InvalidNode, InvalidNode}, "ToS depth-2 ancestors");
    requireVectorEqual(depth2.descendants, {InvalidNode, InvalidNode, InvalidNode, 1}, "ToS depth-2 descendants");

    DepthStableRegionComputer<float> depthStability(tree);
    if constexpr (contract::validationsEnabled) {
        requireThrows<std::logic_error>([&]() { static_cast<void>(depthStability.getVariation(0)); }, "depth variation requires computation");
        requireThrows<std::logic_error>([&]() { static_cast<void>(depthStability.getVariations()); }, "depth variation buffer requires computation");
        requireThrows<std::logic_error>([&]() { static_cast<void>(depthStability.numNodes()); }, "depth selected count requires computation");
        requireThrows<std::logic_error>([&]() { static_cast<void>(depthStability.nodeWithMinimumVariationInWindow(0)); }, "depth minimum requires computation");
        requireThrows<std::logic_error>([&]() { static_cast<void>(depthStability.ancestorInStabilityWindow(0)); }, "depth ancestor requires computation");
        requireThrows<std::logic_error>([&]() { static_cast<void>(depthStability.descendantInStabilityWindow(0)); }, "depth descendant requires computation");
    }
    (void)depthStability.computeByDepth(1);
    requireNear(depthStability.getVariations()[2], 0.5f, 1e-6f, "ToS depth stability uses topology-only area");

    const NodePreservationMask allPreserved(static_cast<std::size_t>(tree.numInternalNodeSlots()), true);
    const NodePreservationMask allPreservedAdjusted = adjustNodePreservationMaskByDepthStability(tos, allPreserved, 1);
    requireVectorEqual(allPreservedAdjusted.decisions(), std::vector<bool>(allPreserved.size(), true),
                       "depth stability must retain an all-preserved mask");

    std::vector<bool> rejectedAtMiddle(static_cast<std::size_t>(tree.numInternalNodeSlots()), true);
    rejectedAtMiddle[2] = false;
    const NodePreservationMask adjusted = adjustNodePreservationMaskByDepthStability(
        tos, NodePreservationMask(rejectedAtMiddle), 1, IncompleteStabilityWindowPolicy::PreserveInputDecision);
    require(!adjusted[2], "an incomplete depth-stability window must retain the input rejection");
    require(adjusted[0] && adjusted[1] && adjusted[3], "depth stability must return preservation polarity");

    if constexpr (contract::validationsEnabled) {
        requireThrows<std::invalid_argument>(
            [&]() { static_cast<void>(adjustNodePreservationMaskByAltitudeStability(tos, allPreserved, AltitudeDifference<std::uint8_t>{1})); },
            "altitude stability must reject an unconstrained-altitude tree");
    }

    DepthStableRegionComputer<float> ownedCopy = [&]() {
        DepthStableRegionComputer<float> source(tree, std::vector<float>{11.0f, 22.0f, 33.0f, 44.0f});
        return DepthStableRegionComputer<float>(source);
    }();
    requireNear(ownedCopy.getAttribute(2), 33.0f, 1e-6f, "copied depth evaluator owns an independent attribute buffer");

    std::vector<float> borrowedAttributes{1.0f, 2.0f, 3.0f, 4.0f};
    DepthStableRegionComputer<float> borrowed(tree, borrowedAttributes.data());
    DepthStableRegionComputer<float> borrowedCopy(borrowed);
    borrowedAttributes[2] = 8.0f;
    requireNear(borrowedCopy.getAttribute(2), 8.0f, 1e-6f, "copied depth evaluator preserves an explicit borrowed view");
}

} // namespace

int main() {
    static_assert(std::is_copy_constructible_v<MSERComputer<std::uint8_t>>);
    static_assert(std::is_move_constructible_v<MSERComputer<std::uint8_t>>);
    static_assert(!std::is_copy_assignable_v<MSERComputer<std::uint8_t>>);
    static_assert(!std::is_move_assignable_v<MSERComputer<std::uint8_t>>);
    static_assert(std::is_copy_constructible_v<DepthStableRegionComputer<float>>);
    static_assert(std::is_move_constructible_v<DepthStableRegionComputer<float>>);
    static_assert(!std::is_copy_assignable_v<DepthStableRegionComputer<float>>);
    static_assert(!std::is_move_assignable_v<DepthStableRegionComputer<float>>);

    try {
        verifyStableTiePolicies();
        verifyClassicalMserContracts();
        verifyDepthStabilityContracts();
    } catch (const std::exception& ex) {
        std::cerr << "Stability computer test failed: " << ex.what() << "\n";
        return 1;
    }
    return 0;
}
