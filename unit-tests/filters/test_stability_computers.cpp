#include "support/TestSupport.hpp"

#include "mmcfilters/filters/AttributeFilters.hpp"
#include "mmcfilters/filters/DepthStableRegionComputer.hpp"
#include "mmcfilters/filters/MSERComputer.hpp"
#include "mmcfilters/trees/detail/TreeStabilityNeighborhood.hpp"

#include <cstdint>
#include <cmath>
#include <iostream>
#include <span>
#include <vector>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

namespace {

WeightedMorphologicalTree<std::uint8_t> makeLinearMaxTreeFixture() {
    const std::vector<NodeId> parent{
        5, 6, 7, 8, 9,
        6, 7, 8, 9, 9};
    const std::vector<std::uint8_t> altitude{
        4, 3, 2, 1, 0,
        4, 3, 2, 1, 0};
    return MorphologicalTreeFactory::createFromHigraParent(
        std::span<const NodeId>(parent),
        std::span<const std::uint8_t>(altitude),
        1,
        5,
        MorphologicalTreeKind::MAX_TREE,
        AdjacencyRelation(1, 5, 1.0));
}

WeightedMorphologicalTree<std::uint8_t> makeImportedTreeOfShapesFixture() {
    const std::vector<NodeId> parent{
        4, 5, 5, 6,
        6, 6, 7, 7};
    const std::vector<std::uint8_t> altitude{
        11, 42, 42, 7,
        11, 42, 7, 99};
    return MorphologicalTreeFactory::createFromHigraParent(
        std::span<const NodeId>(parent),
        std::span<const std::uint8_t>(altitude),
        1,
        4,
        MorphologicalTreeKind::TREE_OF_SHAPES);
}

void verifyClassicalMserContracts() {
    auto linear = makeLinearMaxTreeFixture();
    const MorphologicalTree& tree = linear.topology();

    MSERComputer<std::uint8_t> mser(linear);
    requireThrows<std::invalid_argument>(
        [&]() { static_cast<void>(mser.computeMSER(0)); },
        "MSERComputer must reject zero altitude delta");

    const std::vector<uint8_t> flags = mser.computeMSER(1);
    const std::vector<float>& variation = mser.getVariations();
    requireEqual(static_cast<int>(flags.size()), tree.getNumInternalNodeSlots(), "MSER flags size");
    requireNear(variation[1], 1.0f, 1e-6f, "linear MSER node 1 variation");
    requireNear(variation[2], 2.0f / 3.0f, 1e-6f, "linear MSER node 2 variation");
    requireNear(variation[3], 0.5f, 1e-6f, "linear MSER node 3 variation");

    std::vector<bool> rejectedAtMiddle(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), true);
    rejectedAtMiddle[2] = false;
    const std::vector<bool> adaptive =
        AttributeFilters<std::uint8_t>::getAdaptiveCriterion(linear, rejectedAtMiddle, AltitudeDiff<std::uint8_t>{1});
    require(adaptive[3], "MSER adaptive criterion must choose the minimum-variation ascendant");
    require(!adaptive[1], "MSER adaptive criterion must not choose the maximum-variation descendant");
    require(!adaptive[2], "MSER adaptive criterion must move the rejected node when a more stable neighbour exists");
}

void verifyDepthStabilityContracts() {
    auto tos = makeImportedTreeOfShapesFixture();
    const MorphologicalTree& tree = tos.topology();

    requireThrows<std::invalid_argument>(
        [&]() { MSERComputer<std::uint8_t> invalid(tos); },
        "MSERComputer must reject tree of shapes");
    requireThrows<std::invalid_argument>(
        [&]() { static_cast<void>(detail::computeDepthStabilityNeighborhood(tree, 0)); },
        "depth stability must reject zero depth delta");

    const detail::StabilityNeighborhood depth1 =
        detail::computeDepthStabilityNeighborhood(tree, 1);
    requireVectorEqual(depth1.ascendants, {2, 2, 3, InvalidNode}, "ToS depth-1 ascendants");
    requireVectorEqual(depth1.descendants, {InvalidNode, InvalidNode, 1, 2}, "ToS depth-1 descendants");

    const detail::StabilityNeighborhood depth2 =
        detail::computeDepthStabilityNeighborhood(tree, 2);
    requireVectorEqual(depth2.ascendants, {3, 3, InvalidNode, InvalidNode}, "ToS depth-2 ascendants");
    requireVectorEqual(depth2.descendants, {InvalidNode, InvalidNode, InvalidNode, 1}, "ToS depth-2 descendants");

    DepthStableRegionComputer<float> depthStability(tree);
    (void)depthStability.computeByDepth(1);
    requireNear(depthStability.getVariations()[2], 0.5f, 1e-6f, "ToS depth stability uses topology-only area");

    std::vector<bool> rejectedAtMiddle(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), true);
    rejectedAtMiddle[2] = false;
    const std::vector<bool> adaptive =
        AttributeFilters<std::uint8_t>::getAdaptiveCriterionByDepth(tos, rejectedAtMiddle, 1);
    require(adaptive[2], "depth adaptive criterion prunes the rejected node when neighbours have no finite variation");
}

void verifyResidualTreeContract() {
    auto image = makeImage(1, 2, {0, 1});
    auto residual = MorphologicalTreeFactory::createSelfDualResidualTree(image, 1.0);
    require(
        residual.topology().getTreeType() == MorphologicalTreeKind::SELF_DUAL_RESIDUAL_TREE,
        "residual fixture tree type");
    requireThrows<std::invalid_argument>(
        [&]() { MSERComputer<std::uint8_t> invalid(residual); },
        "MSERComputer must reject self-dual residual trees");

    DepthStableRegionComputer<float> stability(residual.topology());
    requireThrows<std::invalid_argument>(
        [&]() { static_cast<void>(stability.computeByDepth(0)); },
        "residual depth stability must reject zero depth delta");
    static_cast<void>(stability.computeByDepth(1));
}

} // namespace

int main() {
    try {
        verifyClassicalMserContracts();
        verifyDepthStabilityContracts();
        verifyResidualTreeContract();
    } catch (const std::exception& ex) {
        std::cerr << "Stability computer test failed: " << ex.what() << "\n";
        return 1;
    }
    return 0;
}
