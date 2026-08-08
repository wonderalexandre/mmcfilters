#include "support/TestSupport.hpp"

#include "mmcfilters/filters/MSERComputer.hpp"
#include "mmcfilters/filters/AttributeFilters.hpp"
#include "mmcfilters/filters/UltimateAttributeOpening.hpp"
#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "mmcfilters/trees/WeightedTreeView.hpp"
#include "mmcfilters/trees/adjust/CasfComponentTrees.hpp"
#include "mmcfilters/trees/adjust/DualMinMaxTreeIncrementalFilter.hpp"
#include "mmcfilters/trees/adjust/DualMinMaxTreeIncrementalFilterLeaf.hpp"

#include <array>
#include <concepts>
#include <cstdint>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

namespace {

template <class T>
concept HasAdaptiveCriterion = requires(T& filters, std::vector<bool>& criterion) {
    { filters.getAdaptiveCriterion(criterion, 1) } -> std::same_as<std::vector<bool>>;
};

template <class T>
concept HasExecuteWithMSER = requires(T& uao) { uao.executeWithMSER(1, 1); };

template <class T>
concept HasPublicTrustedCommit = requires(T& editor) { editor.commitTrusted(); };

template <class Factory>
concept HasPublicTrustedNativeHierarchy = requires(NativeHierarchyView<std::uint8_t> hierarchy) { Factory::createFromTrustedNativeHierarchy(hierarchy); };

static_assert(std::is_same_v<decltype(MorphologicalTreeFactory::createMaxTree(std::declval<ImageUInt8Ptr>(), 1.5)), WeightedMorphologicalTree<std::uint8_t>>);
static_assert(std::is_same_v<decltype(MorphologicalTreeFactory::createMinTree(std::declval<ImageUInt8Ptr>(), 1.5)), WeightedMorphologicalTree<std::uint8_t>>);
static_assert(std::is_same_v<decltype(MorphologicalTreeFactory::createTreeOfShapes(std::declval<ImageUInt8Ptr>(), ToSInterpolation::SelfDual,
                                                                                   ToSDefaultInfinityRow, ToSDefaultInfinityCol)),
                             WeightedMorphologicalTree<std::uint8_t>>);
static_assert(std::is_constructible_v<MSERComputer<std::uint8_t>, const WeightedMorphologicalTree<std::uint8_t>&>);
static_assert(std::is_constructible_v<MSERComputer<std::int32_t>, const WeightedMorphologicalTree<std::int32_t>&>);
static_assert(!std::is_constructible_v<MSERComputer<std::uint8_t>, const MorphologicalTree&>);
static_assert(!std::is_constructible_v<MSERComputer<std::uint8_t>, const WeightedTreeView<std::uint8_t>&>);
static_assert(!std::is_constructible_v<MSERComputer<std::uint8_t>, const WeightedTreeView<std::int16_t>&>);

static_assert(HasAdaptiveCriterion<AttributeFilters<std::uint8_t>>);
static_assert(HasAdaptiveCriterion<AttributeFilters<std::int16_t>>);
static_assert(HasExecuteWithMSER<UltimateAttributeOpening<std::uint8_t>>);
static_assert(HasExecuteWithMSER<UltimateAttributeOpening<std::int16_t>>);
static_assert(!HasPublicTrustedCommit<WeightedTreeEditor<std::uint8_t>>);
static_assert(!HasPublicTrustedNativeHierarchy<MorphologicalTreeFactory>);
static_assert(std::is_move_constructible_v<detail::ValidatedNativeHierarchy<std::uint8_t>>);
static_assert(!std::is_copy_constructible_v<detail::ValidatedNativeHierarchy<std::uint8_t>>);
static_assert(!std::is_copy_constructible_v<detail::NativeTopologyProof>);

static_assert(std::is_same_v<decltype(std::declval<adjust::CasfComponentTrees<std::uint8_t>&>().minTree()), const WeightedMorphologicalTree<std::uint8_t>&>);
static_assert(std::is_same_v<decltype(std::declval<adjust::CasfComponentTrees<std::uint8_t>&>().maxTree()), const WeightedMorphologicalTree<std::uint8_t>&>);
static_assert(std::is_same_v<decltype(std::declval<adjust::CasfComponentTrees<std::int32_t>&>().minTree()), const WeightedMorphologicalTree<std::int32_t>&>);
static_assert(std::is_same_v<decltype(std::declval<adjust::CasfComponentTrees<float>&>().filter(std::declval<const std::vector<double>&>())), ImageFloatPtr>);

static_assert(std::is_constructible_v<adjust::DualMinMaxTreeIncrementalFilter<std::uint8_t>, WeightedMorphologicalTree<std::uint8_t>*,
                                      WeightedMorphologicalTree<std::uint8_t>*, RegularGridAdjacency2D&>);
static_assert(!std::is_constructible_v<adjust::DualMinMaxTreeIncrementalFilter<std::uint8_t>, WeightedTreeView<std::uint8_t>*, WeightedTreeView<std::uint8_t>*,
                                       RegularGridAdjacency2D&>);
static_assert(std::is_constructible_v<adjust::DualMinMaxTreeIncrementalFilter<std::int32_t>, WeightedMorphologicalTree<std::int32_t>*,
                                      WeightedMorphologicalTree<std::int32_t>*, RegularGridAdjacency2D&>);
static_assert(!std::is_constructible_v<adjust::DualMinMaxTreeIncrementalFilter<std::int32_t>, WeightedMorphologicalTree<std::uint8_t>*,
                                       WeightedMorphologicalTree<std::uint8_t>*, RegularGridAdjacency2D&>);
static_assert(std::is_constructible_v<adjust::DualMinMaxTreeIncrementalFilterLeaf<std::uint8_t>, WeightedMorphologicalTree<std::uint8_t>*,
                                      WeightedMorphologicalTree<std::uint8_t>*, RegularGridAdjacency2D&>);
static_assert(!std::is_constructible_v<adjust::DualMinMaxTreeIncrementalFilterLeaf<std::uint8_t>, WeightedTreeView<std::uint8_t>*,
                                       WeightedTreeView<std::uint8_t>*, RegularGridAdjacency2D&>);
static_assert(std::is_constructible_v<adjust::DualMinMaxTreeIncrementalFilterLeaf<std::int32_t>, WeightedMorphologicalTree<std::int32_t>*,
                                      WeightedMorphologicalTree<std::int32_t>*, RegularGridAdjacency2D&>);
static_assert(!std::is_constructible_v<adjust::DualMinMaxTreeIncrementalFilterLeaf<std::int32_t>, WeightedMorphologicalTree<std::uint8_t>*,
                                       WeightedMorphologicalTree<std::uint8_t>*, RegularGridAdjacency2D&>);

} // namespace

int main() {
    const auto image = makeComponentTreeFixture();
    auto minTree = MorphologicalTreeFactory::createMinTree(image);
    auto maxTree = MorphologicalTreeFactory::createMaxTree(image);
    auto adjacency = *minTree.topology().getUniformGridAdjacency2D();

    adjust::DualMinMaxTreeIncrementalFilter<std::uint8_t> adjust(&minTree, &maxTree, adjacency);
    adjust.pruneMinTreeAndUpdateMaxTree({});
    adjust.pruneMaxTreeAndUpdateMinTree({});
    require(minTree.topology().getRoot() != InvalidNode, "dual adjuster must keep a valid mutable min-tree owner");
    require(maxTree.topology().getRoot() != InvalidNode, "dual adjuster must keep a valid mutable max-tree owner");

    MSERComputer<std::uint8_t> mser(maxTree);
    requireEqual(static_cast<int>(mser.computeMSER(1).size()), maxTree.topology().getNumInternalNodeSlots(),
                 "MSER output must be indexed by the owner tree node slots");

    std::array<std::int32_t, 16> intPixels{
        3, 3, 2, 2, 3, 4, 4, 2, 1, 4, 5, 2, 1, 1, 5, 0,
    };
    auto intImage = ImageInt32::fromExternal(intPixels.data(), 4, 4);
    auto intMaxTree = MorphologicalTreeFactory::createMaxTree(intImage);
    MSERComputer<std::int32_t> intMser(intMaxTree);
    requireEqual(static_cast<int>(intMser.computeMSER(AltitudeDiff<std::int32_t>{1}).size()), intMaxTree.topology().getNumInternalNodeSlots(),
                 "typed int32 MSER output must be indexed by the owner tree node slots");
    auto intArea = AttributeComputation::computeSingleAttribute(intMaxTree, AREA);
    UltimateAttributeOpening<std::int32_t> intUao(intMaxTree, intArea.second);
    intUao.executeWithMSER(4, AltitudeDiff<std::int32_t>{1});
    requireEqual(intUao.getMaxContrastImage()->getSize(), intImage->getSize(), "typed int32 UAO MSER output must match image size");

    std::array<float, 16> floatPixels{
        0.3f, 0.3f, 0.2f, 0.2f, 0.3f, 0.4f, 0.4f, 0.2f, 0.1f, 0.4f, 0.5f, 0.2f, 0.1f, 0.1f, 0.5f, 0.0f,
    };
    auto floatImage = ImageFloat::fromExternal(floatPixels.data(), 4, 4);
    auto floatMaxTree = MorphologicalTreeFactory::createMaxTree(floatImage);
    MSERComputer<float> floatMser(floatMaxTree);
    requireEqual(static_cast<int>(floatMser.computeMSER(0.1f).size()), floatMaxTree.topology().getNumInternalNodeSlots(),
                 "typed float MSER output must be indexed by the owner tree node slots");
    auto floatArea = AttributeComputation::computeSingleAttribute(floatMaxTree, AREA);
    UltimateAttributeOpening<float> floatUao(floatMaxTree, floatArea.second);
    floatUao.executeWithMSER(4, 0.1f);
    requireEqual(floatUao.getMaxContrastImage()->getSize(), floatImage->getSize(), "typed float UAO MSER output must match image size");

    adjust::CasfComponentTrees<std::uint8_t> casf(image, adjust::CasfComponentTreesAttribute::AREA);
    require(&casf.minTree() != &casf.maxTree(), "CASF must own distinct min/max tree states");
    require(casf.minTree().topology().getDescriptiveKind() == MorphologicalTreeKind::MIN_TREE, "CASF min-tree accessor must expose the owner state");
    require(casf.maxTree().topology().getDescriptiveKind() == MorphologicalTreeKind::MAX_TREE, "CASF max-tree accessor must expose the owner state");

    return 0;
}
