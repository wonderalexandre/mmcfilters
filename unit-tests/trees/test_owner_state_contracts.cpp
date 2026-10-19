#include "support/TestSupport.hpp"

#include "mmcfilters/filters/MSERComputer.hpp"
#include "mmcfilters/filters/AttributeFilters.hpp"
#include "mmcfilters/filters/NodePreservationStability.hpp"
#include "mmcfilters/filters/UltimateAttributeOpening.hpp"
#include "mmcfilters/trees/NativeHierarchyOwned.hpp"
#include "mmcfilters/trees/NativeHierarchyView.hpp"
#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "mmcfilters/trees/ValuedMorphologicalTreeView.hpp"
#include "mmcfilters/trees/adjust/CasfComponentTrees.hpp"
#include "mmcfilters/trees/adjust/DualMinMaxTreeIncrementalFilter.hpp"
#include "mmcfilters/trees/detail/DualMinMaxTreeIncrementalFilterLeaf.hpp"

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
concept HasExecuteWithMSER = requires(T& uao) { uao.executeWithMSER(1, 1); };

template <class T>
concept HasPublicTrustedCommit = requires(T& editor) { editor.commitTrusted(); };

template <class Factory>
concept HasPublicTrustedNativeHierarchy = requires(NativeHierarchyView<std::uint8_t> hierarchy) { Factory::createFromTrustedNativeHierarchy(hierarchy); };

static_assert(std::is_same_v<decltype(MorphologicalTreeFactory::createMaxTree(std::declval<ImageUInt8Ptr>(), 1.5)), ValuedMorphologicalTree<std::uint8_t>>);
static_assert(std::is_same_v<decltype(MorphologicalTreeFactory::createMinTree(std::declval<ImageUInt8Ptr>(), 1.5)), ValuedMorphologicalTree<std::uint8_t>>);
static_assert(std::is_same_v<decltype(MorphologicalTreeFactory::createTreeOfShapes(std::declval<ImageUInt8Ptr>(), TopographicConvention{})),
                             ValuedMorphologicalTree<std::uint8_t>>);
static_assert(std::is_same_v<decltype(MorphologicalTreeFactory::createTreeOfShapes<ToSGrayLevel>(std::declval<ImageUInt8Ptr>(), TopographicConvention{})),
                             ValuedMorphologicalTree<ToSGrayLevel>>);
static_assert(std::is_constructible_v<MSERComputer<std::uint8_t>, const ValuedMorphologicalTree<std::uint8_t>&>);
static_assert(std::is_constructible_v<MSERComputer<std::int32_t>, const ValuedMorphologicalTree<std::int32_t>&>);
static_assert(!std::is_constructible_v<MSERComputer<std::uint8_t>, const MorphologicalTree&>);
static_assert(!std::is_constructible_v<MSERComputer<std::uint8_t>, const ValuedMorphologicalTreeView<std::uint8_t>&>);
static_assert(!std::is_constructible_v<MSERComputer<std::uint8_t>, const ValuedMorphologicalTreeView<std::int16_t>&>);

static_assert(std::is_same_v<decltype(computeNodePreservationMask(std::declval<std::span<const float>>(), 1.0f)),
                             NodePreservationMask>);
static_assert(std::is_same_v<decltype(adjustNodePreservationMaskByAltitudeStability(
                                 std::declval<const ValuedMorphologicalTree<std::uint8_t>&>(),
                                 std::declval<const NodePreservationMask&>(), AltitudeDifference<std::uint8_t>{1})),
                             NodePreservationMask>);
static_assert(std::is_same_v<decltype(adjustNodePreservationMaskByDepthStability(
                                 std::declval<const MorphologicalTree&>(), std::declval<const NodePreservationMask&>(), 1)),
                             NodePreservationMask>);
static_assert(HasExecuteWithMSER<UltimateAttributeOpening<std::uint8_t>>);
static_assert(HasExecuteWithMSER<UltimateAttributeOpening<std::int16_t>>);
static_assert(!HasPublicTrustedCommit<ValuedMorphologicalTreeEditor<std::uint8_t>>);
template <class Factory>
concept HasPublicNativeHierarchyImport = requires(NativeHierarchyView<std::uint8_t> hierarchy) { Factory::createFromNativeHierarchy(hierarchy); };
template <class Factory>
concept HasPublicNativeTopologyImport = requires(std::span<const NodeId> parent, std::span<const std::uint8_t> altitude) {
    Factory::createFromNativeTopology(parent, parent, altitude, NodeId{0}, MorphologicalTreeSemantics{});
};
template <class Factory>
concept HasPublicGridNativeTopologyImport = requires(std::span<const NodeId> parent, std::span<const std::uint8_t> altitude) {
    Factory::createFromNativeTopology(parent, parent, altitude, NodeId{0}, 1, 1, MorphologicalTreeSemantics{});
};
static_assert(!HasPublicNativeHierarchyImport<MorphologicalTreeFactory>);
static_assert(!HasPublicNativeTopologyImport<MorphologicalTreeFactory>);
static_assert(!HasPublicGridNativeTopologyImport<MorphologicalTreeFactory>);
static_assert(!HasPublicTrustedNativeHierarchy<MorphologicalTreeFactory>);
static_assert(!std::is_constructible_v<ValuedMorphologicalTreeView<std::uint8_t>, const ValuedMorphologicalTree<std::uint8_t>&>);
static_assert(!std::is_default_constructible_v<detail::MorphologicalTreeConstructionTag>);
static_assert(std::is_nothrow_move_constructible_v<NativeHierarchyOwned<std::uint8_t>>);
static_assert(std::is_nothrow_move_assignable_v<NativeHierarchyOwned<std::uint8_t>>);
static_assert(!std::is_copy_constructible_v<NativeHierarchyOwned<std::uint8_t>>);
static_assert(!std::is_copy_assignable_v<NativeHierarchyOwned<std::uint8_t>>);

static_assert(std::is_same_v<decltype(std::declval<adjust::CasfComponentTrees<std::uint8_t>&>().minTree()), const ValuedMorphologicalTree<std::uint8_t>&>);
static_assert(std::is_same_v<decltype(std::declval<adjust::CasfComponentTrees<std::uint8_t>&>().maxTree()), const ValuedMorphologicalTree<std::uint8_t>&>);
static_assert(std::is_same_v<decltype(std::declval<adjust::CasfComponentTrees<std::int32_t>&>().minTree()), const ValuedMorphologicalTree<std::int32_t>&>);
static_assert(std::is_same_v<decltype(std::declval<adjust::CasfComponentTrees<float>&>().filter(std::declval<const std::vector<double>&>())), ImageFloatPtr>);

static_assert(std::is_constructible_v<adjust::DualMinMaxTreeIncrementalFilter<std::uint8_t>, ValuedMorphologicalTree<std::uint8_t>*,
                                      ValuedMorphologicalTree<std::uint8_t>*, RegularGridAdjacency2D&>);
static_assert(!std::is_constructible_v<adjust::DualMinMaxTreeIncrementalFilter<std::uint8_t>, ValuedMorphologicalTreeView<std::uint8_t>*, ValuedMorphologicalTreeView<std::uint8_t>*,
                                       RegularGridAdjacency2D&>);
static_assert(std::is_constructible_v<adjust::DualMinMaxTreeIncrementalFilter<std::int32_t>, ValuedMorphologicalTree<std::int32_t>*,
                                      ValuedMorphologicalTree<std::int32_t>*, RegularGridAdjacency2D&>);
static_assert(!std::is_constructible_v<adjust::DualMinMaxTreeIncrementalFilter<std::int32_t>, ValuedMorphologicalTree<std::uint8_t>*,
                                       ValuedMorphologicalTree<std::uint8_t>*, RegularGridAdjacency2D&>);
static_assert(std::is_constructible_v<detail::DualMinMaxTreeIncrementalFilterLeaf<std::uint8_t>, ValuedMorphologicalTree<std::uint8_t>*,
                                      ValuedMorphologicalTree<std::uint8_t>*, RegularGridAdjacency2D&>);
static_assert(!std::is_constructible_v<detail::DualMinMaxTreeIncrementalFilterLeaf<std::uint8_t>, ValuedMorphologicalTreeView<std::uint8_t>*,
                                       ValuedMorphologicalTreeView<std::uint8_t>*, RegularGridAdjacency2D&>);
static_assert(std::is_constructible_v<detail::DualMinMaxTreeIncrementalFilterLeaf<std::int32_t>, ValuedMorphologicalTree<std::int32_t>*,
                                      ValuedMorphologicalTree<std::int32_t>*, RegularGridAdjacency2D&>);
static_assert(!std::is_constructible_v<detail::DualMinMaxTreeIncrementalFilterLeaf<std::int32_t>, ValuedMorphologicalTree<std::uint8_t>*,
                                       ValuedMorphologicalTree<std::uint8_t>*, RegularGridAdjacency2D&>);

void requireNativeHierarchyBufferTransfer() {
    const NodeId* parentData = nullptr;
    const NodeId* smallestNodeData = nullptr;
    const std::uint8_t* altitudeData = nullptr;
    auto hierarchy = [&] {
        std::vector<NodeId> parent{0, 0};
        std::vector<NodeId> smallestNodeMap{0, 1, 1};
        std::vector<std::uint8_t> altitudes{2, 5};
        parentData = parent.data();
        smallestNodeData = smallestNodeMap.data();
        altitudeData = altitudes.data();
        return makeNativeHierarchyOwned(
            std::move(parent), std::move(smallestNodeMap), std::move(altitudes), NodeId{0}, GridDomain2D{1, 3},
            makeMorphologicalTreeSemantics(MorphologicalTreeKind::MaxTree,
                                          SharedAdjacencyContext{RegularGridAdjacency2D::horizontalLine(1, 3, 2)}));
    }();
    const auto requireOriginalBuffers = [&](const NativeHierarchyOwned<std::uint8_t>& owner) {
        require(owner.parent.data() == parentData, "native hierarchy transfers the parent allocation without copying");
        require(owner.smallestNodeMap.data() == smallestNodeData, "native hierarchy transfers the smallest-node allocation without copying");
        require(owner.nodeAltitudes.data() == altitudeData, "native hierarchy transfers the altitude allocation without copying");
    };
    requireOriginalBuffers(hierarchy);
    auto moved = std::move(hierarchy);
    requireOriginalBuffers(moved);
    auto assigned = makeNativeHierarchyOwned(std::vector<NodeId>{0}, std::vector<NodeId>{0}, std::vector<std::uint8_t>{1}, NodeId{0},
                                                    std::nullopt, makeMorphologicalTreeSemantics(MorphologicalTreeKind::Generic));
    assigned = std::move(moved);
    requireOriginalBuffers(assigned);
    requireEqual(assigned.root, NodeId{0}, "native hierarchy move assignment retains the root");
    require(assigned.gridDomain2D.has_value(), "native hierarchy move assignment retains the grid domain");
    requireEqual(assigned.gridDomain2D->rows, 1, "native hierarchy move assignment retains grid rows");
    requireEqual(assigned.gridDomain2D->columns, 3, "native hierarchy move assignment retains grid columns");
    require(assigned.semantics.kind == MorphologicalTreeKind::MaxTree, "native hierarchy move assignment retains its kind");
    const auto* context = std::get_if<SharedAdjacencyContext>(&assigned.semantics.constructionContext);
    require(context != nullptr, "native hierarchy move assignment retains its owned adjacency");
    requireVectorEqual(collectPixelIds(context->adjacency.getNeighborIndices(0)), std::vector<PixelId>{1, 2},
                       "native hierarchy move assignment retains the construction stencil after source destruction");

    auto imported = [&] {
        auto source = std::move(assigned);
        std::vector<NodeId> higraParent{3, 4, 4, 3, 3};
        std::vector<std::uint8_t> higraAltitude{2, 5, 5, 2, 5};
        auto tree = MorphologicalTreeFactory::createFromHigraParent(
            std::span<const NodeId>(higraParent), std::span<const std::uint8_t>(higraAltitude), 1, 3,
            MorphologicalTreeKind::MaxTree, RegularGridAdjacency2D(1, 3, 2.0));
        higraParent[0] = InvalidNode;
        higraAltitude[0] = 255;
        return tree;
    }();
    requireEqual(imported.topology().parent(1), NodeId{0}, "Higra import owns its topology after source mutation and destruction");
    requireVectorEqual(collectImageValues(imported.reconstructFromNodeAltitudes()), std::vector<std::uint8_t>{2, 5, 5},
                       "Higra import owns its pixel assignments and altitudes after source mutation and destruction");
    requireVectorEqual(collectPixelIds(imported.topology().sharedAdjacencyContext()->adjacency.getNeighborIndices(0)), std::vector<PixelId>{1, 2},
                       "Higra import owns its construction adjacency after source destruction");
}

} // namespace

int main() {
    requireNativeHierarchyBufferTransfer();

    const auto image = makeComponentTreeFixture();
    auto minTree = MorphologicalTreeFactory::createMinTree(image);
    auto maxTree = MorphologicalTreeFactory::createMaxTree(image);
    auto adjacency = minTree.topology().sharedAdjacencyContext()->adjacency;

    adjust::DualMinMaxTreeIncrementalFilter<std::uint8_t> adjust(&minTree, &maxTree, adjacency);
    adjust.pruneMinTreeAndUpdateMaxTree({});
    adjust.pruneMaxTreeAndUpdateMinTree({});
    require(minTree.topology().root() != InvalidNode, "dual adjuster must keep a valid mutable min-tree owner");
    require(maxTree.topology().root() != InvalidNode, "dual adjuster must keep a valid mutable max-tree owner");

    MSERComputer<std::uint8_t> mser(maxTree);
    requireEqual(static_cast<int>(mser.computeMSER(1).size()), maxTree.topology().numInternalNodeSlots(),
                 "MSER output must be indexed by the owner tree node slots");

    std::array<std::int32_t, 16> intPixels{
        3, 3, 2, 2, 3, 4, 4, 2, 1, 4, 5, 2, 1, 1, 5, 0,
    };
    auto intImage = ImageInt32::fromExternal(intPixels.data(), 4, 4);
    auto intMaxTree = MorphologicalTreeFactory::createMaxTree(intImage);
    MSERComputer<std::int32_t> intMser(intMaxTree);
    requireEqual(static_cast<int>(intMser.computeMSER(AltitudeDifference<std::int32_t>{1}).size()), intMaxTree.topology().numInternalNodeSlots(),
                 "typed int32 MSER output must be indexed by the owner tree node slots");
    auto intArea = AttributeComputation::computeSingleAttribute(intMaxTree, Area);
    UltimateAttributeOpening<std::int32_t> intUao(intMaxTree, intArea.second);
    intUao.executeWithMSER(4, AltitudeDifference<std::int32_t>{1});
    requireEqual(intUao.getMaxContrastImage()->getSize(), intImage->getSize(), "typed int32 UAO MSER output must match image size");

    std::array<float, 16> floatPixels{
        0.3f, 0.3f, 0.2f, 0.2f, 0.3f, 0.4f, 0.4f, 0.2f, 0.1f, 0.4f, 0.5f, 0.2f, 0.1f, 0.1f, 0.5f, 0.0f,
    };
    auto floatImage = ImageFloat::fromExternal(floatPixels.data(), 4, 4);
    auto floatMaxTree = MorphologicalTreeFactory::createMaxTree(floatImage);
    MSERComputer<float> floatMser(floatMaxTree);
    requireEqual(static_cast<int>(floatMser.computeMSER(0.1f).size()), floatMaxTree.topology().numInternalNodeSlots(),
                 "typed float MSER output must be indexed by the owner tree node slots");
    auto floatArea = AttributeComputation::computeSingleAttribute(floatMaxTree, Area);
    UltimateAttributeOpening<float> floatUao(floatMaxTree, floatArea.second);
    floatUao.executeWithMSER(4, 0.1f);
    requireEqual(floatUao.getMaxContrastImage()->getSize(), floatImage->getSize(), "typed float UAO MSER output must match image size");

    adjust::CasfComponentTrees<std::uint8_t> casf(image, adjust::CasfComponentTreesAttribute::Area);
    require(&casf.minTree() != &casf.maxTree(), "CASF must own distinct min/max tree states");
    require(casf.minTree().topology().kind() == MorphologicalTreeKind::MinTree, "CASF min-tree accessor must expose the owner state");
    require(casf.maxTree().topology().kind() == MorphologicalTreeKind::MaxTree, "CASF max-tree accessor must expose the owner state");

    return 0;
}
