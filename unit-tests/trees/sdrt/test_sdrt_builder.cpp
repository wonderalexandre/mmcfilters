#include "../../support/TestSupport.hpp"

#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "mmcfilters/trees/adjust/DualMinMaxTreeIncrementalFilterLeaf.hpp"
#include "mmcfilters/trees/sdrt/SelfDualResidualTreeBuilder.hpp"

#include <cstdint>
#include <iostream>
#include <initializer_list>
#include <span>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace mmcfilters;
using namespace mmcfilters::sdrt;
using namespace mmcfilters::unit_tests;

namespace sdrt_detail = mmcfilters::sdrt::detail;

namespace {

void requireValidSdrtTopology(const MorphologicalTree& tree, const std::string& label) {
    tree.validateConnectedRootedTree();
    require(tree.getTreeType() == MorphologicalTreeKind::SELF_DUAL_RESIDUAL_TREE, label + " tree type");
    require(tree.getRoot() != InvalidNode, label + " root");
}

template<AltitudeValue T>
sdrt_detail::ResidualPolarity polarityOf(const WeightedMorphologicalTree<T>& tree, NodeId nodeId) {
    const auto& topology = tree.topology();
    require(!topology.isRoot(nodeId), "SDRT polarity requires a non-root node");
    const NodeId parentId = topology.getNodeParent(nodeId);
    const AltitudeDiff<T> signedDelta =
        static_cast<AltitudeDiff<T>>(tree.getAltitude(nodeId)) -
        static_cast<AltitudeDiff<T>>(tree.getAltitude(parentId));
    require(signedDelta != AltitudeDiff<T>{}, "SDRT test event must be non-flat");
    return signedDelta > AltitudeDiff<T>{} ? sdrt_detail::ResidualPolarity::Max : sdrt_detail::ResidualPolarity::Min;
}

template<AltitudeValue T>
T absorptionOf(const WeightedMorphologicalTree<T>& tree, NodeId nodeId) {
    const auto& topology = tree.topology();
    require(!topology.isRoot(nodeId), "SDRT absorption requires a non-root node");
    return tree.getAltitude(topology.getNodeParent(nodeId));
}

int supportSizeOf(const MorphologicalTree& tree, NodeId nodeId) {
    int area = 0;
    for (NodeId pixelId : tree.getConnectedComponent(nodeId)) {
        static_cast<void>(pixelId);
        ++area;
    }
    return area;
}

template<class T>
std::vector<T> toVector(std::span<const T> values) {
    return std::vector<T>(values.begin(), values.end());
}

template<AltitudeValue T>
ImagePtr<T> makeTypedImage(int rows, int cols, std::initializer_list<T> values) {
    requireEqual(static_cast<int>(values.size()), rows * cols, "typed image buffer size");
    auto image = Image<T>::create(rows, cols);
    int index = 0;
    for (T value : values) {
        (*image)[index++] = value;
    }
    return image;
}

template<AltitudeValue T>
void buildWithFactorySeeds(
    SelfDualResidualTreeBuilder<T>& builder,
    const ImagePtr<T>& image,
    double radius = 1.5) {
    auto minTree = MorphologicalTreeFactory::createMinTree(image, radius);
    auto maxTree = MorphologicalTreeFactory::createMaxTree(image, radius);
    builder.build(image, std::move(minTree), std::move(maxTree));
}

static_assert(std::is_same_v<
    decltype(std::declval<SelfDualResidualTreeBuilder<std::int32_t>&>().getAltitude()),
    std::span<const std::int32_t>>);
static_assert(std::is_same_v<
    decltype(MorphologicalTreeFactory::createSelfDualResidualTree(std::declval<ImageInt32Ptr>(), 1.5)),
    WeightedMorphologicalTree<std::int32_t>>);
static_assert(std::is_same_v<
    decltype(MorphologicalTreeFactory::createSelfDualResidualTree(std::declval<ImageFloatPtr>(), 1.5)),
    WeightedMorphologicalTree<float>>);

void test_builder_accessors_require_completed_build() {
    SelfDualResidualTreeBuilder<std::uint8_t> builder;
    requireThrows<std::logic_error>(
        [&]() {
            static_cast<void>(builder.getRoot());
        },
        "SDRT builder access before build");

    const auto image = makeImage(1, 1, {3});
    buildWithFactorySeeds(builder, image);
    requireEqual(builder.getRows(), 1, "builder accessor rows");
    requireEqual(builder.getCols(), 1, "builder accessor cols");
    requireEqual(builder.getRoot(), 0, "builder accessor root");
    requireVectorEqual(toVector(builder.getNodeParent()), std::vector<NodeId>({0}), "builder accessor parent");
    requireVectorEqual(toVector(builder.getProperPartOwner()), std::vector<NodeId>({0}), "builder accessor owners");
    requireVectorEqual(toVector(builder.getAltitude()), std::vector<std::uint8_t>({3}), "builder accessor altitude");
}

void test_builder_constant_image() {
    const auto image = makeImage(1, 1, {42});
    SelfDualResidualTreeBuilder<std::uint8_t> builder;
    buildWithFactorySeeds(builder, image);
    requireEqual(builder.getRows(), 1, "constant SDRT rows");
    requireEqual(builder.getCols(), 1, "constant SDRT cols");
    requireEqual(static_cast<int>(builder.getNodeParent().size()), 1, "constant SDRT node count");
    requireEqual(static_cast<int>(builder.getProperPartOwner().size()), 1, "constant SDRT proper part count");
    requireEqual(static_cast<int>(builder.getAltitude().size()), 1, "constant SDRT altitude count");
    requireEqual(builder.getRoot(), 0, "constant SDRT root");
    requireEqual(builder.getNodeParent()[0], 0, "constant SDRT root parent");
    requireEqual(builder.getProperPartOwner()[0], 0, "constant SDRT pixel owner");
    requireEqual(builder.getAltitude()[0], 42, "constant SDRT root altitude");

    auto weighted = MorphologicalTreeFactory::createSelfDualResidualTree(image);
    requireValidSdrtTopology(weighted.topology(), "constant weighted SDRT");
    requireVectorEqual(
        collectImageValues(weighted.reconstructionImage()),
        collectImageValues(image),
        "constant SDRT reconstruction");
}

void test_leaf_uses_dense_backend_for_uint8_altitudes() {
    using Leaf = mmcfilters::adjust::DualMinMaxTreeIncrementalFilterLeaf<std::uint8_t>;
    require(Leaf::usesDenseLevelBackend(), "uint8 SDRT leaf dense backend");
    requireEqual(Leaf::denseLevelBackendMaxBits(), 8, "uint8 SDRT dense backend threshold");
}

void test_leaf_uses_sparse_backend_for_wide_or_real_altitudes() {
    using IntLeaf = mmcfilters::adjust::DualMinMaxTreeIncrementalFilterLeaf<std::int32_t>;
    using FloatLeaf = mmcfilters::adjust::DualMinMaxTreeIncrementalFilterLeaf<float>;
    require(!IntLeaf::usesDenseLevelBackend(), "int32 SDRT leaf sparse backend");
    require(!FloatLeaf::usesDenseLevelBackend(), "float SDRT leaf sparse backend");
}

void test_builder_small_nonconstant_image() {
    const auto image = makeImage(1, 2, {4, 9});
    auto weighted = MorphologicalTreeFactory::createSelfDualResidualTree(image);
    requireValidSdrtTopology(weighted.topology(), "small weighted SDRT");
    weighted.validateAltitudeBufferShape();
    weighted.validateMonotoneAltitude();
    requireEqual(weighted.topology().getNumTotalProperParts(), 2, "small SDRT proper parts");
    require(weighted.topology().getNumInternalNodeSlots() >= 2, "small SDRT internal node count");
    requireVectorEqual(
        collectImageValues(weighted.reconstructionImage()),
        collectImageValues(image),
        "small SDRT reconstruction");
}

void test_builder_bright_peak_matches_original_sdrt_case() {
    const auto image = makeImage(3, 3, {
        2, 2, 2,
        2, 9, 2,
        2, 2, 2,
    });
    auto weighted = MorphologicalTreeFactory::createSelfDualResidualTree(image, 1.0);

    requireValidSdrtTopology(weighted.topology(), "bright peak SDRT");
    requireEqual(weighted.topology().getNumInternalNodeSlots(), 2, "bright peak node count");
    require(polarityOf(weighted, 1) == sdrt_detail::ResidualPolarity::Max, "bright peak polarity");
    requireEqual(weighted.getAltitude(1), 9, "bright peak event altitude");
    requireEqual(absorptionOf(weighted, 1), 2, "bright peak absorption");
    requireEqual(supportSizeOf(weighted.topology(), 1), 1, "bright peak support size");
    requireEqual(weighted.topology().getProperPartOwner(4), 1, "bright peak owner");
    requireVectorEqual(
        collectImageValues(weighted.reconstructionImage()),
        collectImageValues(image),
        "bright peak reconstruction");
}

void test_builder_dark_peak_matches_original_sdrt_case() {
    const auto image = makeImage(3, 3, {
        7, 7, 7,
        7, 1, 7,
        7, 7, 7,
    });
    auto weighted = MorphologicalTreeFactory::createSelfDualResidualTree(image, 1.0);

    requireValidSdrtTopology(weighted.topology(), "dark peak SDRT");
    requireEqual(weighted.topology().getNumInternalNodeSlots(), 2, "dark peak node count");
    require(polarityOf(weighted, 1) == sdrt_detail::ResidualPolarity::Min, "dark peak polarity");
    requireEqual(weighted.getAltitude(1), 1, "dark peak event altitude");
    requireEqual(absorptionOf(weighted, 1), 7, "dark peak absorption");
    requireEqual(supportSizeOf(weighted.topology(), 1), 1, "dark peak support size");
    requireEqual(weighted.topology().getProperPartOwner(4), 1, "dark peak owner");
    requireVectorEqual(
        collectImageValues(weighted.reconstructionImage()),
        collectImageValues(image),
        "dark peak reconstruction");
}

void test_builder_fixed_tie_order_matches_original_sdrt_case() {
    const auto image = makeImage(3, 5, {
        2, 2, 2, 2, 2,
        2, 9, 2, 1, 2,
        2, 2, 2, 2, 2,
    });

    auto weighted = MorphologicalTreeFactory::createSelfDualResidualTree(image, 1.0);
    requireValidSdrtTopology(weighted.topology(), "fixed-tie SDRT");
    requireEqual(weighted.topology().getNumInternalNodeSlots(), 3, "fixed-tie node count");
    require(polarityOf(weighted, 1) == sdrt_detail::ResidualPolarity::Max, "fixed-tie first polarity");
    require(polarityOf(weighted, 2) == sdrt_detail::ResidualPolarity::Min, "fixed-tie second polarity");
    requireEqual(weighted.topology().getProperPartOwner(6), 1, "fixed-tie bright owner");
    requireEqual(weighted.topology().getProperPartOwner(8), 2, "fixed-tie dark owner");
    requireVectorEqual(
        collectImageValues(weighted.reconstructionImage()),
        collectImageValues(image),
        "fixed-tie reconstruction");
}

void test_builder_two_level_bright_matches_original_sdrt_case() {
    const auto image = makeImage(3, 5, {
        2, 2, 2, 2, 2,
        2, 5, 5, 5, 2,
        2, 5, 9, 5, 2,
    });
    auto weighted = MorphologicalTreeFactory::createSelfDualResidualTree(image, 1.0);

    requireValidSdrtTopology(weighted.topology(), "two-level bright SDRT");
    requireEqual(weighted.topology().getNumInternalNodeSlots(), 3, "two-level bright node count");
    require(polarityOf(weighted, 1) == sdrt_detail::ResidualPolarity::Max, "two-level first polarity");
    require(polarityOf(weighted, 2) == sdrt_detail::ResidualPolarity::Max, "two-level second polarity");
    requireEqual(weighted.topology().getProperPartOwner(12), 1, "two-level peak owner");
    requireEqual(supportSizeOf(weighted.topology(), 1), 1, "two-level inner support size");
    requireEqual(supportSizeOf(weighted.topology(), 2), 6, "two-level outer support size");
    requireVectorEqual(
        collectImageValues(weighted.reconstructionImage()),
        collectImageValues(image),
        "two-level bright reconstruction");
}

void test_builder_dual_area_update_matches_original_sdrt_case() {
    const auto image = makeImage(1, 5, {4, 2, 2, 5, 5});
    auto weighted = MorphologicalTreeFactory::createSelfDualResidualTree(image, 1.0);

    requireValidSdrtTopology(weighted.topology(), "dual-area SDRT");
    requireEqual(weighted.topology().getNumInternalNodeSlots(), 3, "dual-area node count");
    require(polarityOf(weighted, 1) == sdrt_detail::ResidualPolarity::Max, "dual-area first polarity");
    require(polarityOf(weighted, 2) == sdrt_detail::ResidualPolarity::Max, "dual-area second polarity");
    requireEqual(weighted.getAltitude(1), 4, "dual-area first altitude");
    requireEqual(absorptionOf(weighted, 1), 2, "dual-area first absorption");
    requireEqual(weighted.getAltitude(2), 5, "dual-area second altitude");
    requireEqual(absorptionOf(weighted, 2), 2, "dual-area second absorption");
    requireEqual(supportSizeOf(weighted.topology(), 1), 1, "dual-area first support size");
    requireEqual(supportSizeOf(weighted.topology(), 2), 2, "dual-area second support size");
    requireEqual(weighted.topology().getProperPartOwner(0), 1, "dual-area left owner");
    requireEqual(weighted.topology().getProperPartOwner(3), 2, "dual-area right owner 3");
    requireEqual(weighted.topology().getProperPartOwner(4), 2, "dual-area right owner 4");
    requireVectorEqual(
        collectImageValues(weighted.reconstructionImage()),
        collectImageValues(image),
        "dual-area reconstruction");
}

void test_builder_int32_altitude_type_reconstructs_input() {
    const auto image = makeTypedImage<std::int32_t>(1, 3, {-1000, 0, 1000});
    SelfDualResidualTreeBuilder<std::int32_t> builder;
    buildWithFactorySeeds(builder, image, 1.0);

    requireEqual(builder.getRows(), 1, "int32 SDRT rows");
    requireEqual(builder.getCols(), 3, "int32 SDRT cols");
    requireEqual(static_cast<int>(builder.getProperPartOwner().size()), 3, "int32 SDRT proper part count");
    requireEqual(static_cast<int>(builder.getAltitude().size()), static_cast<int>(builder.getNodeParent().size()), "int32 SDRT altitude count");

    auto weighted = MorphologicalTreeFactory::createSelfDualResidualTree(image, 1.0);
    requireValidSdrtTopology(weighted.topology(), "int32 weighted SDRT");
    weighted.validateAltitudeBufferShape();
    requireVectorEqual(
        collectImageValues(weighted.reconstructionImage()),
        collectImageValues(image),
        "int32 SDRT reconstruction");
}

void test_builder_float_altitude_type_reconstructs_input() {
    const auto image = makeTypedImage<float>(1, 3, {-1.0f, 0.5f, 2.5f});
    SelfDualResidualTreeBuilder<float> builder;
    buildWithFactorySeeds(builder, image, 1.0);

    requireEqual(builder.getRows(), 1, "float SDRT rows");
    requireEqual(builder.getCols(), 3, "float SDRT cols");
    requireEqual(static_cast<int>(builder.getProperPartOwner().size()), 3, "float SDRT proper part count");
    requireEqual(static_cast<int>(builder.getAltitude().size()), static_cast<int>(builder.getNodeParent().size()), "float SDRT altitude count");

    auto weighted = MorphologicalTreeFactory::createSelfDualResidualTree(image, 1.0);
    requireValidSdrtTopology(weighted.topology(), "float weighted SDRT");
    weighted.validateAltitudeBufferShape();
    requireVectorEqual(
        collectImageValues(weighted.reconstructionImage()),
        collectImageValues(image),
        "float SDRT reconstruction");
}

} // namespace

int main() {
    try {
        test_builder_accessors_require_completed_build();
        test_leaf_uses_dense_backend_for_uint8_altitudes();
        test_leaf_uses_sparse_backend_for_wide_or_real_altitudes();
        test_builder_constant_image();
        test_builder_small_nonconstant_image();
        test_builder_bright_peak_matches_original_sdrt_case();
        test_builder_dark_peak_matches_original_sdrt_case();
        test_builder_fixed_tie_order_matches_original_sdrt_case();
        test_builder_two_level_bright_matches_original_sdrt_case();
        test_builder_dual_area_update_matches_original_sdrt_case();
        test_builder_int32_altitude_type_reconstructs_input();
        test_builder_float_altitude_type_reconstructs_input();
    } catch (const std::exception& ex) {
        std::cerr << "SDRT builder test failed: " << ex.what() << "\n";
        return 1;
    }
    return 0;
}
