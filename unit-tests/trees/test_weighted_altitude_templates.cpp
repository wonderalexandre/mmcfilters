#include "support/TestSupport.hpp"

#include "mmcfilters/trees/detail/TreeAltitudeDeltaNeighborhood.hpp"
#include "mmcfilters/attributes/computers/GrayLevelStatsComputer.hpp"
#include "mmcfilters/attributes/computers/VolumeComputer.hpp"
#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "mmcfilters/trees/TreeAltitudeAlgorithms.hpp"
#include "mmcfilters/trees/WeightedMorphologicalTree.hpp"
#include "mmcfilters/trees/WeightedTreeView.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

template class mmcfilters::WeightedMorphologicalTree<float>;
template class mmcfilters::WeightedTreeEditor<float>;

static_assert(std::is_same_v<typename WeightedMorphologicalTree<float>::altitude_type, float>);
static_assert(std::is_same_v<decltype(std::declval<const WeightedMorphologicalTree<float>&>().getAltitudeBuffer()), const std::vector<float>&>);
static_assert(std::is_same_v<decltype(std::declval<const WeightedMorphologicalTree<float>&>().asView()), WeightedTreeView<float>>);
static_assert(std::is_constructible_v<WeightedTreeView<std::int32_t>, const MorphologicalTree&, const AltitudeBuffer<std::int32_t>&>);
static_assert(!std::is_constructible_v<WeightedTreeView<std::int32_t>, const MorphologicalTree&, AltitudeBuffer<std::int32_t>&&>);
static_assert(std::is_same_v<decltype(MorphologicalTreeFactory::createMaxTree(std::declval<ImageInt32Ptr>(), 1.5)), WeightedMorphologicalTree<std::int32_t>>);
static_assert(std::is_same_v<decltype(MorphologicalTreeFactory::createMinTree(std::declval<ImageFloatPtr>(), 1.5)), WeightedMorphologicalTree<float>>);
static_assert(std::is_same_v<decltype(MorphologicalTreeFactory::createFromHigraParent<double>(
                                 std::declval<std::span<const NodeId>>(), std::declval<std::span<const double>>(), 4, 4, MorphologicalTreeKind::MAX_TREE,
                                 std::declval<std::optional<RegularGridAdjacency2D>>())),
                             WeightedMorphologicalTree<double>>);
static_assert(std::is_same_v<decltype(std::declval<const WeightedMorphologicalTree<std::int32_t>&>().reconstructionImage()), ImageInt32Ptr>);
static_assert(std::is_same_v<decltype(std::declval<const WeightedMorphologicalTree<float>&>().reconstructionImage()), ImageFloatPtr>);
static_assert(
    std::is_same_v<decltype(AttributeComputation::projectNodeValuesToExportedHigra(
                       std::declval<const WeightedMorphologicalTree<float>&>(), std::declval<const AttributeNames&>(), std::declval<std::span<const float>>())),
                   std::vector<float>>);
static_assert(std::is_same_v<decltype(AttributeComputation::projectNodeValuesToExportedHigra<double>(std::declval<const WeightedMorphologicalTree<float>&>(),
                                                                                                     std::declval<const AttributeNames&>(),
                                                                                                     std::declval<std::span<const double>>())),
                             std::vector<double>>);
static_assert(std::is_same_v<decltype(AttributeComputation::computeAttributes(std::declval<const WeightedMorphologicalTree<std::uint8_t>&>(),
                                                                              std::declval<const std::vector<AttributeOrGroup>&>())),
                             ComputedAttributeData<float>>);
static_assert(std::is_same_v<decltype(AttributeComputation::computeAttributes<double>(std::declval<const WeightedMorphologicalTree<std::uint8_t>&>(),
                                                                                      std::declval<const std::vector<AttributeOrGroup>&>())),
                             ComputedAttributeData<double>>);
static_assert(std::is_same_v<decltype(AttributeComputation::computeTopologyAttributes<double>(std::declval<const MorphologicalTree&>(),
                                                                                              std::declval<const std::vector<AttributeOrGroup>&>())),
                             ComputedAttributeData<double>>);
static_assert(std::is_same_v<decltype(AttributeComputation::computeSingleAttributeWithDelta(std::declval<const WeightedMorphologicalTree<std::int32_t>&>(),
                                                                                            AREA, std::declval<AltitudeDiff<std::int32_t>>(), 2)),
                             ComputedAttributeDataWithDelta<float>>);
static_assert(std::is_same_v<decltype(AttributeComputation::computeSingleAttributeWithDelta(std::declval<const WeightedMorphologicalTree<float>&>(), LEVEL,
                                                                                            std::declval<AltitudeDiff<float>>(), 2)),
                             ComputedAttributeDataWithDelta<float>>);
static_assert(std::is_same_v<decltype(AttributeComputation::computeSingleAttributeWithDelta<double>(
                                 std::declval<const WeightedMorphologicalTree<std::int32_t>&>(), AREA, std::declval<AltitudeDiff<std::int32_t>>(), 2)),
                             ComputedAttributeDataWithDelta<double>>);
static_assert(
    std::is_same_v<decltype(AttributeComputation::computeAttributeMapping<double>(std::declval<const WeightedMorphologicalTree<std::uint8_t>&>(), LEVEL)),
                   ImagePtr<double>>);

template <class T> std::vector<T> makeGenericAltitude(const WeightedMorphologicalTree<std::uint8_t>& weighted) {
    std::vector<T> altitude(static_cast<std::size_t>(weighted.topology().getNumInternalNodeSlots()), T{});

    for (NodeId nodeId : weighted.topology().getAliveNodeIds()) {
        const std::uint8_t base = weighted.getAltitude(nodeId);
        if constexpr (std::is_floating_point_v<T>) {
            altitude[static_cast<std::size_t>(nodeId)] = static_cast<T>(base) + static_cast<T>(0.25);
        } else {
            altitude[static_cast<std::size_t>(nodeId)] = static_cast<T>(base);
        }
    }

    return altitude;
}

template <class T> std::vector<T> makeEquivalentAltitude(const WeightedMorphologicalTree<std::uint8_t>& weighted) {
    std::vector<T> altitude(static_cast<std::size_t>(weighted.topology().getNumInternalNodeSlots()), T{});

    for (NodeId nodeId : weighted.topology().getAliveNodeIds()) {
        altitude[static_cast<std::size_t>(nodeId)] = static_cast<T>(weighted.getAltitude(nodeId));
    }

    return altitude;
}

NodeId findNonRootResidueSample(const WeightedMorphologicalTree<std::uint8_t>& weighted) {
    const auto& tree = weighted.topology();
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        if (tree.isRoot(nodeId)) {
            continue;
        }
        const NodeId parentNodeId = tree.getNodeParent(nodeId);
        if (weighted.getAltitude(nodeId) != weighted.getAltitude(parentNodeId)) {
            return nodeId;
        }
    }
    throw std::runtime_error("fixture must contain a non-root node with non-zero residue");
}

template <class T> void checkGenericStaticAltitudeAccess(const WeightedMorphologicalTree<std::uint8_t>& weighted, NodeId sampleNodeId) {
    const auto& tree = weighted.topology();
    const std::vector<T> altitude = makeGenericAltitude<T>(weighted);
    const std::span<const T> view(altitude);

    static_assert(std::is_same_v<decltype(TreeAltitudeAlgorithms::getAltitude(view, sampleNodeId)), T>);

    TreeAltitudeAlgorithms::validateAltitudeBufferShape(tree, view);
    require(TreeAltitudeAlgorithms::getAltitude(view, sampleNodeId) == altitude[static_cast<std::size_t>(sampleNodeId)],
            "templated getAltitude must preserve the altitude value type");

    if constexpr (contract::validationsEnabled) {
        requireThrows<std::runtime_error>(
            [&]() {
                const std::vector<T> wrongSize(static_cast<std::size_t>(tree.getNumInternalNodeSlots() - 1), T{});
                TreeAltitudeAlgorithms::validateAltitudeBufferShape(tree, std::span<const T>(wrongSize));
            },
            "templated validateAltitudeBufferShape must reject wrong size");

        requireThrows<std::invalid_argument>([&]() { static_cast<void>(TreeAltitudeAlgorithms::getAltitude(view, InvalidNode)); },
                                             "templated getAltitude must reject invalid node ids");
    }
}

template <class T> void checkWeightedTreeViewContract(const WeightedMorphologicalTree<std::uint8_t>& weighted, NodeId sampleNodeId, const std::string& label) {
    const auto& tree = weighted.topology();
    const std::vector<T> altitude = makeGenericAltitude<T>(weighted);
    const WeightedTreeView<T> view(tree, std::span<const T>(altitude));
    auto inferredSpanView = WeightedTreeView(tree, std::span<const T>(altitude));
    auto inferredBufferView = WeightedTreeView(tree, altitude);

    static_assert(std::is_same_v<typename WeightedTreeView<T>::altitude_type, T>);
    static_assert(std::is_same_v<decltype(inferredSpanView), WeightedTreeView<T>>);
    static_assert(std::is_same_v<decltype(inferredBufferView), WeightedTreeView<T>>);

    require(&view.topology() == &tree, label + " view must reference the original topology");
    require(&inferredSpanView.topology() == &tree, label + " span view must reference the original topology");
    require(&inferredBufferView.topology() == &tree, label + " buffer view must reference the original topology");
    requireEqual(static_cast<int>(view.altitude().size()), tree.getNumInternalNodeSlots(), label + " view altitude span size");
    requireEqual(static_cast<int>(inferredSpanView.altitude().size()), tree.getNumInternalNodeSlots(), label + " inferred span view altitude span size");
    requireEqual(static_cast<int>(inferredBufferView.altitude().size()), tree.getNumInternalNodeSlots(), label + " inferred buffer view altitude span size");
    require(view.getAltitude(sampleNodeId) == altitude[static_cast<std::size_t>(sampleNodeId)],
            label + " view getAltitude must preserve the altitude value type");
    require(inferredSpanView.getAltitude(sampleNodeId) == altitude[static_cast<std::size_t>(sampleNodeId)],
            label + " inferred span view getAltitude must preserve the altitude value type");
    require(inferredBufferView.getAltitude(sampleNodeId) == altitude[static_cast<std::size_t>(sampleNodeId)],
            label + " inferred buffer view getAltitude must preserve the altitude value type");

    const auto expectedResidue = TreeAltitudeAlgorithms::getNodeResidue(tree, std::span<const T>(altitude), sampleNodeId);
    if constexpr (std::is_floating_point_v<T>) {
        requireNear(view.getNodeResidue(sampleNodeId), expectedResidue, static_cast<T>(1.0e-6), label + " view residue must match static generic helper");
    } else {
        requireEqual(view.getNodeResidue(sampleNodeId), expectedResidue, label + " view residue must match static generic helper");
    }

    const auto fromView = AttributeComputation::computeAttributesFromAltitudeSpan(view, {LEVEL, GRAY_HEIGHT});
    const auto fromSpan = AttributeComputation::computeAttributesFromAltitudeSpan(inferredSpanView, {LEVEL, GRAY_HEIGHT});
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        for (Attribute attribute : {LEVEL, GRAY_HEIGHT}) {
            requireNear(fromView.second[fromView.first.linearIndex(nodeId, attribute)], fromSpan.second[fromSpan.first.linearIndex(nodeId, attribute)], 1.0e-5f,
                        label + " view attribute pipeline must match inferred span view");
        }
    }

    if constexpr (contract::validationsEnabled) {
        requireThrows<std::runtime_error>(
            [&]() {
                const std::vector<T> wrongSize(static_cast<std::size_t>(tree.getNumInternalNodeSlots() - 1), T{});
                static_cast<void>(WeightedTreeView<T>(tree, std::span<const T>(wrongSize)));
            },
            label + " view must reject wrong altitude size");

        requireThrows<std::invalid_argument>([&]() { static_cast<void>(view.getAltitude(InvalidNode)); }, label + " view must reject invalid node ids");
    }
}

void requireComputedAttributesNear(const ComputedAttributeData<float>& actual, const ComputedAttributeData<float>& expected, const MorphologicalTree& tree,
                                   std::initializer_list<Attribute> attributes, const std::string& label) {
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        for (Attribute attribute : attributes) {
            requireNear(actual.second[actual.first.linearIndex(nodeId, attribute)], expected.second[expected.first.linearIndex(nodeId, attribute)], 1.0e-5f,
                        label + " " + AttributeNames::toString(attribute) + " node " + std::to_string(nodeId));
        }
    }
}

void requireFloatEquivalent(float actual, float expected, const std::string& label) {
    if (std::isnan(expected)) {
        require(std::isnan(actual), label + " expected NaN");
        return;
    }
    if (std::isinf(expected)) {
        require(std::isinf(actual) && std::signbit(actual) == std::signbit(expected), label + " expected infinity");
        return;
    }
    requireNear(actual, expected, 1.0e-5f, label);
}

void requireComputedAttributesEquivalent(const ComputedAttributeData<float>& actual, const ComputedAttributeData<float>& expected,
                                         const MorphologicalTree& tree, const std::vector<Attribute>& attributes, const std::string& label) {
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        for (Attribute attribute : attributes) {
            requireFloatEquivalent(actual.second[actual.first.linearIndex(nodeId, attribute)], expected.second[expected.first.linearIndex(nodeId, attribute)],
                                   label + " " + AttributeNames::toString(attribute) + " node " + std::to_string(nodeId));
        }
    }
}

void requireDeltaAttributeNear(const ComputedAttributeDataWithDelta<float>& actual, const ComputedAttributeDataWithDelta<float>& expected,
                               const MorphologicalTree& tree, Attribute attribute, int delta, const std::string& label) {
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        for (int d = -delta; d <= delta; ++d) {
            requireNear(actual.second[actual.first.linearIndex(nodeId, attribute, d)], expected.second[expected.first.linearIndex(nodeId, attribute, d)],
                        1.0e-5f, label + " node " + std::to_string(nodeId) + " delta " + std::to_string(d));
        }
    }
}

void requireImageNear(const ImageFloatPtr& actual, const ImageFloatPtr& expected, const std::string& label) {
    requireEqual(actual->getSize(), expected->getSize(), label + " image size");
    for (int i = 0; i < actual->getSize(); ++i) {
        requireNear((*actual)[i], (*expected)[i], 1.0e-5f, label + " pixel " + std::to_string(i));
    }
}

template <class T> ImagePtr<T> makeTypedComponentTreeFixture() {
    auto image = Image<T>::create(4, 4);
    const std::array<int, 16> values{
        3, 3, 2, 2, 3, 4, 4, 2, 1, 4, 5, 2, 1, 1, 5, 0,
    };
    for (std::size_t i = 0; i < values.size(); ++i) {
        (*image)[static_cast<int>(i)] = static_cast<T>(values[i]);
    }
    return image;
}

template <class T> ImagePtr<T> makeTypedComponentTreeFixtureWithValue(std::size_t index, T value) {
    auto image = makeTypedComponentTreeFixture<T>();
    (*image)[static_cast<int>(index)] = value;
    return image;
}

template <class T>
void requireTypedOwnerMatchesCanonicalImageTree(const WeightedMorphologicalTree<T>& typed, const WeightedMorphologicalTree<std::uint8_t>& canonical,
                                                const ImageUInt8Ptr& expectedReconstruction, const std::string& label) {
    requireEqual(typed.topology().getNumInternalNodeSlots(), canonical.topology().getNumInternalNodeSlots(), label + " internal node slot count");
    requireEqual(typed.topology().getNumTotalProperParts(), canonical.topology().getNumTotalProperParts(), label + " proper part count");
    requireEqual(typed.topology().getRoot(), canonical.topology().getRoot(), label + " root");
    requireEqual(static_cast<int>(typed.getAltitudeBuffer().size()), typed.topology().getNumInternalNodeSlots(), label + " typed altitude buffer size");

    for (NodeId nodeId : typed.topology().getAliveNodeIds()) {
        requireNear(static_cast<double>(typed.getAltitude(nodeId)), static_cast<double>(canonical.getAltitude(nodeId)), 1.0e-6, label + " node altitude");
    }

    std::vector<T> expectedTypedReconstruction;
    expectedTypedReconstruction.reserve(static_cast<std::size_t>(expectedReconstruction->getSize()));
    for (auto value : collectImageValues(expectedReconstruction)) {
        expectedTypedReconstruction.push_back(static_cast<T>(value));
    }
    requireVectorEqual(collectImageValues(typed.reconstructionImage()), expectedTypedReconstruction, label + " typed reconstruction");
}

void checkTypedHigraImportFactory(ImageUInt8Ptr image) {
    auto canonicalMax = MorphologicalTreeFactory::createMaxTree(image, 1.5);

    const auto [higraParent, higraAltitude] = canonicalMax.exportHigraHierarchy();
    std::vector<float> floatHigraAltitude(higraAltitude.begin(), higraAltitude.end());
    auto importedFloat = MorphologicalTreeFactory::createFromHigraParent<float>(
        std::span<const NodeId>(higraParent), std::span<const float>(floatHigraAltitude), image->getNumRows(), image->getNumCols(),
        MorphologicalTreeKind::MAX_TREE, RegularGridAdjacency2D(image->getNumRows(), image->getNumCols(), 1.5));
    static_assert(std::is_same_v<decltype(importedFloat), WeightedMorphologicalTree<float>>);
    importedFloat.validateMonotoneAltitude();

    const auto [roundtripParent, roundtripAltitude] = importedFloat.exportHigraHierarchy();
    requireVectorEqual(roundtripParent, higraParent, "typed Higra factory parent roundtrip");
    requireEqual(roundtripAltitude.size(), floatHigraAltitude.size(), "typed Higra factory altitude roundtrip size");
    for (std::size_t i = 0; i < roundtripAltitude.size(); ++i) {
        requireNear(roundtripAltitude[i], floatHigraAltitude[i], 1.0e-6f, "typed Higra factory altitude roundtrip");
    }
    std::vector<float> expectedFloatReconstruction;
    expectedFloatReconstruction.reserve(static_cast<std::size_t>(image->getSize()));
    for (auto value : collectImageValues(image)) {
        expectedFloatReconstruction.push_back(static_cast<float>(value));
    }
    requireVectorEqual(collectImageValues(importedFloat.reconstructionImage()), expectedFloatReconstruction, "typed Higra factory reconstruction");

    requireThrows<std::invalid_argument>(
        [&]() {
            const std::vector<double> wrongAltitude(higraParent.size() - 1, 0.0);
            static_cast<void>(MorphologicalTreeFactory::createFromHigraParent<double>(
                std::span<const NodeId>(higraParent), std::span<const double>(wrongAltitude), image->getNumRows(), image->getNumCols(),
                MorphologicalTreeKind::MAX_TREE, RegularGridAdjacency2D(image->getNumRows(), image->getNumCols(), 1.5)));
        },
        "typed Higra factory must reject altitude size mismatch");
}

void checkFiniteFloatAltitudeValidation(ImageUInt8Ptr image) {
    requireThrows<std::invalid_argument>(
        [&]() {
            static_cast<void>(
                MorphologicalTreeFactory::createMaxTree(makeTypedComponentTreeFixtureWithValue<float>(3, std::numeric_limits<float>::quiet_NaN()), 1.5));
        },
        "Image<float> max-tree factory must reject NaN");

    requireThrows<std::invalid_argument>(
        [&]() {
            static_cast<void>(
                MorphologicalTreeFactory::createMinTree(makeTypedComponentTreeFixtureWithValue<float>(4, std::numeric_limits<float>::infinity()), 1.5));
        },
        "Image<float> min-tree factory must reject +inf");

    requireThrows<std::invalid_argument>(
        [&]() {
            static_cast<void>(
                MorphologicalTreeFactory::createMaxTree(makeTypedComponentTreeFixtureWithValue<float>(5, -std::numeric_limits<float>::infinity()), 1.5));
        },
        "Image<float> max-tree factory must reject -inf");

    auto floatMax = MorphologicalTreeFactory::createMaxTree(makeTypedComponentTreeFixture<float>(), 1.5);
    const NodeId sampleNodeId = floatMax.topology().getRoot();

    std::vector<float> invalidBuffer = floatMax.getAltitudeBuffer();
    invalidBuffer[static_cast<std::size_t>(sampleNodeId)] = std::numeric_limits<float>::quiet_NaN();
    requireThrows<std::invalid_argument>([&]() { floatMax.setAltitudeBuffer(std::move(invalidBuffer)); },
                                         "WeightedMorphologicalTree<float>::setAltitudeBuffer must reject NaN");

    requireThrows<std::invalid_argument>([&]() { floatMax.setAltitude(sampleNodeId, std::numeric_limits<float>::infinity()); },
                                         "WeightedMorphologicalTree<float>::setAltitude must reject +inf");

    {
        auto editor = floatMax.edit();
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(editor.createDetachedNode(-std::numeric_limits<float>::infinity())); },
                                             "WeightedTreeEditor<float>::createDetachedNode must reject -inf");
        editor.commit();
    }

    {
        auto editor = floatMax.edit();
        requireThrows<std::invalid_argument>([&]() { editor.setNodeAltitude(sampleNodeId, std::numeric_limits<float>::quiet_NaN()); },
                                             "WeightedTreeEditor<float>::setNodeAltitude must reject NaN");
        editor.commit();
    }

    const auto canonicalMax = MorphologicalTreeFactory::createMaxTree(image, 1.5);
    const auto [higraParent, higraAltitude] = canonicalMax.exportHigraHierarchy();
    std::vector<float> floatHigraAltitude(higraAltitude.begin(), higraAltitude.end());
    floatHigraAltitude.back() = std::numeric_limits<float>::infinity();
    requireThrows<std::invalid_argument>(
        [&]() {
            static_cast<void>(MorphologicalTreeFactory::createFromHigraParent<float>(
                std::span<const NodeId>(higraParent), std::span<const float>(floatHigraAltitude), image->getNumRows(), image->getNumCols(),
                MorphologicalTreeKind::MAX_TREE, RegularGridAdjacency2D(image->getNumRows(), image->getNumCols(), 1.5)));
        },
        "typed Higra import must reject non-finite float altitude");
}

void checkImageFactoriesDeduceAltitudeType(ImageUInt8Ptr image) {
    const auto canonicalMax = MorphologicalTreeFactory::createMaxTree(image, 1.5);
    const auto canonicalMin = MorphologicalTreeFactory::createMinTree(image, 1.5);

    auto int32Max = MorphologicalTreeFactory::createMaxTree(makeTypedComponentTreeFixture<std::int32_t>(), 1.5);
    static_assert(std::is_same_v<decltype(int32Max), WeightedMorphologicalTree<std::int32_t>>);
    int32Max.validateMonotoneAltitude();
    requireTypedOwnerMatchesCanonicalImageTree(int32Max, canonicalMax, image, "Image<int32_t> max-tree factory");
    auto int32Min = MorphologicalTreeFactory::createMinTree(makeTypedComponentTreeFixture<std::int32_t>(), 1.5);
    static_assert(std::is_same_v<decltype(int32Min), WeightedMorphologicalTree<std::int32_t>>);
    int32Min.validateMonotoneAltitude();
    requireTypedOwnerMatchesCanonicalImageTree(int32Min, canonicalMin, image, "Image<int32_t> min-tree factory");
    {
        const std::vector<AttributeOrGroup> requests{AttributeGroup::SHAPE, AttributeGroup::TREE_TOPOLOGY};
        const auto typedAttrs = AttributeComputation::computeAttributes(int32Max, requests);
        const auto canonicalAttrs = AttributeComputation::computeAttributes(canonicalMax, requests);
        std::vector<Attribute> expectedAttrs = ATTRIBUTE_GROUPS.at(AttributeGroup::SHAPE);
        const auto topologyAttrs = ATTRIBUTE_GROUPS.at(AttributeGroup::TREE_TOPOLOGY);
        expectedAttrs.insert(expectedAttrs.end(), topologyAttrs.begin(), topologyAttrs.end());
        requireComputedAttributesEquivalent(typedAttrs, canonicalAttrs, int32Max.topology(), expectedAttrs,
                                            "Image<int32_t> owner must compute semantic shape and topology groups");
        requireImageNear(AttributeComputation::computeAttributeMapping(int32Max, CONTOUR_PERIMETER),
                         AttributeComputation::computeAttributeMapping(canonicalMax, CONTOUR_PERIMETER),
                         "Image<int32_t> owner mapping overload must support contour attributes");
    }

    auto floatMax = MorphologicalTreeFactory::createMaxTree(makeTypedComponentTreeFixture<float>(), 1.5);
    static_assert(std::is_same_v<decltype(floatMax), WeightedMorphologicalTree<float>>);
    floatMax.validateMonotoneAltitude();
    requireTypedOwnerMatchesCanonicalImageTree(floatMax, canonicalMax, image, "Image<float> max-tree factory");

    auto floatMin = MorphologicalTreeFactory::createMinTree(makeTypedComponentTreeFixture<float>(), 1.5);
    static_assert(std::is_same_v<decltype(floatMin), WeightedMorphologicalTree<float>>);
    floatMin.validateMonotoneAltitude();
    requireTypedOwnerMatchesCanonicalImageTree(floatMin, canonicalMin, image, "Image<float> min-tree factory");
    {
        const std::vector<AttributeOrGroup> requests{LEVEL, AttributeGroup::MOMENTS};
        const auto typedAttrs = AttributeComputation::computeAttributes(floatMin, requests);
        const auto canonicalAttrs = AttributeComputation::computeAttributes(canonicalMin, requests);
        std::vector<Attribute> expectedAttrs{LEVEL};
        const auto momentAttrs = ATTRIBUTE_GROUPS.at(AttributeGroup::MOMENTS);
        expectedAttrs.insert(expectedAttrs.end(), momentAttrs.begin(), momentAttrs.end());
        requireComputedAttributesEquivalent(typedAttrs, canonicalAttrs, floatMin.topology(), expectedAttrs,
                                            "Image<float> owner must compute mixed generic attributes");

        const auto singleLevel = AttributeComputation::computeSingleAttribute(floatMin, LEVEL);
        const auto canonicalSingleLevel = AttributeComputation::computeSingleAttribute(canonicalMin, LEVEL);
        requireComputedAttributesEquivalent(singleLevel, canonicalSingleLevel, floatMin.topology(), {LEVEL},
                                            "Image<float> owner single-attribute overload must compute LEVEL");

        requireImageNear(AttributeComputation::computeAttributeMapping(floatMin, LEVEL), AttributeComputation::computeAttributeMapping(canonicalMin, LEVEL),
                         "Image<float> owner mapping overload must support LEVEL");
    }
}

template <class T> void checkTypedOwnerExportedHigraAttributeProjection(const std::string& label) {
    auto image = makeTypedComponentTreeFixture<T>();
    for (int i = 0; i < image->getSize(); ++i) {
        if constexpr (std::is_floating_point_v<T>) {
            (*image)[i] = static_cast<T>((*image)[i] + static_cast<T>(0.25));
        } else {
            (*image)[i] = static_cast<T>((*image)[i] + static_cast<T>(1000));
        }
    }

    auto typed = MorphologicalTreeFactory::createMaxTree(image, 1.5);
    const std::vector<AttributeOrGroup> requests{AREA,           VOLUME,      RELATIVE_VOLUME, LEVEL,       MEAN_LEVEL,
                                                 VARIANCE_LEVEL, GRAY_HEIGHT, MAX_DIST,        BOX_COL_MIN, BOX_ROW_MIN};
    const auto computed = AttributeComputation::computeAttributes(typed, requests);
    const auto projected = AttributeComputation::projectNodeValuesToExportedHigra(typed, computed.first, computed.second);

    const auto expectedSize = static_cast<std::size_t>(typed.topology().getNumTotalProperParts() + typed.topology().getNumNodes()) *
                              static_cast<std::size_t>(computed.first.NUM_ATTRIBUTES);
    requireEqual(projected.size(), expectedSize, label + " exported-Higra typed projection size");

    const NodeId sampleProperPart = 10;
    const NodeId sampleOwner = typed.topology().getProperPartOwner(sampleProperPart);
    const auto [sampleRow, sampleCol] = ImageUtils::to2D(sampleProperPart, typed.topology().getNumColsOfGridDomain2D());
    const float sampleAltitude = static_cast<float>(typed.getAltitude(sampleOwner));
    auto projectedAt = [&](Attribute attribute) { return projected[computed.first.linearIndex(sampleProperPart, attribute)]; };

    requireNear(projectedAt(AREA), 1.0f, 1.0e-5f, label + " unit AREA");
    requireNear(projectedAt(VOLUME), sampleAltitude, 1.0e-5f, label + " unit VOLUME must preserve typed altitude");
    requireNear(projectedAt(RELATIVE_VOLUME), 1.0f, 1.0e-5f, label + " unit RELATIVE_VOLUME");
    requireNear(projectedAt(LEVEL), sampleAltitude, 1.0e-5f, label + " unit LEVEL must preserve typed altitude");
    requireNear(projectedAt(MEAN_LEVEL), sampleAltitude, 1.0e-5f, label + " unit MEAN_LEVEL must preserve typed altitude");
    requireNear(projectedAt(VARIANCE_LEVEL), 0.0f, 1.0e-5f, label + " unit VARIANCE_LEVEL");
    requireNear(projectedAt(GRAY_HEIGHT), 0.0f, 1.0e-5f, label + " unit GRAY_HEIGHT");
    requireNear(projectedAt(MAX_DIST), 0.0f, 1.0e-5f, label + " unit MAX_DIST");
    requireNear(projectedAt(BOX_COL_MIN), static_cast<float>(sampleCol), 1.0e-5f, label + " unit BOX_COL_MIN");
    requireNear(projectedAt(BOX_ROW_MIN), static_cast<float>(sampleRow), 1.0e-5f, label + " unit BOX_ROW_MIN");
}

void checkWeightedTreeViewIncrementalAttributePipeline(const WeightedMorphologicalTree<std::uint8_t>& weighted) {
    const auto& tree = weighted.topology();
    const auto ownerView = weighted.asView();

    const auto singleBaseline = AttributeComputation::computeSingleAttribute(weighted, MAX_DIST);
    const auto singleFromView = AttributeComputation::computeSingleAttribute(ownerView, MAX_DIST);
    requireComputedAttributesNear(singleFromView, singleBaseline, tree, {MAX_DIST}, "owner view single attribute");

    const auto multiBaseline = AttributeComputation::computeAttributes(weighted, {AREA, MAX_DIST});
    const auto multiFromView = AttributeComputation::computeAttributes(ownerView, {AREA, MAX_DIST});
    requireComputedAttributesNear(multiFromView, multiBaseline, tree, {AREA, MAX_DIST}, "owner view multi attribute");

    requireImageNear(AttributeComputation::computeAttributeMapping(ownerView, LEVEL), AttributeComputation::computeAttributeMapping(weighted, LEVEL),
                     "owner view attribute mapping");

    const std::vector<std::uint8_t> externalAltitude = makeEquivalentAltitude<std::uint8_t>(weighted);
    const WeightedTreeView<std::uint8_t> externalView(tree, std::span<const std::uint8_t>(externalAltitude));

    const auto externalRegular = AttributeComputation::computeAttributes(externalView, {LEVEL, GRAY_HEIGHT});
    const auto externalFromAltitudeSpan = AttributeComputation::computeAttributesFromAltitudeSpan(externalView, {LEVEL, GRAY_HEIGHT});
    requireComputedAttributesNear(externalRegular, externalFromAltitudeSpan, tree, {LEVEL, GRAY_HEIGHT},
                                  "external view regular overload must use attribute-pipeline subset");

    const auto externalSingle = AttributeComputation::computeSingleAttribute(externalView, LEVEL);
    requireComputedAttributesNear(externalSingle, externalFromAltitudeSpan, tree, {LEVEL}, "external view single overload must use attribute-pipeline subset");

    const auto externalMaxDist = AttributeComputation::computeSingleAttribute(externalView, MAX_DIST);
    requireComputedAttributesNear(externalMaxDist, singleBaseline, tree, {MAX_DIST}, "external view attribute-pipeline path must support MAX_DIST");

    const auto externalMixed = AttributeComputation::computeAttributes(externalView, {LEVEL, MAX_DIST});
    const auto baselineMixed = AttributeComputation::computeAttributes(weighted, {LEVEL, MAX_DIST});
    requireComputedAttributesNear(externalMixed, baselineMixed, tree, {LEVEL, MAX_DIST},
                                  "external view attribute-pipeline path must support mixed MAX_DIST requests");

    requireImageNear(AttributeComputation::computeAttributeMapping(externalView, LEVEL), AttributeComputation::computeAttributeMapping(weighted, LEVEL),
                     "external view attribute-pipeline attribute mapping");
}

void checkFloatResidueAndMonotonicity(const WeightedMorphologicalTree<std::uint8_t>& weighted, NodeId sampleNodeId) {
    const auto& tree = weighted.topology();
    std::vector<float> altitude = makeGenericAltitude<float>(weighted);
    const std::span<const float> view(altitude);

    const NodeId parentNodeId = tree.getNodeParent(sampleNodeId);
    const float expectedResidue = altitude[static_cast<std::size_t>(sampleNodeId)] - altitude[static_cast<std::size_t>(parentNodeId)];
    requireNear(TreeAltitudeAlgorithms::getNodeResidue(tree, view, sampleNodeId), expectedResidue, 1.0e-6f,
                "templated getNodeResidue must preserve float arithmetic");

    TreeAltitudeAlgorithms::validateMonotoneAltitude(tree, view);

    std::vector<float> nonMonotone = altitude;
    nonMonotone[static_cast<std::size_t>(parentNodeId)] = nonMonotone[static_cast<std::size_t>(sampleNodeId)] + 1.0f;
    requireThrows<std::runtime_error>([&]() { TreeAltitudeAlgorithms::validateMonotoneAltitude(tree, std::span<const float>(nonMonotone)); },
                                      "templated validateMonotoneAltitude must reject non-monotone float altitude");
}

void checkUnsignedResidueDoesNotUnderflow(const WeightedMorphologicalTree<std::uint8_t>& weighted, NodeId sampleNodeId) {
    const auto& tree = weighted.topology();
    const NodeId parentNodeId = tree.getNodeParent(sampleNodeId);
    require(parentNodeId != InvalidNode && parentNodeId != sampleNodeId, "unsigned residue sample must be a non-root node");

    std::vector<std::uint32_t> altitude(static_cast<std::size_t>(tree.getNumInternalNodeSlots()), 0u);
    altitude[static_cast<std::size_t>(parentNodeId)] = 10u;
    altitude[static_cast<std::size_t>(sampleNodeId)] = 3u;

    const std::span<const std::uint32_t> view(altitude);
    const AltitudeDiff<std::uint32_t> expectedResidue = -7;
    requireEqual(TreeAltitudeAlgorithms::getNodeResidue(tree, view, sampleNodeId), expectedResidue,
                 "templated getNodeResidue must not underflow unsigned altitude subtraction");

    const WeightedTreeView<std::uint32_t> weightedView(tree, view);
    requireEqual(weightedView.getNodeResidue(sampleNodeId), expectedResidue,
                 "WeightedTreeView::getNodeResidue must not underflow unsigned altitude subtraction");
}

void checkFloatHigraExport(const WeightedMorphologicalTree<std::uint8_t>& weighted) {
    const auto& tree = weighted.topology();
    const std::vector<float> altitude = makeGenericAltitude<float>(weighted);
    const auto [parent, exportedAltitude] = TreeAltitudeAlgorithms::exportHigraHierarchy(tree, std::span<const float>(altitude));

    static_assert(std::is_same_v<std::remove_cvref_t<decltype(exportedAltitude)>, std::vector<float>>);

    requireEqual(parent.size(), exportedAltitude.size(), "templated Higra export parent/altitude size");
    requireEqual(static_cast<int>(parent.size()), tree.getNumTotalProperParts() + tree.getNumNodes(), "templated Higra export compact domain size");

    bool foundFractionalAltitude = false;
    for (float value : exportedAltitude) {
        if (std::abs((value - std::floor(value)) - 0.25f) < 1.0e-6f) {
            foundFractionalAltitude = true;
            break;
        }
    }
    require(foundFractionalAltitude, "templated Higra export must preserve fractional float altitudes");

    const NodeId numLeaves = tree.getNumTotalProperParts();
    for (NodeId leafId = 0; leafId < numLeaves; ++leafId) {
        const NodeId ownerHigraId = parent[static_cast<std::size_t>(leafId)];
        require(ownerHigraId >= numLeaves && ownerHigraId < static_cast<NodeId>(parent.size()),
                "templated Higra export leaf parent must be an internal Higra node");
        requireNear(exportedAltitude[static_cast<std::size_t>(leafId)], exportedAltitude[static_cast<std::size_t>(ownerHigraId)], 1.0e-6f,
                    "templated Higra export leaf altitude must match owner altitude");
    }
}

void checkInt32HigraExport(const WeightedMorphologicalTree<std::uint8_t>& weighted) {
    const auto& tree = weighted.topology();
    const std::vector<std::int32_t> altitude = makeGenericAltitude<std::int32_t>(weighted);
    const auto [parent, exportedAltitude] = TreeAltitudeAlgorithms::exportHigraHierarchy(tree, std::span<const std::int32_t>(altitude));

    static_assert(std::is_same_v<std::remove_cvref_t<decltype(exportedAltitude)>, std::vector<std::int32_t>>);

    requireEqual(parent.size(), exportedAltitude.size(), "int32 templated Higra export parent/altitude size");
    requireEqual(static_cast<int>(parent.size()), tree.getNumTotalProperParts() + tree.getNumNodes(), "int32 templated Higra export compact domain size");

    const NodeId numLeaves = tree.getNumTotalProperParts();
    for (NodeId leafId = 0; leafId < numLeaves; ++leafId) {
        const NodeId ownerHigraId = parent[static_cast<std::size_t>(leafId)];
        require(ownerHigraId >= numLeaves && ownerHigraId < static_cast<NodeId>(parent.size()),
                "int32 templated Higra export leaf parent must be an internal Higra node");
        requireEqual(exportedAltitude[static_cast<std::size_t>(leafId)], exportedAltitude[static_cast<std::size_t>(ownerHigraId)],
                     "int32 templated Higra export leaf altitude must match owner altitude");
    }
}

template <class T> void checkGenericAltitudeDelta(const WeightedMorphologicalTree<std::uint8_t>& weighted, const std::string& label) {
    const auto& tree = weighted.topology();
    const std::vector<T> altitude = makeGenericAltitude<T>(weighted);
    const std::span<const T> view(altitude);
    const auto baseline = detail::computeAscendantsAndDescendantsByAltitude(tree, std::span<const std::uint8_t>(weighted.getAltitudeBuffer()),
                                                                            static_cast<AltitudeDiff<std::uint8_t>>(2));
    const auto fromAltitudeSpan = detail::computeAscendantsAndDescendantsByAltitude(tree, view, static_cast<AltitudeDiff<T>>(2));

    requireVectorEqual(fromAltitudeSpan.first, baseline.first, label + " generic delta ascendants");
    requireVectorEqual(fromAltitudeSpan.second, baseline.second, label + " generic delta descendants");

    NodeId sampleNodeId = InvalidNode;
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        if (baseline.first[static_cast<std::size_t>(nodeId)] != InvalidNode) {
            sampleNodeId = nodeId;
            break;
        }
    }
    require(sampleNodeId != InvalidNode, label + " generic delta fixture needs one ascendant sample");
    requireEqual(detail::findAscendantByAltitudeDelta(tree, view, sampleNodeId, static_cast<AltitudeDiff<T>>(2)),
                 baseline.first[static_cast<std::size_t>(sampleNodeId)], label + " generic findAscendantByAltitudeDelta");

    if constexpr (contract::validationsEnabled) {
        requireThrows<std::runtime_error>(
            [&]() {
                const std::vector<T> wrongSize(static_cast<std::size_t>(tree.getNumInternalNodeSlots() - 1), T{});
                static_cast<void>(detail::computeAscendantsAndDescendantsByAltitude(tree, std::span<const T>(wrongSize),
                                                                                     static_cast<AltitudeDiff<T>>(2)));
            },
            label + " generic delta must reject wrong altitude size");

        requireThrows<std::invalid_argument>(
            [&]() { static_cast<void>(detail::findAscendantByAltitudeDelta(tree, view, InvalidNode, static_cast<AltitudeDiff<T>>(2))); },
            label + " generic delta must reject invalid node ids");
    }
}

template <class T> void checkTypedOwnerDeltaAttributeApi(const std::string& label) {
    auto canonicalImage = makeTypedComponentTreeFixture<std::uint8_t>();
    auto canonical = MorphologicalTreeFactory::createMaxTree(canonicalImage, 1.5);

    auto typedImage = makeTypedComponentTreeFixture<T>();
    AltitudeDiff<T> deltaStep{};
    if constexpr (std::is_floating_point_v<T>) {
        for (int i = 0; i < typedImage->getSize(); ++i) {
            (*typedImage)[i] = static_cast<T>((*typedImage)[i] * static_cast<T>(0.25));
        }
        deltaStep = static_cast<AltitudeDiff<T>>(0.25);
    } else {
        for (int i = 0; i < typedImage->getSize(); ++i) {
            (*typedImage)[i] = static_cast<T>((*typedImage)[i] + static_cast<T>(1000));
        }
        deltaStep = static_cast<AltitudeDiff<T>>(1);
    }

    auto typed = MorphologicalTreeFactory::createMaxTree(typedImage, 1.5);
    auto typedDelta = AttributeComputation::computeSingleAttributeWithDelta(typed, AREA, deltaStep, 2, "last-padding");
    auto canonicalDelta = AttributeComputation::computeSingleAttributeWithDelta(canonical, AREA, AltitudeDiff<std::uint8_t>{1}, 2, "last-padding");
    requireDeltaAttributeNear(typedDelta, canonicalDelta, typed.topology(), AREA, 2, label + " typed delta AREA must match equivalent canonical layout");

    auto typedLevelDelta = AttributeComputation::computeSingleAttributeWithDelta(typed, LEVEL, deltaStep, 1, "null-padding");
    bool foundAscendantSample = false;
    const auto ascendantAndDescendant = detail::computeAscendantsAndDescendantsByAltitude(typed.topology(), typed.altitudeSpan(), deltaStep);
    const auto& ascendants = ascendantAndDescendant.first;
    for (NodeId nodeId : typed.topology().getAliveNodeIds()) {
        requireNear(typedLevelDelta.second[typedLevelDelta.first.linearIndex(nodeId, LEVEL, 0)], static_cast<float>(typed.getAltitude(nodeId)), 1.0e-5f,
                    label + " typed delta LEVEL center must preserve typed altitude");

        const NodeId ascendant = ascendants[static_cast<std::size_t>(nodeId)];
        if (ascendant != InvalidNode && ascendant != nodeId) {
            requireNear(typedLevelDelta.second[typedLevelDelta.first.linearIndex(nodeId, LEVEL, -1)], static_cast<float>(typed.getAltitude(ascendant)), 1.0e-5f,
                        label + " typed delta LEVEL ascendant must use typed altitude distance");
            foundAscendantSample = true;
            break;
        }
    }
    require(foundAscendantSample, label + " typed delta fixture must contain an ascendant sample");

    if constexpr (contract::validationsEnabled) {
        requireThrows<std::invalid_argument>(
            [&]() { static_cast<void>(AttributeComputation::computeSingleAttributeWithDelta(typed, AREA, deltaStep, -1, "last-padding")); },
            label + " typed delta must reject negative radius");

        if constexpr (std::is_floating_point_v<T>) {
            requireThrows<std::invalid_argument>(
                [&]() {
                    static_cast<void>(AttributeComputation::computeSingleAttributeWithDelta(
                        typed, AREA, std::numeric_limits<AltitudeDiff<T>>::quiet_NaN(), 1, "last-padding"));
                },
                label + " typed delta must reject non-finite floating delta step");
        }
    }
}

template <class T> void checkGenericVolumeKernel(const WeightedMorphologicalTree<std::uint8_t>& weighted, const std::string& label) {
    const auto& tree = weighted.topology();
    const auto areaComputed = AttributeComputation::computeSingleAttribute(weighted, AREA);
    const auto baseline = AttributeComputation::computeAttributes(weighted, {VOLUME, RELATIVE_VOLUME});

    const AttributeNames volumeNames = AttributeNames::fromList({VOLUME, RELATIVE_VOLUME});
    const std::vector<Attribute> requested{VOLUME, RELATIVE_VOLUME};
    const std::array<DependencySourceT<float>, 1> dependencies{{DependencySourceT<float>{&areaComputed.first, areaComputed.second.data()}}};
    std::vector<float> genericBuffer(static_cast<std::size_t>(tree.getNumInternalNodeSlots()) * static_cast<std::size_t>(volumeNames.NUM_ATTRIBUTES), 0.0f);

    const std::vector<T> equivalentAltitude = makeEquivalentAltitude<T>(weighted);
    const auto equivalentContext = AltitudeAttributeComputeContext<float, T>{tree, std::span<const T>(equivalentAltitude), std::span<float>(genericBuffer),
                                                                             volumeNames, requested,
                                                                             std::span<const DependencySourceT<float>>(dependencies)};
    const auto volumeRequest = mmcfilters::attributes::computers::detail::VolumeRequest::from(requested);
    mmcfilters::attributes::computers::detail::kernel::computeVolume(equivalentContext, volumeRequest, &dependencies[0]);

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        requireNear(genericBuffer[volumeNames.linearIndex(nodeId, VOLUME)], baseline.second[baseline.first.linearIndex(nodeId, VOLUME)], 1.0e-5f,
                    label + " generic VOLUME must match current path");
        requireNear(genericBuffer[volumeNames.linearIndex(nodeId, RELATIVE_VOLUME)], baseline.second[baseline.first.linearIndex(nodeId, RELATIVE_VOLUME)],
                    1.0e-5f, label + " generic RELATIVE_VOLUME must match current path");
    }

    if constexpr (std::is_floating_point_v<T>) {
        std::vector<float> fractionalBuffer(genericBuffer.size(), 0.0f);
        const std::vector<T> fractionalAltitude = makeGenericAltitude<T>(weighted);
        const auto fractionalContext = AltitudeAttributeComputeContext<float, T>{
            tree, std::span<const T>(fractionalAltitude), std::span<float>(fractionalBuffer), volumeNames, requested,
            std::span<const DependencySourceT<float>>(dependencies)};
        mmcfilters::attributes::computers::detail::kernel::computeVolume(fractionalContext, volumeRequest, &dependencies[0]);

        for (NodeId nodeId : tree.getAliveNodeIds()) {
            const float area = areaComputed.second[areaComputed.first.linearIndex(nodeId, AREA)];
            requireNear(fractionalBuffer[volumeNames.linearIndex(nodeId, VOLUME)], baseline.second[baseline.first.linearIndex(nodeId, VOLUME)] + (area * 0.25f),
                        1.0e-5f, label + " fractional VOLUME must preserve float altitude contribution");
            requireNear(fractionalBuffer[volumeNames.linearIndex(nodeId, RELATIVE_VOLUME)],
                        baseline.second[baseline.first.linearIndex(nodeId, RELATIVE_VOLUME)], 1.0e-5f,
                        label + " uniform fractional offset must preserve RELATIVE_VOLUME");
        }
    }
}

template <class T> void checkGenericGrayLevelStatsKernel(const WeightedMorphologicalTree<std::uint8_t>& weighted, const std::string& label) {
    const auto& tree = weighted.topology();
    const auto areaComputed = AttributeComputation::computeSingleAttribute(weighted, AREA);
    const auto volumeComputed = AttributeComputation::computeAttributes(weighted, {VOLUME, RELATIVE_VOLUME});
    const auto baseline = AttributeComputation::computeAttributes(weighted, {LEVEL, MEAN_LEVEL, VARIANCE_LEVEL, GRAY_HEIGHT});

    const AttributeNames grayNames = AttributeNames::fromList({LEVEL, MEAN_LEVEL, VARIANCE_LEVEL, GRAY_HEIGHT});
    const std::vector<Attribute> requested{LEVEL, MEAN_LEVEL, VARIANCE_LEVEL, GRAY_HEIGHT};
    const std::array<DependencySourceT<float>, 2> dependencies{{DependencySourceT<float>{&volumeComputed.first, volumeComputed.second.data()},
                                                                DependencySourceT<float>{&areaComputed.first, areaComputed.second.data()}}};
    std::vector<float> genericBuffer(static_cast<std::size_t>(tree.getNumInternalNodeSlots()) * static_cast<std::size_t>(grayNames.NUM_ATTRIBUTES), 0.0f);

    const std::vector<T> equivalentAltitude = makeEquivalentAltitude<T>(weighted);
    const auto equivalentContext = AltitudeAttributeComputeContext<float, T>{tree, std::span<const T>(equivalentAltitude), std::span<float>(genericBuffer),
                                                                             grayNames, requested,
                                                                             std::span<const DependencySourceT<float>>(dependencies)};
    const auto grayRequest = mmcfilters::attributes::computers::detail::GrayLevelStatsRequest::from(requested);
    mmcfilters::attributes::computers::detail::kernel::computeGrayLevelStats(equivalentContext, grayRequest, &dependencies[0], &dependencies[1]);

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        for (Attribute attribute : requested) {
            requireNear(genericBuffer[grayNames.linearIndex(nodeId, attribute)], baseline.second[baseline.first.linearIndex(nodeId, attribute)], 1.0e-5f,
                        label + " generic gray-level attribute must match current path");
        }
    }

    if constexpr (std::is_floating_point_v<T>) {
        const AttributeNames volumeNames = AttributeNames::fromList({VOLUME, RELATIVE_VOLUME});
        const std::vector<Attribute> volumeRequested{VOLUME, RELATIVE_VOLUME};
        const std::array<DependencySourceT<float>, 1> volumeDependencies{{DependencySourceT<float>{&areaComputed.first, areaComputed.second.data()}}};
        std::vector<float> fractionalVolumeBuffer(
            static_cast<std::size_t>(tree.getNumInternalNodeSlots()) * static_cast<std::size_t>(volumeNames.NUM_ATTRIBUTES), 0.0f);
        const std::vector<T> fractionalAltitude = makeGenericAltitude<T>(weighted);
        const auto fractionalVolumeContext = AltitudeAttributeComputeContext<float, T>{
            tree, std::span<const T>(fractionalAltitude), std::span<float>(fractionalVolumeBuffer), volumeNames, volumeRequested,
            std::span<const DependencySourceT<float>>(volumeDependencies)};
        const auto fractionalVolumeRequest = mmcfilters::attributes::computers::detail::VolumeRequest::from(volumeRequested);
        mmcfilters::attributes::computers::detail::kernel::computeVolume(fractionalVolumeContext, fractionalVolumeRequest, &volumeDependencies[0]);

        const std::array<DependencySourceT<float>, 2> fractionalDependencies{
            {DependencySourceT<float>{&volumeNames, fractionalVolumeBuffer.data()}, DependencySourceT<float>{&areaComputed.first, areaComputed.second.data()}}};
        std::vector<float> fractionalGrayBuffer(genericBuffer.size(), 0.0f);
        const auto fractionalGrayContext = AltitudeAttributeComputeContext<float, T>{
            tree, std::span<const T>(fractionalAltitude), std::span<float>(fractionalGrayBuffer), grayNames, requested,
            std::span<const DependencySourceT<float>>(fractionalDependencies)};
        mmcfilters::attributes::computers::detail::kernel::computeGrayLevelStats(fractionalGrayContext, grayRequest, &fractionalDependencies[0],
                                                                                 &fractionalDependencies[1]);

        for (NodeId nodeId : tree.getAliveNodeIds()) {
            requireNear(fractionalGrayBuffer[grayNames.linearIndex(nodeId, LEVEL)], baseline.second[baseline.first.linearIndex(nodeId, LEVEL)] + 0.25f, 1.0e-5f,
                        label + " fractional LEVEL must preserve float altitude");
            requireNear(fractionalGrayBuffer[grayNames.linearIndex(nodeId, MEAN_LEVEL)],
                        baseline.second[baseline.first.linearIndex(nodeId, MEAN_LEVEL)] + 0.25f, 1.0e-5f,
                        label + " fractional MEAN_LEVEL must preserve uniform altitude offset");
            requireNear(fractionalGrayBuffer[grayNames.linearIndex(nodeId, VARIANCE_LEVEL)],
                        baseline.second[baseline.first.linearIndex(nodeId, VARIANCE_LEVEL)], 1.0e-4f,
                        label + " uniform fractional offset must preserve VARIANCE_LEVEL");
            requireNear(fractionalGrayBuffer[grayNames.linearIndex(nodeId, GRAY_HEIGHT)], baseline.second[baseline.first.linearIndex(nodeId, GRAY_HEIGHT)],
                        1.0e-5f, label + " uniform fractional offset must preserve GRAY_HEIGHT");
        }
    }
}

template <class T> void checkAttributePipelineApi(const WeightedMorphologicalTree<std::uint8_t>& weighted, const std::string& label) {
    const auto& tree = weighted.topology();
    const std::vector<AttributeOrGroup> requests{AREA, VOLUME, RELATIVE_VOLUME, LEVEL, MEAN_LEVEL, VARIANCE_LEVEL, GRAY_HEIGHT, MAX_DIST};
    const std::vector<Attribute> expectedAttributes{AREA, VOLUME, RELATIVE_VOLUME, LEVEL, GRAY_HEIGHT, MEAN_LEVEL, VARIANCE_LEVEL, MAX_DIST};
    const auto baseline = AttributeComputation::computeAttributes(weighted, requests);
    const std::vector<T> equivalentAltitude = makeEquivalentAltitude<T>(weighted);
    const WeightedTreeView<T> equivalentView(tree, std::span<const T>(equivalentAltitude));
    const auto fromAltitudeSpan = AttributeComputation::computeAttributesFromAltitudeSpan(equivalentView, requests);
    const auto viewOverload = AttributeComputation::computeAttributes(equivalentView, requests);

    for (NodeId nodeId : tree.getAliveNodeIds()) {
        for (Attribute attribute : expectedAttributes) {
            requireNear(fromAltitudeSpan.second[fromAltitudeSpan.first.linearIndex(nodeId, attribute)],
                        baseline.second[baseline.first.linearIndex(nodeId, attribute)], 1.0e-5f,
                        label + " attribute pipeline must match current path for " + AttributeNames::toString(attribute));
            requireNear(viewOverload.second[viewOverload.first.linearIndex(nodeId, attribute)],
                        fromAltitudeSpan.second[fromAltitudeSpan.first.linearIndex(nodeId, attribute)], 1.0e-5f,
                        label + " altitude-span view overload must match equivalent view API for " + AttributeNames::toString(attribute));
        }
    }

    const auto singleViewOverload = AttributeComputation::computeSingleAttribute(equivalentView, LEVEL);
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        requireNear(singleViewOverload.second[singleViewOverload.first.linearIndex(nodeId, LEVEL)], baseline.second[baseline.first.linearIndex(nodeId, LEVEL)],
                    1.0e-5f, label + " altitude-span single view overload must match current LEVEL");
    }

    requireImageNear(AttributeComputation::computeAttributeMapping(equivalentView, LEVEL), AttributeComputation::computeAttributeMapping(weighted, LEVEL),
                     label + " altitude-span view mapping overload must match current LEVEL mapping");

    const auto grayLevelFromAltitudeSpan = AttributeComputation::computeAttributesFromAltitudeSpan(equivalentView, {AttributeGroup::GRAY_LEVEL});
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        for (Attribute attribute : {VOLUME, RELATIVE_VOLUME, LEVEL, GRAY_HEIGHT, MEAN_LEVEL, VARIANCE_LEVEL}) {
            requireNear(grayLevelFromAltitudeSpan.second[grayLevelFromAltitudeSpan.first.linearIndex(nodeId, attribute)],
                        baseline.second[baseline.first.linearIndex(nodeId, attribute)], 1.0e-5f,
                        label + " altitude-span GRAY_LEVEL group must match current path for " + AttributeNames::toString(attribute));
        }
    }

    if constexpr (std::is_floating_point_v<T>) {
        const auto areaComputed = AttributeComputation::computeSingleAttribute(weighted, AREA);
        const std::vector<T> fractionalAltitude = makeGenericAltitude<T>(weighted);
        const WeightedTreeView<T> fractionalView(tree, std::span<const T>(fractionalAltitude));
        const auto fractional = AttributeComputation::computeAttributesFromAltitudeSpan(fractionalView, requests);

        for (NodeId nodeId : tree.getAliveNodeIds()) {
            const float area = areaComputed.second[areaComputed.first.linearIndex(nodeId, AREA)];
            requireNear(fractional.second[fractional.first.linearIndex(nodeId, VOLUME)],
                        baseline.second[baseline.first.linearIndex(nodeId, VOLUME)] + (area * 0.25f), 1.0e-5f,
                        label + " altitude-span API fractional VOLUME must preserve altitude contribution");
            requireNear(fractional.second[fractional.first.linearIndex(nodeId, RELATIVE_VOLUME)],
                        baseline.second[baseline.first.linearIndex(nodeId, RELATIVE_VOLUME)], 1.0e-5f,
                        label + " altitude-span API uniform fractional offset must preserve RELATIVE_VOLUME");
            requireNear(fractional.second[fractional.first.linearIndex(nodeId, LEVEL)], baseline.second[baseline.first.linearIndex(nodeId, LEVEL)] + 0.25f,
                        1.0e-5f, label + " altitude-span API fractional LEVEL must preserve altitude");
            requireNear(fractional.second[fractional.first.linearIndex(nodeId, MEAN_LEVEL)],
                        baseline.second[baseline.first.linearIndex(nodeId, MEAN_LEVEL)] + 0.25f, 1.0e-5f,
                        label + " altitude-span API fractional MEAN_LEVEL must preserve altitude");
            requireNear(fractional.second[fractional.first.linearIndex(nodeId, VARIANCE_LEVEL)],
                        baseline.second[baseline.first.linearIndex(nodeId, VARIANCE_LEVEL)], 1.0e-4f,
                        label + " altitude-span API uniform fractional offset must preserve VARIANCE_LEVEL");
            requireNear(fractional.second[fractional.first.linearIndex(nodeId, GRAY_HEIGHT)], baseline.second[baseline.first.linearIndex(nodeId, GRAY_HEIGHT)],
                        1.0e-5f, label + " altitude-span API uniform fractional offset must preserve GRAY_HEIGHT");
            requireNear(fractional.second[fractional.first.linearIndex(nodeId, MAX_DIST)], baseline.second[baseline.first.linearIndex(nodeId, MAX_DIST)],
                        1.0e-5f, label + " altitude-span API uniform fractional offset must preserve MAX_DIST");
        }
    }
}

template <class T> void checkTopologyOnlyAttributePipelineApi(const WeightedMorphologicalTree<std::uint8_t>& weighted, const std::string& label) {
    const auto& tree = weighted.topology();
    const std::vector<T> equivalentAltitude = makeEquivalentAltitude<T>(weighted);
    const WeightedTreeView<T> equivalentView(tree, std::span<const T>(equivalentAltitude));

    for (AttributeGroup group : {AttributeGroup::GRAY_LEVEL, AttributeGroup::SHAPE, AttributeGroup::MOMENTS, AttributeGroup::BOUNDARY,
                                 AttributeGroup::TREE_TOPOLOGY, AttributeGroup::ALL}) {
        const auto baseline = AttributeComputation::computeAttributes(weighted, {group});
        const auto fromAltitudeSpan = AttributeComputation::computeAttributesFromAltitudeSpan(equivalentView, {group});
        requireComputedAttributesEquivalent(fromAltitudeSpan, baseline, tree, ATTRIBUTE_GROUPS.at(group),
                                            label + " attribute-pipeline group must match current path");
    }

    const auto singleBoundaryBaseline = AttributeComputation::computeSingleAttribute(weighted, AttributeGroup::BOUNDARY);
    const auto singleBoundaryFromAltitudeSpan = AttributeComputation::computeSingleAttribute(equivalentView, AttributeGroup::BOUNDARY);
    requireComputedAttributesEquivalent(singleBoundaryFromAltitudeSpan, singleBoundaryBaseline, tree, ATTRIBUTE_GROUPS.at(AttributeGroup::BOUNDARY),
                                        label + " attribute-pipeline single group overload must support BOUNDARY");

    requireImageNear(AttributeComputation::computeAttributeMapping(equivalentView, CONTOUR_PERIMETER),
                     AttributeComputation::computeAttributeMapping(weighted, CONTOUR_PERIMETER),
                     label + " altitude-span view mapping overload must support contour attributes");
}

void checkAttributePipelineMaxDistAndRejectsInvalidInputs(const WeightedMorphologicalTree<std::uint8_t>& weighted) {
    const auto& tree = weighted.topology();
    const std::vector<float> altitude = makeEquivalentAltitude<float>(weighted);

    const auto baselineMaxDist = AttributeComputation::computeSingleAttribute(weighted, MAX_DIST);
    const WeightedTreeView<float> view(tree, std::span<const float>(altitude));
    const auto maxDistFromAltitudeSpan = AttributeComputation::computeAttributesFromAltitudeSpan(view, {MAX_DIST});
    requireComputedAttributesNear(maxDistFromAltitudeSpan, baselineMaxDist, tree, {MAX_DIST}, "altitude-span attribute API must support MAX_DIST");

    const std::vector<std::int32_t> int32Altitude = makeEquivalentAltitude<std::int32_t>(weighted);
    const WeightedTreeView<std::int32_t> int32View(tree, std::span<const std::int32_t>(int32Altitude));
    const auto maxDistFromInt32AltitudeSpan = AttributeComputation::computeAttributesFromAltitudeSpan(int32View, {MAX_DIST});
    requireComputedAttributesNear(maxDistFromInt32AltitudeSpan, baselineMaxDist, tree, {MAX_DIST},
                                  "altitude-span attribute API must support MAX_DIST with int32 altitude");

    if constexpr (contract::validationsEnabled) {
        requireThrows<std::runtime_error>(
            [&]() {
                const std::vector<float> wrongSize(static_cast<std::size_t>(tree.getNumInternalNodeSlots() - 1), 0.0f);
                static_cast<void>(WeightedTreeView<float>(tree, std::span<const float>(wrongSize)));
            },
            "weighted view API must reject wrong altitude size");
    }
}

int main() {
    auto image = makeComponentTreeFixture();
    auto weighted = makeWeightedComponentTree(image, true);
    auto minWeighted = makeWeightedComponentTree(image, false);
    const NodeId sampleNodeId = findNonRootResidueSample(*weighted);

    checkGenericStaticAltitudeAccess<std::uint8_t>(*weighted, sampleNodeId);
    checkGenericStaticAltitudeAccess<std::int32_t>(*weighted, sampleNodeId);
    checkGenericStaticAltitudeAccess<float>(*weighted, sampleNodeId);
    checkGenericStaticAltitudeAccess<double>(*weighted, sampleNodeId);
    checkWeightedTreeViewContract<std::uint8_t>(*weighted, sampleNodeId, "uint8");
    checkWeightedTreeViewContract<std::int32_t>(*weighted, sampleNodeId, "int32");
    checkWeightedTreeViewContract<float>(*weighted, sampleNodeId, "float");
    checkWeightedTreeViewContract<double>(*weighted, sampleNodeId, "double");
    const auto ownedView = weighted->asView();
    require(&ownedView.topology() == &weighted->topology(), "WeightedMorphologicalTree<std::uint8_t>::asView must borrow the owned topology");
    requireEqual(static_cast<int>(ownedView.altitude().size()), weighted->topology().getNumInternalNodeSlots(),
                 "WeightedMorphologicalTree<std::uint8_t>::asView altitude span size");
    checkWeightedTreeViewIncrementalAttributePipeline(*weighted);
    checkFloatResidueAndMonotonicity(*weighted, sampleNodeId);
    checkUnsignedResidueDoesNotUnderflow(*weighted, sampleNodeId);
    checkFloatHigraExport(*weighted);
    checkInt32HigraExport(*weighted);
    checkGenericAltitudeDelta<std::int32_t>(*weighted, "max-tree int32");
    checkGenericAltitudeDelta<float>(*weighted, "max-tree float");
    checkGenericAltitudeDelta<double>(*weighted, "max-tree double");
    checkGenericAltitudeDelta<float>(*minWeighted, "min-tree float");
    checkGenericAltitudeDelta<double>(*minWeighted, "min-tree double");
    checkGenericVolumeKernel<float>(*weighted, "max-tree float");
    checkGenericVolumeKernel<double>(*weighted, "max-tree double");
    checkGenericGrayLevelStatsKernel<float>(*weighted, "max-tree float");
    checkGenericGrayLevelStatsKernel<double>(*weighted, "max-tree double");
    checkGenericGrayLevelStatsKernel<float>(*minWeighted, "min-tree float");
    checkGenericGrayLevelStatsKernel<double>(*minWeighted, "min-tree double");
    checkAttributePipelineApi<float>(*weighted, "max-tree float");
    checkAttributePipelineApi<double>(*weighted, "max-tree double");
    checkAttributePipelineApi<float>(*minWeighted, "min-tree float");
    checkAttributePipelineApi<double>(*minWeighted, "min-tree double");
    checkTopologyOnlyAttributePipelineApi<std::int32_t>(*weighted, "max-tree int32");
    checkTopologyOnlyAttributePipelineApi<float>(*weighted, "max-tree float");
    checkTopologyOnlyAttributePipelineApi<double>(*minWeighted, "min-tree double");
    checkAttributePipelineMaxDistAndRejectsInvalidInputs(*weighted);
    checkImageFactoriesDeduceAltitudeType(image);
    checkTypedOwnerExportedHigraAttributeProjection<std::int32_t>("int32 owner");
    checkTypedOwnerExportedHigraAttributeProjection<float>("float owner");
    checkTypedOwnerDeltaAttributeApi<std::int32_t>("int32 owner");
    checkTypedOwnerDeltaAttributeApi<float>("float owner");
    checkTypedHigraImportFactory(image);
    if constexpr (contract::validationsEnabled) {
        checkFiniteFloatAltitudeValidation(image);
    }

    return 0;
}
