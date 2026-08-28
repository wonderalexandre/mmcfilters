#include "support/TestSupport.hpp"

#include "mmcfilters/trees/detail/TreeAttributeSamplingNeighborhood.hpp"
#include "mmcfilters/attributes/computers/GrayLevelStatsComputer.hpp"
#include "mmcfilters/attributes/computers/VolumeComputer.hpp"
#include "mmcfilters/trees/MorphologicalTreeFactory.hpp"
#include "mmcfilters/trees/TreeAltitudeAlgorithms.hpp"
#include "mmcfilters/trees/ValuedMorphologicalTree.hpp"
#include "mmcfilters/trees/ValuedMorphologicalTreeView.hpp"

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

template class mmcfilters::ValuedMorphologicalTree<float>;
template class mmcfilters::ValuedMorphologicalTreeEditor<float>;

static_assert(std::is_same_v<typename ValuedMorphologicalTree<float>::AltitudeType, float>);
static_assert(std::is_same_v<decltype(std::declval<const ValuedMorphologicalTree<float>&>().nodeAltitudes()), const std::vector<float>&>);
static_assert(std::is_same_v<decltype(std::declval<const ValuedMorphologicalTree<float>&>().asView()), ValuedMorphologicalTreeView<float>>);
static_assert(std::is_constructible_v<ValuedMorphologicalTreeView<std::int32_t>, const MorphologicalTree&, const NodeAltitudeBuffer<std::int32_t>&>);
static_assert(!std::is_constructible_v<ValuedMorphologicalTreeView<std::int32_t>, const MorphologicalTree&, NodeAltitudeBuffer<std::int32_t>&&>);
static_assert(std::is_same_v<decltype(MorphologicalTreeFactory::createMaxTree(std::declval<ImageInt32Ptr>(), 1.5)), ValuedMorphologicalTree<std::int32_t>>);
static_assert(std::is_same_v<decltype(MorphologicalTreeFactory::createMinTree(std::declval<ImageFloatPtr>(), 1.5)), ValuedMorphologicalTree<float>>);
static_assert(std::is_same_v<decltype(MorphologicalTreeFactory::createFromHigraParent<double>(
                                 std::declval<std::span<const NodeId>>(), std::declval<std::span<const double>>(), 4, 4, MorphologicalTreeKind::MaxTree,
                                 std::declval<std::optional<RegularGridAdjacency2D>>())),
                             ValuedMorphologicalTree<double>>);
static_assert(std::is_same_v<decltype(std::declval<const ValuedMorphologicalTree<std::int32_t>&>().reconstructFromNodeAltitudes()), ImageInt32Ptr>);
static_assert(std::is_same_v<decltype(std::declval<const ValuedMorphologicalTree<float>&>().reconstructFromNodeAltitudes()), ImageFloatPtr>);
static_assert(
    std::is_same_v<decltype(AttributeComputation::projectNodeValuesToExportedHigra(
                       std::declval<const ValuedMorphologicalTree<float>&>(), std::declval<const AttributeNames&>(), std::declval<std::span<const float>>())),
                   std::vector<float>>);
static_assert(std::is_same_v<decltype(AttributeComputation::projectNodeValuesToExportedHigra<double>(std::declval<const ValuedMorphologicalTree<float>&>(),
                                                                                                     std::declval<const AttributeNames&>(),
                                                                                                     std::declval<std::span<const double>>())),
                             std::vector<double>>);
static_assert(std::is_same_v<decltype(AttributeComputation::computeAttributes(std::declval<const ValuedMorphologicalTree<std::uint8_t>&>(),
                                                                              std::declval<const std::vector<AttributeOrGroup>&>())),
                             ComputedAttributeData<float>>);
static_assert(std::is_same_v<decltype(AttributeComputation::computeAttributes<double>(std::declval<const ValuedMorphologicalTree<std::uint8_t>&>(),
                                                                                      std::declval<const std::vector<AttributeOrGroup>&>())),
                             ComputedAttributeData<double>>);
static_assert(std::is_same_v<decltype(AttributeComputation::computeTopologyAttributes<double>(std::declval<const MorphologicalTree&>(),
                                                                                              std::declval<const std::vector<AttributeOrGroup>&>())),
                             ComputedAttributeData<double>>);
static_assert(std::is_same_v<decltype(AttributeComputation::computeSampledNodeAttribute(std::declval<const ValuedMorphologicalTree<std::int32_t>&>(),
                                                                                            Area, std::declval<AltitudeDifference<std::int32_t>>(), 2)),
                             SampledNodeAttributeData<float>>);
static_assert(std::is_same_v<decltype(AttributeComputation::computeSampledNodeAttribute(std::declval<const ValuedMorphologicalTree<float>&>(), MeanGrayLevel,
                                                                                            std::declval<AltitudeDifference<float>>(), 2)),
                             SampledNodeAttributeData<float>>);
static_assert(std::is_same_v<decltype(AttributeComputation::computeSampledNodeAttribute<double>(
                                 std::declval<const ValuedMorphologicalTree<std::int32_t>&>(), Area, std::declval<AltitudeDifference<std::int32_t>>(), 2)),
                             SampledNodeAttributeData<double>>);
static_assert(
    std::is_same_v<decltype(AttributeComputation::computeAttributeMapping<double>(std::declval<const ValuedMorphologicalTree<std::uint8_t>&>(), MeanGrayLevel)),
                   ImagePtr<double>>);

template <class T> std::vector<T> makeGenericAltitude(const ValuedMorphologicalTree<std::uint8_t>& valuedTree) {
    std::vector<T> altitude(static_cast<std::size_t>(valuedTree.topology().numInternalNodeSlots()), T{});

    for (NodeId nodeId : valuedTree.topology().aliveNodeIds()) {
        const std::uint8_t base = valuedTree.nodeAltitude(nodeId);
        if constexpr (std::is_floating_point_v<T>) {
            altitude[static_cast<std::size_t>(nodeId)] = static_cast<T>(base) + static_cast<T>(0.25);
        } else {
            altitude[static_cast<std::size_t>(nodeId)] = static_cast<T>(base);
        }
    }

    return altitude;
}

template <class T> std::vector<T> makeEquivalentAltitude(const ValuedMorphologicalTree<std::uint8_t>& valuedTree) {
    std::vector<T> altitude(static_cast<std::size_t>(valuedTree.topology().numInternalNodeSlots()), T{});

    for (NodeId nodeId : valuedTree.topology().aliveNodeIds()) {
        altitude[static_cast<std::size_t>(nodeId)] = static_cast<T>(valuedTree.nodeAltitude(nodeId));
    }

    return altitude;
}

NodeId findNonRootResidueSample(const ValuedMorphologicalTree<std::uint8_t>& valuedTree) {
    const auto& tree = valuedTree.topology();
    for (NodeId nodeId : tree.aliveNodeIds()) {
        if (tree.isRoot(nodeId)) {
            continue;
        }
        const NodeId parentNodeId = tree.parent(nodeId);
        if (valuedTree.nodeAltitude(nodeId) != valuedTree.nodeAltitude(parentNodeId)) {
            return nodeId;
        }
    }
    throw std::runtime_error("fixture must contain a non-root node with non-zero residue");
}

void checkRootResidueUsesFixedZeroBaseline() {
    auto image = ImageUInt8::create(1, 1, std::uint8_t{37});
    const auto maxTree = MorphologicalTreeFactory::createMaxTree(image, 1.5);
    const auto minTree = MorphologicalTreeFactory::createMinTree(image, 1.5);

    for (const auto* valuedTree : {&maxTree, &minTree}) {
        const NodeId root = valuedTree->topology().root();
        requireEqual(valuedTree->nodeAltitude(root), std::uint8_t{37}, "single-node tree root altitude");
        requireEqual(valuedTree->nodeResidue(root), AltitudeDifference<std::uint8_t>{37},
                     "root residue must use the fixed zero reconstruction baseline");
    }
}

template <class T> void checkGenericStaticAltitudeAccess(const ValuedMorphologicalTree<std::uint8_t>& valuedTree, NodeId sampleNodeId) {
    const auto& tree = valuedTree.topology();
    const std::vector<T> altitude = makeGenericAltitude<T>(valuedTree);
    const std::span<const T> view(altitude);

    static_assert(std::is_same_v<decltype(TreeAltitudeAlgorithms::nodeAltitude(view, sampleNodeId)), T>);

    TreeAltitudeAlgorithms::validateNodeAltitudeBufferShape(tree, view);
    require(TreeAltitudeAlgorithms::nodeAltitude(view, sampleNodeId) == altitude[static_cast<std::size_t>(sampleNodeId)],
            "templated nodeAltitude must preserve the altitude value type");

    if constexpr (contract::validationsEnabled) {
        requireThrows<std::runtime_error>(
            [&]() {
                const std::vector<T> wrongSize(static_cast<std::size_t>(tree.numInternalNodeSlots() - 1), T{});
                TreeAltitudeAlgorithms::validateNodeAltitudeBufferShape(tree, std::span<const T>(wrongSize));
            },
            "templated validateNodeAltitudeBufferShape must reject wrong size");

        requireThrows<std::invalid_argument>([&]() { static_cast<void>(TreeAltitudeAlgorithms::nodeAltitude(view, InvalidNode)); },
                                             "templated nodeAltitude must reject invalid node ids");
    }
}

template <class T> void checkValuedMorphologicalTreeViewContract(const ValuedMorphologicalTree<std::uint8_t>& valuedTree, NodeId sampleNodeId, const std::string& label) {
    const auto& tree = valuedTree.topology();
    const std::vector<T> altitude = makeGenericAltitude<T>(valuedTree);
    const ValuedMorphologicalTreeView<T> view(tree, std::span<const T>(altitude));
    auto inferredSpanView = ValuedMorphologicalTreeView(tree, std::span<const T>(altitude));
    auto inferredBufferView = ValuedMorphologicalTreeView(tree, altitude);

    static_assert(std::is_same_v<typename ValuedMorphologicalTreeView<T>::AltitudeType, T>);
    static_assert(std::is_same_v<decltype(inferredSpanView), ValuedMorphologicalTreeView<T>>);
    static_assert(std::is_same_v<decltype(inferredBufferView), ValuedMorphologicalTreeView<T>>);

    require(&view.topology() == &tree, label + " view must reference the original topology");
    require(&inferredSpanView.topology() == &tree, label + " span view must reference the original topology");
    require(&inferredBufferView.topology() == &tree, label + " buffer view must reference the original topology");
    requireEqual(static_cast<int>(view.nodeAltitudes().size()), tree.numInternalNodeSlots(), label + " view altitude span size");
    requireEqual(static_cast<int>(inferredSpanView.nodeAltitudes().size()), tree.numInternalNodeSlots(), label + " inferred span view altitude span size");
    requireEqual(static_cast<int>(inferredBufferView.nodeAltitudes().size()), tree.numInternalNodeSlots(), label + " inferred buffer view altitude span size");
    require(view.nodeAltitude(sampleNodeId) == altitude[static_cast<std::size_t>(sampleNodeId)],
            label + " view nodeAltitude must preserve the altitude value type");
    require(inferredSpanView.nodeAltitude(sampleNodeId) == altitude[static_cast<std::size_t>(sampleNodeId)],
            label + " inferred span view nodeAltitude must preserve the altitude value type");
    require(inferredBufferView.nodeAltitude(sampleNodeId) == altitude[static_cast<std::size_t>(sampleNodeId)],
            label + " inferred buffer view nodeAltitude must preserve the altitude value type");

    const auto expectedResidue = TreeAltitudeAlgorithms::nodeResidue(tree, std::span<const T>(altitude), sampleNodeId);
    if constexpr (std::is_floating_point_v<T>) {
        requireNear(view.nodeResidue(sampleNodeId), expectedResidue, static_cast<T>(1.0e-6), label + " view residue must match static generic helper");
    } else {
        requireEqual(view.nodeResidue(sampleNodeId), expectedResidue, label + " view residue must match static generic helper");
    }

    const auto fromView = AttributeComputation::computeAttributesFromAltitudeSpan(view, {MeanGrayLevel, GrayLevelHeight});
    const auto fromSpan = AttributeComputation::computeAttributesFromAltitudeSpan(inferredSpanView, {MeanGrayLevel, GrayLevelHeight});
    for (NodeId nodeId : tree.aliveNodeIds()) {
        for (Attribute attribute : {MeanGrayLevel, GrayLevelHeight}) {
            requireNear(fromView.second[fromView.first.linearIndex(nodeId, attribute)], fromSpan.second[fromSpan.first.linearIndex(nodeId, attribute)], 1.0e-5f,
                        label + " view attribute pipeline must match inferred span view");
        }
    }

    if constexpr (contract::validationsEnabled) {
        requireThrows<std::runtime_error>(
            [&]() {
                const std::vector<T> wrongSize(static_cast<std::size_t>(tree.numInternalNodeSlots() - 1), T{});
                static_cast<void>(ValuedMorphologicalTreeView<T>(tree, std::span<const T>(wrongSize)));
            },
            label + " view must reject wrong altitude size");

        requireThrows<std::invalid_argument>([&]() { static_cast<void>(view.nodeAltitude(InvalidNode)); }, label + " view must reject invalid node ids");
    }
}

void requireComputedAttributesNear(const ComputedAttributeData<float>& actual, const ComputedAttributeData<float>& expected, const MorphologicalTree& tree,
                                   std::initializer_list<Attribute> attributes, const std::string& label) {
    for (NodeId nodeId : tree.aliveNodeIds()) {
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
    for (NodeId nodeId : tree.aliveNodeIds()) {
        for (Attribute attribute : attributes) {
            requireFloatEquivalent(actual.second[actual.first.linearIndex(nodeId, attribute)], expected.second[expected.first.linearIndex(nodeId, attribute)],
                                   label + " " + AttributeNames::toString(attribute) + " node " + std::to_string(nodeId));
        }
    }
}

void requireSampledNodeAttributeNear(const SampledNodeAttributeData<float>& actual, const SampledNodeAttributeData<float>& expected,
                                     const MorphologicalTree& tree, Attribute attribute, int samplingRadius, const std::string& label) {
    for (NodeId nodeId : tree.aliveNodeIds()) {
        for (int sampleOffset = -samplingRadius; sampleOffset <= samplingRadius; ++sampleOffset) {
            requireNear(actual.second[actual.first.linearIndex(nodeId, attribute, sampleOffset)],
                        expected.second[expected.first.linearIndex(nodeId, attribute, sampleOffset)], 1.0e-5f,
                        label + " node " + std::to_string(nodeId) + " sample offset " + std::to_string(sampleOffset));
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
void requireTypedOwnerMatchesCanonicalImageTree(const ValuedMorphologicalTree<T>& typed, const ValuedMorphologicalTree<std::uint8_t>& canonical,
                                                const ImageUInt8Ptr& expectedReconstruction, const std::string& label) {
    requireEqual(typed.topology().numInternalNodeSlots(), canonical.topology().numInternalNodeSlots(), label + " internal node slot count");
    requireEqual(typed.topology().numPixels(), canonical.topology().numPixels(), label + " proper part count");
    requireEqual(typed.topology().root(), canonical.topology().root(), label + " root");
    requireEqual(static_cast<int>(typed.nodeAltitudes().size()), typed.topology().numInternalNodeSlots(), label + " typed altitude buffer size");

    for (NodeId nodeId : typed.topology().aliveNodeIds()) {
        requireNear(static_cast<double>(typed.nodeAltitude(nodeId)), static_cast<double>(canonical.nodeAltitude(nodeId)), 1.0e-6, label + " node altitude");
    }

    std::vector<T> expectedTypedReconstruction;
    expectedTypedReconstruction.reserve(static_cast<std::size_t>(expectedReconstruction->getSize()));
    for (auto value : collectImageValues(expectedReconstruction)) {
        expectedTypedReconstruction.push_back(static_cast<T>(value));
    }
    requireVectorEqual(collectImageValues(typed.reconstructFromNodeAltitudes()), expectedTypedReconstruction, label + " typed reconstruction");
}

void checkTypedHigraImportFactory(ImageUInt8Ptr image) {
    auto canonicalMax = MorphologicalTreeFactory::createMaxTree(image, 1.5);

    const auto [higraParent, higraAltitude] = canonicalMax.exportHigraHierarchy();
    std::vector<float> floatHigraAltitude(higraAltitude.begin(), higraAltitude.end());
    auto importedFloat = MorphologicalTreeFactory::createFromHigraParent<float>(
        std::span<const NodeId>(higraParent), std::span<const float>(floatHigraAltitude), image->getNumRows(), image->getNumColumns(),
        MorphologicalTreeKind::MaxTree, RegularGridAdjacency2D(image->getNumRows(), image->getNumColumns(), 1.5));
    static_assert(std::is_same_v<decltype(importedFloat), ValuedMorphologicalTree<float>>);
    importedFloat.validateMonotoneNodeAltitudes();

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
    requireVectorEqual(collectImageValues(importedFloat.reconstructFromNodeAltitudes()), expectedFloatReconstruction, "typed Higra factory reconstruction");

    requireThrows<std::invalid_argument>(
        [&]() {
            const std::vector<double> wrongAltitude(higraParent.size() - 1, 0.0);
            static_cast<void>(MorphologicalTreeFactory::createFromHigraParent<double>(
                std::span<const NodeId>(higraParent), std::span<const double>(wrongAltitude), image->getNumRows(), image->getNumColumns(),
                MorphologicalTreeKind::MaxTree, RegularGridAdjacency2D(image->getNumRows(), image->getNumColumns(), 1.5)));
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
    const NodeId sampleNodeId = floatMax.topology().root();

    std::vector<float> invalidBuffer = floatMax.nodeAltitudes();
    invalidBuffer[static_cast<std::size_t>(sampleNodeId)] = std::numeric_limits<float>::quiet_NaN();
    requireThrows<std::invalid_argument>([&]() { floatMax.setNodeAltitudes(std::move(invalidBuffer)); },
                                         "ValuedMorphologicalTree<float>::setNodeAltitudes must reject NaN");

    requireThrows<std::invalid_argument>([&]() { floatMax.setNodeAltitude(sampleNodeId, std::numeric_limits<float>::infinity()); },
                                         "ValuedMorphologicalTree<float>::setNodeAltitude must reject +inf");

    {
        auto editor = floatMax.edit();
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(editor.createDetachedNode(-std::numeric_limits<float>::infinity())); },
                                             "ValuedMorphologicalTreeEditor<float>::createDetachedNode must reject -inf");
        editor.commit();
    }

    {
        auto editor = floatMax.edit();
        requireThrows<std::invalid_argument>([&]() { editor.setNodeAltitude(sampleNodeId, std::numeric_limits<float>::quiet_NaN()); },
                                             "ValuedMorphologicalTreeEditor<float>::setNodeAltitude must reject NaN");
        editor.commit();
    }

    const auto canonicalMax = MorphologicalTreeFactory::createMaxTree(image, 1.5);
    const auto [higraParent, higraAltitude] = canonicalMax.exportHigraHierarchy();
    std::vector<float> floatHigraAltitude(higraAltitude.begin(), higraAltitude.end());
    floatHigraAltitude.back() = std::numeric_limits<float>::infinity();
    requireThrows<std::invalid_argument>(
        [&]() {
            static_cast<void>(MorphologicalTreeFactory::createFromHigraParent<float>(
                std::span<const NodeId>(higraParent), std::span<const float>(floatHigraAltitude), image->getNumRows(), image->getNumColumns(),
                MorphologicalTreeKind::MaxTree, RegularGridAdjacency2D(image->getNumRows(), image->getNumColumns(), 1.5)));
        },
        "typed Higra import must reject non-finite float altitude");
}

void checkImageFactoriesDeduceAltitudeType(ImageUInt8Ptr image) {
    const auto canonicalMax = MorphologicalTreeFactory::createMaxTree(image, 1.5);
    const auto canonicalMin = MorphologicalTreeFactory::createMinTree(image, 1.5);

    auto int32Max = MorphologicalTreeFactory::createMaxTree(makeTypedComponentTreeFixture<std::int32_t>(), 1.5);
    static_assert(std::is_same_v<decltype(int32Max), ValuedMorphologicalTree<std::int32_t>>);
    int32Max.validateMonotoneNodeAltitudes();
    requireTypedOwnerMatchesCanonicalImageTree(int32Max, canonicalMax, image, "Image<int32_t> max-tree factory");
    auto int32Min = MorphologicalTreeFactory::createMinTree(makeTypedComponentTreeFixture<std::int32_t>(), 1.5);
    static_assert(std::is_same_v<decltype(int32Min), ValuedMorphologicalTree<std::int32_t>>);
    int32Min.validateMonotoneNodeAltitudes();
    requireTypedOwnerMatchesCanonicalImageTree(int32Min, canonicalMin, image, "Image<int32_t> min-tree factory");
    {
        const std::vector<AttributeOrGroup> requests{AttributeGroup::Shape, AttributeGroup::TreeTopology};
        const auto typedAttrs = AttributeComputation::computeAttributes(int32Max, requests);
        const auto canonicalAttrs = AttributeComputation::computeAttributes(canonicalMax, requests);
        std::vector<Attribute> expectedAttrs = ATTRIBUTE_GROUPS.at(AttributeGroup::Shape);
        const auto topologyAttrs = ATTRIBUTE_GROUPS.at(AttributeGroup::TreeTopology);
        expectedAttrs.insert(expectedAttrs.end(), topologyAttrs.begin(), topologyAttrs.end());
        requireComputedAttributesEquivalent(typedAttrs, canonicalAttrs, int32Max.topology(), expectedAttrs,
                                            "Image<int32_t> owner must compute semantic shape and topology groups");
        requireImageNear(AttributeComputation::computeAttributeMapping(int32Max, ContourPerimeter),
                         AttributeComputation::computeAttributeMapping(canonicalMax, ContourPerimeter),
                         "Image<int32_t> owner mapping overload must support contour attributes");
    }

    auto floatMax = MorphologicalTreeFactory::createMaxTree(makeTypedComponentTreeFixture<float>(), 1.5);
    static_assert(std::is_same_v<decltype(floatMax), ValuedMorphologicalTree<float>>);
    floatMax.validateMonotoneNodeAltitudes();
    requireTypedOwnerMatchesCanonicalImageTree(floatMax, canonicalMax, image, "Image<float> max-tree factory");

    auto floatMin = MorphologicalTreeFactory::createMinTree(makeTypedComponentTreeFixture<float>(), 1.5);
    static_assert(std::is_same_v<decltype(floatMin), ValuedMorphologicalTree<float>>);
    floatMin.validateMonotoneNodeAltitudes();
    requireTypedOwnerMatchesCanonicalImageTree(floatMin, canonicalMin, image, "Image<float> min-tree factory");
    {
        const std::vector<AttributeOrGroup> requests{MeanGrayLevel, AttributeGroup::Moments};
        const auto typedAttrs = AttributeComputation::computeAttributes(floatMin, requests);
        const auto canonicalAttrs = AttributeComputation::computeAttributes(canonicalMin, requests);
        std::vector<Attribute> expectedAttrs{MeanGrayLevel};
        const auto momentAttrs = ATTRIBUTE_GROUPS.at(AttributeGroup::Moments);
        expectedAttrs.insert(expectedAttrs.end(), momentAttrs.begin(), momentAttrs.end());
        requireComputedAttributesEquivalent(typedAttrs, canonicalAttrs, floatMin.topology(), expectedAttrs,
                                            "Image<float> owner must compute mixed generic attributes");

        const auto singleMeanGrayLevel = AttributeComputation::computeSingleAttribute(floatMin, MeanGrayLevel);
        const auto canonicalSingleMeanGrayLevel = AttributeComputation::computeSingleAttribute(canonicalMin, MeanGrayLevel);
        requireComputedAttributesEquivalent(singleMeanGrayLevel, canonicalSingleMeanGrayLevel, floatMin.topology(), {MeanGrayLevel},
                                            "Image<float> owner single-attribute overload must compute MEAN_GRAY_LEVEL");

        requireImageNear(AttributeComputation::computeAttributeMapping(floatMin, MeanGrayLevel),
                         AttributeComputation::computeAttributeMapping(canonicalMin, MeanGrayLevel),
                         "Image<float> owner mapping overload must support MEAN_GRAY_LEVEL");
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
    const std::vector<AttributeOrGroup> requests{Area,           Volume,      RelativeVolume, MeanGrayLevel,
                                                 GrayLevelVariance, GrayLevelHeight, MaxDistExact,        BoxColumnMin, BoxRowMin};
    const auto computed = AttributeComputation::computeAttributes(typed, requests);
    const auto projected = AttributeComputation::projectNodeValuesToExportedHigra(typed, computed.first, computed.second);

    const auto expectedSize = static_cast<std::size_t>(typed.topology().numPixels() + typed.topology().numNodes()) *
                              static_cast<std::size_t>(computed.first.NUM_ATTRIBUTES);
    requireEqual(projected.size(), expectedSize, label + " exported-Higra typed projection size");

    const NodeId sampleProperPart = 10;
    const NodeId sampleSmallestNode = typed.topology().smallestNode(sampleProperPart);
    const auto [sampleRow, sampleColumn] = ImageUtils::to2D(sampleProperPart, typed.topology().numColumns());
    const float sampleAltitude = static_cast<float>(typed.nodeAltitude(sampleSmallestNode));
    auto projectedAt = [&](Attribute attribute) { return projected[computed.first.linearIndex(sampleProperPart, attribute)]; };

    requireNear(projectedAt(Area), 1.0f, 1.0e-5f, label + " unit AREA");
    requireNear(projectedAt(Volume), sampleAltitude, 1.0e-5f, label + " unit VOLUME must preserve typed altitude");
    requireNear(projectedAt(RelativeVolume), 1.0f, 1.0e-5f, label + " unit RELATIVE_VOLUME");
    requireNear(projectedAt(MeanGrayLevel), sampleAltitude, 1.0e-5f, label + " unit MeanGrayLevel must preserve typed altitude");
    requireNear(projectedAt(GrayLevelVariance), 0.0f, 1.0e-5f, label + " unit GrayLevelVariance");
    requireNear(projectedAt(GrayLevelHeight), 0.0f, 1.0e-5f, label + " unit GrayLevelHeight");
    requireNear(projectedAt(MaxDistExact), 0.0f, 1.0e-5f, label + " unit MAX_DIST_EXACT");
    requireNear(projectedAt(BoxColumnMin), static_cast<float>(sampleColumn), 1.0e-5f, label + " unit BOX_COLUMN_MIN");
    requireNear(projectedAt(BoxRowMin), static_cast<float>(sampleRow), 1.0e-5f, label + " unit BOX_ROW_MIN");
}

void checkValuedMorphologicalTreeViewIncrementalAttributePipeline(const ValuedMorphologicalTree<std::uint8_t>& valuedTree) {
    const auto& tree = valuedTree.topology();
    const auto ownerView = valuedTree.asView();

    const auto singleBaseline = AttributeComputation::computeSingleAttribute(valuedTree, MaxDistExact);
    const auto singleFromView = AttributeComputation::computeSingleAttribute(ownerView, MaxDistExact);
    requireComputedAttributesNear(singleFromView, singleBaseline, tree, {MaxDistExact}, "owner view single attribute");

    const auto multiBaseline = AttributeComputation::computeAttributes(valuedTree, {Area, MaxDistExact});
    const auto multiFromView = AttributeComputation::computeAttributes(ownerView, {Area, MaxDistExact});
    requireComputedAttributesNear(multiFromView, multiBaseline, tree, {Area, MaxDistExact}, "owner view multi attribute");

    requireImageNear(AttributeComputation::computeAttributeMapping(ownerView, MeanGrayLevel),
                     AttributeComputation::computeAttributeMapping(valuedTree, MeanGrayLevel),
                     "owner view attribute mapping");

    const std::vector<std::uint8_t> externalAltitude = makeEquivalentAltitude<std::uint8_t>(valuedTree);
    const ValuedMorphologicalTreeView<std::uint8_t> externalView(tree, std::span<const std::uint8_t>(externalAltitude));

    const auto externalRegular = AttributeComputation::computeAttributes(externalView, {MeanGrayLevel, GrayLevelHeight});
    const auto externalFromAltitudeSpan = AttributeComputation::computeAttributesFromAltitudeSpan(externalView, {MeanGrayLevel, GrayLevelHeight});
    requireComputedAttributesNear(externalRegular, externalFromAltitudeSpan, tree, {MeanGrayLevel, GrayLevelHeight},
                                  "external view regular overload must use attribute-pipeline subset");

    const auto externalSingle = AttributeComputation::computeSingleAttribute(externalView, MeanGrayLevel);
    requireComputedAttributesNear(externalSingle, externalFromAltitudeSpan, tree, {MeanGrayLevel},
                                  "external view single overload must use attribute-pipeline subset");

    const auto externalMaxDist = AttributeComputation::computeSingleAttribute(externalView, MaxDistExact);
    requireComputedAttributesNear(externalMaxDist, singleBaseline, tree, {MaxDistExact}, "external view attribute-pipeline path must support MAX_DIST_EXACT");

    const auto externalMixed = AttributeComputation::computeAttributes(externalView, {MeanGrayLevel, MaxDistExact});
    const auto baselineMixed = AttributeComputation::computeAttributes(valuedTree, {MeanGrayLevel, MaxDistExact});
    requireComputedAttributesNear(externalMixed, baselineMixed, tree, {MeanGrayLevel, MaxDistExact},
                                  "external view attribute-pipeline path must support mixed MAX_DIST_EXACT requests");

    requireImageNear(AttributeComputation::computeAttributeMapping(externalView, MeanGrayLevel),
                     AttributeComputation::computeAttributeMapping(valuedTree, MeanGrayLevel),
                     "external view attribute-pipeline attribute mapping");
}

void checkFloatResidueAndMonotonicity(const ValuedMorphologicalTree<std::uint8_t>& valuedTree, NodeId sampleNodeId) {
    const auto& tree = valuedTree.topology();
    std::vector<float> altitude = makeGenericAltitude<float>(valuedTree);
    const std::span<const float> view(altitude);

    const NodeId parentNodeId = tree.parent(sampleNodeId);
    const float expectedResidue = altitude[static_cast<std::size_t>(sampleNodeId)] - altitude[static_cast<std::size_t>(parentNodeId)];
    requireNear(TreeAltitudeAlgorithms::nodeResidue(tree, view, sampleNodeId), expectedResidue, 1.0e-6f,
                "templated nodeResidue must preserve float arithmetic");

    TreeAltitudeAlgorithms::validateMonotoneNodeAltitudes(tree, view);

    std::vector<float> nonMonotone = altitude;
    nonMonotone[static_cast<std::size_t>(parentNodeId)] = nonMonotone[static_cast<std::size_t>(sampleNodeId)] + 1.0f;
    requireThrows<std::runtime_error>([&]() { TreeAltitudeAlgorithms::validateMonotoneNodeAltitudes(tree, std::span<const float>(nonMonotone)); },
                                      "templated validateMonotoneNodeAltitudes must reject non-monotone float altitude");
}

void checkUnsignedResidueDoesNotUnderflow(const ValuedMorphologicalTree<std::uint8_t>& valuedTree, NodeId sampleNodeId) {
    const auto& tree = valuedTree.topology();
    const NodeId parentNodeId = tree.parent(sampleNodeId);
    require(parentNodeId != InvalidNode && parentNodeId != sampleNodeId, "unsigned residue sample must be a non-root node");

    std::vector<std::uint32_t> altitude(static_cast<std::size_t>(tree.numInternalNodeSlots()), 0u);
    altitude[static_cast<std::size_t>(parentNodeId)] = 10u;
    altitude[static_cast<std::size_t>(sampleNodeId)] = 3u;

    const std::span<const std::uint32_t> view(altitude);
    const AltitudeDifference<std::uint32_t> expectedResidue = -7;
    requireEqual(TreeAltitudeAlgorithms::nodeResidue(tree, view, sampleNodeId), expectedResidue,
                 "templated nodeResidue must not underflow unsigned altitude subtraction");

    const ValuedMorphologicalTreeView<std::uint32_t> valuedTreeView(tree, view);
    requireEqual(valuedTreeView.nodeResidue(sampleNodeId), expectedResidue,
                 "ValuedMorphologicalTreeView::nodeResidue must not underflow unsigned altitude subtraction");
}

void checkFloatHigraExport(const ValuedMorphologicalTree<std::uint8_t>& valuedTree) {
    const auto& tree = valuedTree.topology();
    const std::vector<float> altitude = makeGenericAltitude<float>(valuedTree);
    const auto [parent, exportedAltitude] = TreeAltitudeAlgorithms::exportHigraHierarchy(tree, std::span<const float>(altitude));

    static_assert(std::is_same_v<std::remove_cvref_t<decltype(exportedAltitude)>, std::vector<float>>);

    requireEqual(parent.size(), exportedAltitude.size(), "templated Higra export parent/altitude size");
    requireEqual(static_cast<int>(parent.size()), tree.numPixels() + tree.numNodes(), "templated Higra export compact domain size");

    bool foundFractionalAltitude = false;
    for (float value : exportedAltitude) {
        if (std::abs((value - std::floor(value)) - 0.25f) < 1.0e-6f) {
            foundFractionalAltitude = true;
            break;
        }
    }
    require(foundFractionalAltitude, "templated Higra export must preserve fractional float altitudes");

    const int numLeaves = tree.numPixels();
    for (NodeId leafId = 0; leafId < numLeaves; ++leafId) {
        const NodeId smallestNodeHigraId = parent[static_cast<std::size_t>(leafId)];
        require(smallestNodeHigraId >= numLeaves && smallestNodeHigraId < static_cast<NodeId>(parent.size()),
                "templated Higra export leaf parent must be an internal Higra node");
        requireNear(exportedAltitude[static_cast<std::size_t>(leafId)], exportedAltitude[static_cast<std::size_t>(smallestNodeHigraId)], 1.0e-6f,
                    "templated Higra export leaf altitude must match smallest-node altitude");
    }
}

void checkInt32HigraExport(const ValuedMorphologicalTree<std::uint8_t>& valuedTree) {
    const auto& tree = valuedTree.topology();
    const std::vector<std::int32_t> altitude = makeGenericAltitude<std::int32_t>(valuedTree);
    const auto [parent, exportedAltitude] = TreeAltitudeAlgorithms::exportHigraHierarchy(tree, std::span<const std::int32_t>(altitude));

    static_assert(std::is_same_v<std::remove_cvref_t<decltype(exportedAltitude)>, std::vector<std::int32_t>>);

    requireEqual(parent.size(), exportedAltitude.size(), "int32 templated Higra export parent/altitude size");
    requireEqual(static_cast<int>(parent.size()), tree.numPixels() + tree.numNodes(), "int32 templated Higra export compact domain size");

    const int numLeaves = tree.numPixels();
    for (NodeId leafId = 0; leafId < numLeaves; ++leafId) {
        const NodeId smallestNodeHigraId = parent[static_cast<std::size_t>(leafId)];
        require(smallestNodeHigraId >= numLeaves && smallestNodeHigraId < static_cast<NodeId>(parent.size()),
                "int32 templated Higra export leaf parent must be an internal Higra node");
        requireEqual(exportedAltitude[static_cast<std::size_t>(leafId)], exportedAltitude[static_cast<std::size_t>(smallestNodeHigraId)],
                     "int32 templated Higra export leaf altitude must match smallest-node altitude");
    }
}

template <class T> void checkGenericAltitudeSampling(const ValuedMorphologicalTree<std::uint8_t>& valuedTree, const std::string& label) {
    const auto& tree = valuedTree.topology();
    const std::vector<T> altitude = makeGenericAltitude<T>(valuedTree);
    const std::span<const T> view(altitude);
    const auto baseline = detail::computeNodeAttributeSamplingNeighborhood(
        tree, std::span<const std::uint8_t>(valuedTree.nodeAltitudes()), static_cast<AltitudeDifference<std::uint8_t>>(2),
        NodeAttributeSamplingPolicy::LargestSupportDescendant);
    const auto fromAltitudeSpan = detail::computeNodeAttributeSamplingNeighborhood(
        tree, view, static_cast<AltitudeDifference<T>>(2), NodeAttributeSamplingPolicy::LargestSupportDescendant);

    requireVectorEqual(fromAltitudeSpan.ancestors, baseline.ancestors, label + " generic sampled ancestors");
    requireVectorEqual(fromAltitudeSpan.representativeDescendants, baseline.representativeDescendants,
                       label + " generic representative descendants");

    NodeId sampleNodeId = InvalidNode;
    for (NodeId nodeId : tree.aliveNodeIds()) {
        if (baseline.ancestors[static_cast<std::size_t>(nodeId)] != InvalidNode) {
            sampleNodeId = nodeId;
            break;
        }
    }
    require(sampleNodeId != InvalidNode, label + " generic sampling fixture needs one ancestor sample");
    requireEqual(detail::findAncestorByAltitudeDistance(tree, view, sampleNodeId, static_cast<AltitudeDifference<T>>(2)),
                 baseline.ancestors[static_cast<std::size_t>(sampleNodeId)], label + " generic findAncestorByAltitudeDistance");

    if constexpr (contract::validationsEnabled) {
        requireThrows<std::runtime_error>(
            [&]() {
                const std::vector<T> wrongSize(static_cast<std::size_t>(tree.numInternalNodeSlots() - 1), T{});
                static_cast<void>(detail::computeNodeAttributeSamplingNeighborhood(
                    tree, std::span<const T>(wrongSize), static_cast<AltitudeDifference<T>>(2),
                    NodeAttributeSamplingPolicy::LargestSupportDescendant));
            },
            label + " generic sampling must reject wrong altitude size");

        requireThrows<std::invalid_argument>(
            [&]() { static_cast<void>(detail::findAncestorByAltitudeDistance(tree, view, InvalidNode, static_cast<AltitudeDifference<T>>(2))); },
            label + " generic sampling must reject invalid node ids");
    }
}

template <class T> void checkTypedOwnerSampledAttributeApi(const std::string& label) {
    auto canonicalImage = makeTypedComponentTreeFixture<std::uint8_t>();
    auto canonical = MorphologicalTreeFactory::createMaxTree(canonicalImage, 1.5);

    auto typedImage = makeTypedComponentTreeFixture<T>();
    AltitudeDifference<T> altitudeStep{};
    if constexpr (std::is_floating_point_v<T>) {
        for (int i = 0; i < typedImage->getSize(); ++i) {
            (*typedImage)[i] = static_cast<T>((*typedImage)[i] * static_cast<T>(0.25));
        }
        altitudeStep = static_cast<AltitudeDifference<T>>(0.25);
    } else {
        for (int i = 0; i < typedImage->getSize(); ++i) {
            (*typedImage)[i] = static_cast<T>((*typedImage)[i] + static_cast<T>(1000));
        }
        altitudeStep = static_cast<AltitudeDifference<T>>(1);
    }

    auto typed = MorphologicalTreeFactory::createMaxTree(typedImage, 1.5);
    auto typedSamples = AttributeComputation::computeSampledNodeAttribute(typed, Area, altitudeStep, 2);
    auto canonicalSamples = AttributeComputation::computeSampledNodeAttribute(canonical, Area, AltitudeDifference<std::uint8_t>{1}, 2);
    requireSampledNodeAttributeNear(typedSamples, canonicalSamples, typed.topology(), Area, 2,
                                    label + " typed sampled AREA must match equivalent canonical layout");

    const auto typedMeanGrayLevel = AttributeComputation::computeSingleAttribute(typed, MeanGrayLevel);
    auto typedMeanGrayLevelSamples = AttributeComputation::computeSampledNodeAttribute(
        typed, MeanGrayLevel, altitudeStep, 1, NodeAttributeSamplingPolicy::LargestSupportDescendant, MissingNodeAttributeSamplePolicy::NotANumber);
    bool foundAncestorSample = false;
    const auto samplingNeighborhood = detail::computeNodeAttributeSamplingNeighborhood(
        typed.topology(), typed.nodeAltitudeSpan(), altitudeStep, NodeAttributeSamplingPolicy::LargestSupportDescendant);
    const auto& ancestors = samplingNeighborhood.ancestors;
    for (NodeId nodeId : typed.topology().aliveNodeIds()) {
        requireNear(typedMeanGrayLevelSamples.second[typedMeanGrayLevelSamples.first.linearIndex(nodeId, MeanGrayLevel, 0)],
                    typedMeanGrayLevel.second[typedMeanGrayLevel.first.linearIndex(nodeId, MeanGrayLevel)], 1.0e-5f,
                    label + " typed MEAN_GRAY_LEVEL center sample must preserve the node attribute");

        const NodeId ancestor = ancestors[static_cast<std::size_t>(nodeId)];
        if (ancestor != InvalidNode && ancestor != nodeId) {
            requireNear(typedMeanGrayLevelSamples.second[typedMeanGrayLevelSamples.first.linearIndex(nodeId, MeanGrayLevel, -1)],
                        typedMeanGrayLevel.second[typedMeanGrayLevel.first.linearIndex(ancestor, MeanGrayLevel)], 1.0e-5f,
                        label + " typed MEAN_GRAY_LEVEL ancestor sample must use typed altitude distance");
            foundAncestorSample = true;
            break;
        }
    }
    require(foundAncestorSample, label + " typed sampling fixture must contain an ancestor sample");

    if constexpr (contract::validationsEnabled) {
        requireThrows<std::invalid_argument>(
            [&]() { static_cast<void>(AttributeComputation::computeSampledNodeAttribute(typed, Area, altitudeStep, -1)); },
            label + " typed sampling must reject negative radius");

        requireThrows<std::invalid_argument>(
            [&]() { static_cast<void>(AttributeComputation::computeSampledNodeAttribute(typed, Area, AltitudeDifference<T>{}, 1)); },
            label + " typed sampling must reject zero altitude step");

        if constexpr (std::is_floating_point_v<T>) {
            requireThrows<std::invalid_argument>(
                [&]() {
                    static_cast<void>(AttributeComputation::computeSampledNodeAttribute(
                        typed, Area, std::numeric_limits<AltitudeDifference<T>>::quiet_NaN(), 1));
                },
                label + " typed sampling must reject a non-finite altitude step");
        }
    }
}

template <class T> void checkGenericVolumeKernel(const ValuedMorphologicalTree<std::uint8_t>& valuedTree, const std::string& label) {
    const auto& tree = valuedTree.topology();
    const auto areaComputed = AttributeComputation::computeSingleAttribute(valuedTree, Area);
    const auto baseline = AttributeComputation::computeAttributes(valuedTree, {Volume, RelativeVolume});

    const AttributeNames volumeNames = AttributeNames::fromList({Volume, RelativeVolume});
    const std::vector<Attribute> requested{Volume, RelativeVolume};
    const std::array<DependencySourceT<float>, 1> dependencies{{DependencySourceT<float>{&areaComputed.first, areaComputed.second.data()}}};
    std::vector<float> genericBuffer(static_cast<std::size_t>(tree.numInternalNodeSlots()) * static_cast<std::size_t>(volumeNames.NUM_ATTRIBUTES), 0.0f);

    const std::vector<T> equivalentAltitude = makeEquivalentAltitude<T>(valuedTree);
    const auto equivalentContext = AltitudeAttributeComputeContext<float, T>{tree, std::span<const T>(equivalentAltitude), std::span<float>(genericBuffer),
                                                                             volumeNames, requested,
                                                                             std::span<const DependencySourceT<float>>(dependencies)};
    const auto volumeRequest = mmcfilters::attributes::computers::detail::VolumeRequest::from(requested);
    mmcfilters::attributes::computers::detail::kernel::computeVolume(equivalentContext, volumeRequest, &dependencies[0]);

    for (NodeId nodeId : tree.aliveNodeIds()) {
        requireNear(genericBuffer[volumeNames.linearIndex(nodeId, Volume)], baseline.second[baseline.first.linearIndex(nodeId, Volume)], 1.0e-5f,
                    label + " generic VOLUME must match current path");
        requireNear(genericBuffer[volumeNames.linearIndex(nodeId, RelativeVolume)], baseline.second[baseline.first.linearIndex(nodeId, RelativeVolume)],
                    1.0e-5f, label + " generic RELATIVE_VOLUME must match current path");
    }

    if constexpr (std::is_floating_point_v<T>) {
        std::vector<float> fractionalBuffer(genericBuffer.size(), 0.0f);
        const std::vector<T> fractionalAltitude = makeGenericAltitude<T>(valuedTree);
        const auto fractionalContext = AltitudeAttributeComputeContext<float, T>{
            tree, std::span<const T>(fractionalAltitude), std::span<float>(fractionalBuffer), volumeNames, requested,
            std::span<const DependencySourceT<float>>(dependencies)};
        mmcfilters::attributes::computers::detail::kernel::computeVolume(fractionalContext, volumeRequest, &dependencies[0]);

        for (NodeId nodeId : tree.aliveNodeIds()) {
            const float area = areaComputed.second[areaComputed.first.linearIndex(nodeId, Area)];
            requireNear(fractionalBuffer[volumeNames.linearIndex(nodeId, Volume)], baseline.second[baseline.first.linearIndex(nodeId, Volume)] + (area * 0.25f),
                        1.0e-5f, label + " fractional VOLUME must preserve float altitude contribution");
            requireNear(fractionalBuffer[volumeNames.linearIndex(nodeId, RelativeVolume)],
                        baseline.second[baseline.first.linearIndex(nodeId, RelativeVolume)], 1.0e-5f,
                        label + " uniform fractional offset must preserve RELATIVE_VOLUME");
        }
    }
}

template <class T> void checkGenericGrayLevelStatsKernel(const ValuedMorphologicalTree<std::uint8_t>& valuedTree, const std::string& label) {
    const auto& tree = valuedTree.topology();
    const auto areaComputed = AttributeComputation::computeSingleAttribute(valuedTree, Area);
    const auto volumeComputed = AttributeComputation::computeAttributes(valuedTree, {Volume, RelativeVolume});
    const auto baseline = AttributeComputation::computeAttributes(valuedTree, {MeanGrayLevel, GrayLevelVariance, GrayLevelHeight});

    const AttributeNames grayNames = AttributeNames::fromList({MeanGrayLevel, GrayLevelVariance, GrayLevelHeight});
    const std::vector<Attribute> requested{MeanGrayLevel, GrayLevelVariance, GrayLevelHeight};
    const std::array<DependencySourceT<float>, 2> dependencies{{DependencySourceT<float>{&volumeComputed.first, volumeComputed.second.data()},
                                                                DependencySourceT<float>{&areaComputed.first, areaComputed.second.data()}}};
    std::vector<float> genericBuffer(static_cast<std::size_t>(tree.numInternalNodeSlots()) * static_cast<std::size_t>(grayNames.NUM_ATTRIBUTES), 0.0f);

    const std::vector<T> equivalentAltitude = makeEquivalentAltitude<T>(valuedTree);
    const auto equivalentContext = AltitudeAttributeComputeContext<float, T>{tree, std::span<const T>(equivalentAltitude), std::span<float>(genericBuffer),
                                                                             grayNames, requested,
                                                                             std::span<const DependencySourceT<float>>(dependencies)};
    const auto grayRequest = mmcfilters::attributes::computers::detail::GrayLevelStatsRequest::from(requested);
    mmcfilters::attributes::computers::detail::kernel::computeGrayLevelStats(equivalentContext, grayRequest, &dependencies[0], &dependencies[1]);

    for (NodeId nodeId : tree.aliveNodeIds()) {
        for (Attribute attribute : requested) {
            requireNear(genericBuffer[grayNames.linearIndex(nodeId, attribute)], baseline.second[baseline.first.linearIndex(nodeId, attribute)], 1.0e-5f,
                        label + " generic gray-level attribute must match current path");
        }
    }

    if constexpr (std::is_floating_point_v<T>) {
        const AttributeNames volumeNames = AttributeNames::fromList({Volume, RelativeVolume});
        const std::vector<Attribute> volumeRequested{Volume, RelativeVolume};
        const std::array<DependencySourceT<float>, 1> volumeDependencies{{DependencySourceT<float>{&areaComputed.first, areaComputed.second.data()}}};
        std::vector<float> fractionalVolumeBuffer(
            static_cast<std::size_t>(tree.numInternalNodeSlots()) * static_cast<std::size_t>(volumeNames.NUM_ATTRIBUTES), 0.0f);
        const std::vector<T> fractionalAltitude = makeGenericAltitude<T>(valuedTree);
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

        for (NodeId nodeId : tree.aliveNodeIds()) {
            requireNear(fractionalGrayBuffer[grayNames.linearIndex(nodeId, MeanGrayLevel)],
                        baseline.second[baseline.first.linearIndex(nodeId, MeanGrayLevel)] + 0.25f, 1.0e-5f,
                        label + " fractional MeanGrayLevel must preserve uniform altitude offset");
            requireNear(fractionalGrayBuffer[grayNames.linearIndex(nodeId, GrayLevelVariance)],
                        baseline.second[baseline.first.linearIndex(nodeId, GrayLevelVariance)], 1.0e-4f,
                        label + " uniform fractional offset must preserve GrayLevelVariance");
            requireNear(fractionalGrayBuffer[grayNames.linearIndex(nodeId, GrayLevelHeight)], baseline.second[baseline.first.linearIndex(nodeId, GrayLevelHeight)],
                        1.0e-5f, label + " uniform fractional offset must preserve GrayLevelHeight");
        }
    }
}

template <class T> void checkAttributePipelineApi(const ValuedMorphologicalTree<std::uint8_t>& valuedTree, const std::string& label) {
    const auto& tree = valuedTree.topology();
    const std::vector<AttributeOrGroup> requests{Area, Volume, RelativeVolume, MeanGrayLevel, GrayLevelVariance, GrayLevelHeight, MaxDistExact};
    const std::vector<Attribute> expectedAttributes{Area, Volume, RelativeVolume, GrayLevelHeight, MeanGrayLevel, GrayLevelVariance, MaxDistExact};
    const auto baseline = AttributeComputation::computeAttributes(valuedTree, requests);
    const std::vector<T> equivalentAltitude = makeEquivalentAltitude<T>(valuedTree);
    const ValuedMorphologicalTreeView<T> equivalentView(tree, std::span<const T>(equivalentAltitude));
    const auto fromAltitudeSpan = AttributeComputation::computeAttributesFromAltitudeSpan(equivalentView, requests);
    const auto viewOverload = AttributeComputation::computeAttributes(equivalentView, requests);

    for (NodeId nodeId : tree.aliveNodeIds()) {
        for (Attribute attribute : expectedAttributes) {
            requireNear(fromAltitudeSpan.second[fromAltitudeSpan.first.linearIndex(nodeId, attribute)],
                        baseline.second[baseline.first.linearIndex(nodeId, attribute)], 1.0e-5f,
                        label + " attribute pipeline must match current path for " + AttributeNames::toString(attribute));
            requireNear(viewOverload.second[viewOverload.first.linearIndex(nodeId, attribute)],
                        fromAltitudeSpan.second[fromAltitudeSpan.first.linearIndex(nodeId, attribute)], 1.0e-5f,
                        label + " altitude-span view overload must match equivalent view API for " + AttributeNames::toString(attribute));
        }
    }

    const auto singleViewOverload = AttributeComputation::computeSingleAttribute(equivalentView, MeanGrayLevel);
    for (NodeId nodeId : tree.aliveNodeIds()) {
        requireNear(singleViewOverload.second[singleViewOverload.first.linearIndex(nodeId, MeanGrayLevel)],
                    baseline.second[baseline.first.linearIndex(nodeId, MeanGrayLevel)], 1.0e-5f,
                    label + " altitude-span single view overload must match current MEAN_GRAY_LEVEL");
    }

    requireImageNear(AttributeComputation::computeAttributeMapping(equivalentView, MeanGrayLevel),
                     AttributeComputation::computeAttributeMapping(valuedTree, MeanGrayLevel),
                     label + " altitude-span view mapping overload must match current MEAN_GRAY_LEVEL mapping");

    const auto grayLevelFromAltitudeSpan = AttributeComputation::computeAttributesFromAltitudeSpan(equivalentView, {AttributeGroup::GrayLevel});
    for (NodeId nodeId : tree.aliveNodeIds()) {
        for (Attribute attribute : {Volume, RelativeVolume, GrayLevelHeight, MeanGrayLevel, GrayLevelVariance}) {
            requireNear(grayLevelFromAltitudeSpan.second[grayLevelFromAltitudeSpan.first.linearIndex(nodeId, attribute)],
                        baseline.second[baseline.first.linearIndex(nodeId, attribute)], 1.0e-5f,
                        label + " altitude-span GRAY_LEVEL group must match current path for " + AttributeNames::toString(attribute));
        }
    }

    if constexpr (std::is_floating_point_v<T>) {
        const auto areaComputed = AttributeComputation::computeSingleAttribute(valuedTree, Area);
        const std::vector<T> fractionalAltitude = makeGenericAltitude<T>(valuedTree);
        const ValuedMorphologicalTreeView<T> fractionalView(tree, std::span<const T>(fractionalAltitude));
        const auto fractional = AttributeComputation::computeAttributesFromAltitudeSpan(fractionalView, requests);

        for (NodeId nodeId : tree.aliveNodeIds()) {
            const float area = areaComputed.second[areaComputed.first.linearIndex(nodeId, Area)];
            requireNear(fractional.second[fractional.first.linearIndex(nodeId, Volume)],
                        baseline.second[baseline.first.linearIndex(nodeId, Volume)] + (area * 0.25f), 1.0e-5f,
                        label + " altitude-span API fractional VOLUME must preserve altitude contribution");
            requireNear(fractional.second[fractional.first.linearIndex(nodeId, RelativeVolume)],
                        baseline.second[baseline.first.linearIndex(nodeId, RelativeVolume)], 1.0e-5f,
                        label + " altitude-span API uniform fractional offset must preserve RELATIVE_VOLUME");
            requireNear(fractional.second[fractional.first.linearIndex(nodeId, MeanGrayLevel)],
                        baseline.second[baseline.first.linearIndex(nodeId, MeanGrayLevel)] + 0.25f, 1.0e-5f,
                        label + " altitude-span API fractional MeanGrayLevel must preserve altitude");
            requireNear(fractional.second[fractional.first.linearIndex(nodeId, GrayLevelVariance)],
                        baseline.second[baseline.first.linearIndex(nodeId, GrayLevelVariance)], 1.0e-4f,
                        label + " altitude-span API uniform fractional offset must preserve GrayLevelVariance");
            requireNear(fractional.second[fractional.first.linearIndex(nodeId, GrayLevelHeight)], baseline.second[baseline.first.linearIndex(nodeId, GrayLevelHeight)],
                        1.0e-5f, label + " altitude-span API uniform fractional offset must preserve GrayLevelHeight");
            requireNear(fractional.second[fractional.first.linearIndex(nodeId, MaxDistExact)], baseline.second[baseline.first.linearIndex(nodeId, MaxDistExact)],
                        1.0e-5f, label + " altitude-span API uniform fractional offset must preserve MAX_DIST_EXACT");
        }
    }
}

template <class T> void checkTopologyOnlyAttributePipelineApi(const ValuedMorphologicalTree<std::uint8_t>& valuedTree, const std::string& label) {
    const auto& tree = valuedTree.topology();
    const std::vector<T> equivalentAltitude = makeEquivalentAltitude<T>(valuedTree);
    const ValuedMorphologicalTreeView<T> equivalentView(tree, std::span<const T>(equivalentAltitude));

    for (AttributeGroup group : {AttributeGroup::GrayLevel, AttributeGroup::Shape, AttributeGroup::Moments, AttributeGroup::Boundary,
                                 AttributeGroup::TreeTopology, AttributeGroup::All}) {
        const auto baseline = AttributeComputation::computeAttributes(valuedTree, {group});
        const auto fromAltitudeSpan = AttributeComputation::computeAttributesFromAltitudeSpan(equivalentView, {group});
        requireComputedAttributesEquivalent(fromAltitudeSpan, baseline, tree, ATTRIBUTE_GROUPS.at(group),
                                            label + " attribute-pipeline group must match current path");
    }

    const auto singleBoundaryBaseline = AttributeComputation::computeSingleAttribute(valuedTree, AttributeGroup::Boundary);
    const auto singleBoundaryFromAltitudeSpan = AttributeComputation::computeSingleAttribute(equivalentView, AttributeGroup::Boundary);
    requireComputedAttributesEquivalent(singleBoundaryFromAltitudeSpan, singleBoundaryBaseline, tree, ATTRIBUTE_GROUPS.at(AttributeGroup::Boundary),
                                        label + " attribute-pipeline single group overload must support BOUNDARY");

    requireImageNear(AttributeComputation::computeAttributeMapping(equivalentView, ContourPerimeter),
                     AttributeComputation::computeAttributeMapping(valuedTree, ContourPerimeter),
                     label + " altitude-span view mapping overload must support contour attributes");
}

void checkAttributePipelineMaxDistAndRejectsInvalidInputs(const ValuedMorphologicalTree<std::uint8_t>& valuedTree) {
    const auto& tree = valuedTree.topology();
    const std::vector<float> altitude = makeEquivalentAltitude<float>(valuedTree);

    const auto baselineMaxDist = AttributeComputation::computeSingleAttribute(valuedTree, MaxDistExact);
    const ValuedMorphologicalTreeView<float> view(tree, std::span<const float>(altitude));
    const auto maxDistFromAltitudeSpan = AttributeComputation::computeAttributesFromAltitudeSpan(view, {MaxDistExact});
    requireComputedAttributesNear(maxDistFromAltitudeSpan, baselineMaxDist, tree, {MaxDistExact},
                                  "topology-only MAX_DIST_EXACT must be invariant across valued-tree wrappers");

    const std::vector<std::int32_t> int32Altitude = makeEquivalentAltitude<std::int32_t>(valuedTree);
    const ValuedMorphologicalTreeView<std::int32_t> int32View(tree, std::span<const std::int32_t>(int32Altitude));
    const auto maxDistFromInt32AltitudeSpan = AttributeComputation::computeAttributesFromAltitudeSpan(int32View, {MaxDistExact});
    requireComputedAttributesNear(maxDistFromInt32AltitudeSpan, baselineMaxDist, tree, {MaxDistExact},
                                  "topology-only MAX_DIST_EXACT must ignore the altitude scalar type");

    if constexpr (contract::validationsEnabled) {
        requireThrows<std::runtime_error>(
            [&]() {
                const std::vector<float> wrongSize(static_cast<std::size_t>(tree.numInternalNodeSlots() - 1), 0.0f);
                static_cast<void>(ValuedMorphologicalTreeView<float>(tree, std::span<const float>(wrongSize)));
            },
            "valued-tree view API must reject wrong altitude size");
    }
}

int main() {
    auto image = makeComponentTreeFixture();
    auto valuedTree = makeValuedComponentTree(image, true);
    auto minValuedTree = makeValuedComponentTree(image, false);
    const NodeId sampleNodeId = findNonRootResidueSample(*valuedTree);

    checkRootResidueUsesFixedZeroBaseline();

    checkGenericStaticAltitudeAccess<std::uint8_t>(*valuedTree, sampleNodeId);
    checkGenericStaticAltitudeAccess<std::int32_t>(*valuedTree, sampleNodeId);
    checkGenericStaticAltitudeAccess<float>(*valuedTree, sampleNodeId);
    checkGenericStaticAltitudeAccess<double>(*valuedTree, sampleNodeId);
    checkValuedMorphologicalTreeViewContract<std::uint8_t>(*valuedTree, sampleNodeId, "uint8");
    checkValuedMorphologicalTreeViewContract<std::int32_t>(*valuedTree, sampleNodeId, "int32");
    checkValuedMorphologicalTreeViewContract<float>(*valuedTree, sampleNodeId, "float");
    checkValuedMorphologicalTreeViewContract<double>(*valuedTree, sampleNodeId, "double");
    const auto ownedView = valuedTree->asView();
    require(&ownedView.topology() == &valuedTree->topology(), "ValuedMorphologicalTree<std::uint8_t>::asView must borrow the owned topology");
    requireEqual(static_cast<int>(ownedView.nodeAltitudes().size()), valuedTree->topology().numInternalNodeSlots(),
                 "ValuedMorphologicalTree<std::uint8_t>::asView altitude span size");
    checkValuedMorphologicalTreeViewIncrementalAttributePipeline(*valuedTree);
    checkFloatResidueAndMonotonicity(*valuedTree, sampleNodeId);
    checkUnsignedResidueDoesNotUnderflow(*valuedTree, sampleNodeId);
    checkFloatHigraExport(*valuedTree);
    checkInt32HigraExport(*valuedTree);
    checkGenericAltitudeSampling<std::int32_t>(*valuedTree, "max-tree int32");
    checkGenericAltitudeSampling<float>(*valuedTree, "max-tree float");
    checkGenericAltitudeSampling<double>(*valuedTree, "max-tree double");
    checkGenericAltitudeSampling<float>(*minValuedTree, "min-tree float");
    checkGenericAltitudeSampling<double>(*minValuedTree, "min-tree double");
    checkGenericVolumeKernel<float>(*valuedTree, "max-tree float");
    checkGenericVolumeKernel<double>(*valuedTree, "max-tree double");
    checkGenericGrayLevelStatsKernel<float>(*valuedTree, "max-tree float");
    checkGenericGrayLevelStatsKernel<double>(*valuedTree, "max-tree double");
    checkGenericGrayLevelStatsKernel<float>(*minValuedTree, "min-tree float");
    checkGenericGrayLevelStatsKernel<double>(*minValuedTree, "min-tree double");
    checkAttributePipelineApi<float>(*valuedTree, "max-tree float");
    checkAttributePipelineApi<double>(*valuedTree, "max-tree double");
    checkAttributePipelineApi<float>(*minValuedTree, "min-tree float");
    checkAttributePipelineApi<double>(*minValuedTree, "min-tree double");
    checkTopologyOnlyAttributePipelineApi<std::int32_t>(*valuedTree, "max-tree int32");
    checkTopologyOnlyAttributePipelineApi<float>(*valuedTree, "max-tree float");
    checkTopologyOnlyAttributePipelineApi<double>(*minValuedTree, "min-tree double");
    checkAttributePipelineMaxDistAndRejectsInvalidInputs(*valuedTree);
    checkImageFactoriesDeduceAltitudeType(image);
    checkTypedOwnerExportedHigraAttributeProjection<std::int32_t>("int32 owner");
    checkTypedOwnerExportedHigraAttributeProjection<float>("float owner");
    checkTypedOwnerSampledAttributeApi<std::int32_t>("int32 owner");
    checkTypedOwnerSampledAttributeApi<float>("float owner");
    checkTypedHigraImportFactory(image);
    if constexpr (contract::validationsEnabled) {
        checkFiniteFloatAltitudeValidation(image);
    }

    return 0;
}
