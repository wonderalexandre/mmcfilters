#include "support/TestSupport.hpp"

#include "mmcfilters/trees/saliency/HierarchySaliencyMap.hpp"
#include "mmcfilters/trees/sdrt/SaturatedResidualTreeBuilder.hpp"
#include "mmcfilters/trees/sdrt/UnrestrictedResidualTreeBuilder.hpp"
#include "mmcfilters/trees/sdrt/detail/FlatZonePartition.hpp"
#include "mmcfilters/trees/sdrt/detail/ResidualTreeEventAssembler.hpp"
#include "mmcfilters/trees/sdrt/detail/ResidualTreeRegionTypes.hpp"
#include "sdrt_reference/OptimizedUnionFindSelfDualResidualTreeBuilder.hpp"
#include "sdrt_reference/detail/RecomputeFromSupportReference.hpp"
#include "sdrt_reference/oracle/rag/SingleAdjacencySaturatedResidualTreeOracle.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <span>
#include <type_traits>
#include <utility>
#include <vector>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

namespace {

template <class Engine>
concept HasInfinityPixelAccessor = requires(const Engine& engine) { engine.getInfinityPixel(); };

template <class Assembler>
concept LvalueFinalizableAssembler = requires(Assembler& assembler) { assembler.finalize(sdrt::detail::RegionId{0}, 0); };

template <class Assembler>
concept RvalueFinalizableAssembler = requires(Assembler&& assembler) { std::move(assembler).finalize(sdrt::detail::RegionId{0}, 0); };

static_assert(!std::is_same_v<sdrt::detail::MinMaxResidualTreeEngine<std::uint8_t, true>, sdrt::detail::MinMaxResidualTreeEngine<std::uint8_t, false>>,
              "saturated and unrestricted builders must use distinct engine specializations");
static_assert(HasInfinityPixelAccessor<sdrt::detail::MinMaxResidualTreeEngine<std::uint8_t, true>>,
              "saturated construction must expose its configured exterior seed");
static_assert(!HasInfinityPixelAccessor<sdrt::detail::MinMaxResidualTreeEngine<std::uint8_t, false>>,
              "unrestricted construction must not expose exterior-seed state");
static_assert(!LvalueFinalizableAssembler<sdrt::detail::ResidualTreeEventAssembler<int>>,
              "residual event assembly must not be finalizable through an lvalue");
static_assert(RvalueFinalizableAssembler<sdrt::detail::ResidualTreeEventAssembler<int>>,
              "residual event assembly must be finalizable through an rvalue");
static_assert(!std::is_constructible_v<sdrt::UnrestrictedResidualTreeBuilder<std::uint8_t>, RegularGridAdjacency2D, NodeId, sdrt::SaturatedResidualTreeOptions>,
              "unrestricted builders must not accept an exterior seed or saturated policies");
static_assert(!std::is_constructible_v<sdrt::UnrestrictedResidualTreeBuilder<std::uint8_t>, RegularGridAdjacency2D, sdrt::SdrtTiePolicy>,
              "unrestricted builders must require their mode-specific option object");
static_assert(!std::is_constructible_v<sdrt::SaturatedResidualTreeBuilder<std::uint8_t>, RegularGridAdjacency2D, NodeId, sdrt::SdrtTiePolicy>,
              "saturated builders must require their mode-specific option object");

template <class T> std::vector<T> toVector(std::span<const T> values) { return {values.begin(), values.end()}; }

template <class T> std::vector<NodeId> parentBuffer(const WeightedMorphologicalTree<T>& tree) {
    std::vector<NodeId> parent(static_cast<std::size_t>(tree.topology().getNumInternalNodeSlots()));
    for (NodeId node : tree.topology().getAliveNodeIds()) {
        parent[static_cast<std::size_t>(node)] = tree.topology().getNodeParent(node);
    }
    return parent;
}

template <class T> std::vector<NodeId> ownerBuffer(const WeightedMorphologicalTree<T>& tree) {
    std::vector<NodeId> owner(static_cast<std::size_t>(tree.topology().getNumTotalProperParts()));
    for (NodeId pixel = 0; pixel < tree.topology().getNumTotalProperParts(); ++pixel) {
        owner[static_cast<std::size_t>(pixel)] = tree.topology().getProperPartOwner(pixel);
    }
    return owner;
}

template <class T> std::vector<T> altitudeBuffer(const WeightedMorphologicalTree<T>& tree) {
    std::vector<T> altitude(static_cast<std::size_t>(tree.topology().getNumInternalNodeSlots()));
    for (NodeId node : tree.topology().getAliveNodeIds()) {
        altitude[static_cast<std::size_t>(node)] = tree.getAltitude(node);
    }
    return altitude;
}

template <class T> struct ResidualBuffers {
    std::vector<NodeId> parent;
    std::vector<NodeId> owner;
    std::vector<T> altitude;
};

template <class T>
ResidualBuffers<T> buildResidualWithPolicies(const ImagePtr<T>& image, const RegularGridAdjacency2D& adjacency, NodeId infinityPixel,
                                             sdrt::SaturatedMinMaxLcaPolicy lcaPolicy, sdrt::SaturatedMinMaxFallbackPolicy fallbackPolicy) {
    auto minTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto maxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
    const sdrt::SaturatedResidualTreeOptions options{sdrt::SdrtTiePolicy::ContrastInvariantSpatial, lcaPolicy, fallbackPolicy};
    sdrt::SaturatedResidualTreeBuilder<T> builder(adjacency, infinityPixel, options);
    builder.build(image, std::move(minTree), std::move(maxTree));
    return {toVector(builder.getNodeParent()), toVector(builder.getProperPartOwner()), toVector(builder.getAltitude())};
}

ImageUInt8Ptr makeRadixImage(int code) {
    auto image = ImageUInt8::create(2, 3);
    for (NodeId pixel = 0; pixel < 6; ++pixel) {
        (*image)[pixel] = static_cast<std::uint8_t>(code % 3);
        code /= 3;
    }
    return image;
}

template <class T> ImagePtr<T> makeTypedImage(int rows, int cols, std::initializer_list<T> values) {
    requireEqual(static_cast<int>(values.size()), rows * cols, "typed image buffer size");
    auto image = Image<T>::create(rows, cols);
    NodeId pixel = 0;
    for (T value : values) {
        (*image)[pixel++] = value;
    }
    return image;
}

template <class T> void requireExactReconstruction(const WeightedMorphologicalTree<T>& tree, const ImagePtr<T>& image, const std::string& label) {
    for (NodeId pixel = 0; pixel < image->getSize(); ++pixel) {
        const NodeId owner = tree.topology().getProperPartOwner(pixel);
        requireEqual(tree.getAltitude(owner), (*image)[pixel], label);
    }
}

void requireResidualSemantics(const MorphologicalTree& tree, const RegularGridAdjacency2D& adjacency, const std::string& label) {
    require(tree.getDescriptiveKind() == MorphologicalTreeKind::SELF_DUAL_RESIDUAL_TREE, label + ": descriptive kind");
    require(tree.getAltitudeOrder() == AltitudeOrder::UNCONSTRAINED, label + ": unconstrained altitude order");
    const auto* stored = tree.getUniformGridAdjacency2D();
    require(stored != nullptr, label + ": uniform adjacency");
    requireEqual(stored->getNumRows(), adjacency.getNumRows(), label + ": adjacency rows");
    requireEqual(stored->getNumCols(), adjacency.getNumCols(), label + ": adjacency columns");
}

void testFactoryModesAndSaliencyIntegration() {
    const auto image = makeImage(3, 4,
                                 {
                                     0,
                                     1,
                                     2,
                                     1,
                                     1,
                                     4,
                                     1,
                                     0,
                                     2,
                                     1,
                                     3,
                                     1,
                                 });
    const RegularGridAdjacency2D adjacency(3, 4, 1.0);

    const auto unrestricted = MorphologicalTreeFactory::createSelfDualResidualTree(image, adjacency);
    const auto saturated = MorphologicalTreeFactory::createSaturatedSelfDualResidualTree(image, adjacency, NodeId{0});

    requireExactReconstruction(unrestricted, image, "unrestricted residual reconstruction");
    requireExactReconstruction(saturated, image, "saturated residual reconstruction");
    requireResidualSemantics(unrestricted.topology(), adjacency, "unrestricted residual semantics");
    requireResidualSemantics(saturated.topology(), adjacency, "saturated residual semantics");

    const auto unrestrictedSaliency = HierarchySaliencyMap::computeTopologicalLevelEdgeMap(unrestricted.topology());
    const auto saturatedSaliency = HierarchySaliencyMap::computeTopologicalLevelEdgeMap(saturated.topology());
    require(!unrestrictedSaliency.empty(), "unrestricted residual topological saliency");
    require(!saturatedSaliency.empty(), "saturated residual topological saliency");
    requireThrows<std::invalid_argument>([&] { static_cast<void>(HierarchySaliencyMap::computeNormalizedAltitudeEdgeMap(unrestricted)); },
                                         "residual gray-level altitude is not a hierarchy saliency valuation");
}

void testContrastInversion() {
    const auto image = makeImage(1, 6, {0, 2, 0, 2, 1, 3});
    const auto inverted = makeImage(1, 6, {3, 1, 3, 1, 2, 0});
    const RegularGridAdjacency2D adjacency(1, 6, 1.0);

    const auto primary = MorphologicalTreeFactory::createSaturatedSelfDualResidualTree(image, adjacency, NodeId{0});
    const auto conjugate = MorphologicalTreeFactory::createSaturatedSelfDualResidualTree(inverted, adjacency, NodeId{0});

    const MorphologicalTree& lhs = primary.topology();
    const MorphologicalTree& rhs = conjugate.topology();
    requireEqual(lhs.getNumNodes(), rhs.getNumNodes(), "contrast node count");
    for (NodeId node : lhs.getAliveNodeIds()) {
        requireEqual(lhs.getNodeParent(node), rhs.getNodeParent(node), "contrast parent");
        requireEqual(static_cast<int>(primary.getAltitude(node)) + static_cast<int>(conjugate.getAltitude(node)), 3, "contrast altitude");
    }
    for (NodeId pixel = 0; pixel < image->getSize(); ++pixel) {
        requireEqual(lhs.getProperPartOwner(pixel), rhs.getProperPartOwner(pixel), "contrast proper-part owner");
    }
}

void testDirectBuilderPolicies() {
    const auto image = makeImage(2, 3,
                                 {
                                     0,
                                     2,
                                     1,
                                     2,
                                     0,
                                     0,
                                 });
    const RegularGridAdjacency2D adjacency(2, 3, 1.0);
    auto minTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto maxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
    sdrt::UnrestrictedResidualTreeBuilder<std::uint8_t> builder(adjacency,
                                                                sdrt::UnrestrictedResidualTreeOptions{sdrt::SdrtTiePolicy::ContrastInvariantSpatial});
    builder.build(image, std::move(minTree), std::move(maxTree));
    requireEqual(builder.getStatistics().rejectedExtrema, std::size_t{0}, "unrestricted builder rejects no extrema");
    requireEqual(builder.getStatistics().complementTraversalCertificates, std::size_t{0}, "unrestricted builder performs no complement traversal");
}

void testResidualTreeEventAssemblerInIsolation() {
    using sdrt::detail::RegionId;
    const std::vector<RegionId> regionByPixel{0, 0, 1, 2};
    sdrt::detail::ResidualTreeEventAssembler<int> assembler(3, regionByPixel);
    sdrt::reference::PersistentSupportHistory<int> reference(4);

    requireEqual(assembler.emitEvent(RegionId{0}, 5), NodeId{1}, "assembler first event id");
    const std::array<NodeId, 2> firstSupport{0, 1};
    reference.recordEvent(5, firstSupport);
    const std::array<RegionId, 1> firstAbsorption{0};
    assembler.consume(RegionId{1}, firstAbsorption);
    requireEqual(assembler.emitEvent(RegionId{1}, 3), NodeId{2}, "assembler second event id");
    const std::array<NodeId, 3> secondSupport{0, 1, 2};
    reference.recordEvent(3, secondSupport);
    const std::array<RegionId, 1> secondAbsorption{1};
    assembler.consume(RegionId{2}, secondAbsorption);
    requireEqual(assembler.numEvents(), std::size_t{2}, "assembler event count");

    auto output = std::move(assembler).finalize(RegionId{2}, 1);
    const auto referenceOutput = reference.finalize(1);
    requireVectorEqual(output.nodeParent, referenceOutput.nodeParent, "assembler reference parents");
    requireVectorEqual(output.properPartOwner, referenceOutput.properPartOwner, "assembler reference owners");
    requireVectorEqual(output.altitude, referenceOutput.altitude, "assembler reference altitudes");
    requireVectorEqual(output.nodeParent, std::vector<NodeId>{0, 2, 0}, "assembler parent materialization");
    requireVectorEqual(output.properPartOwner, std::vector<NodeId>{1, 1, 2, 0}, "assembler owner materialization");
    requireVectorEqual(output.altitude, std::vector<int>{1, 5, 3}, "assembler altitude materialization");

    sdrt::detail::ResidualTreeEventAssembler<int> invalidAssembler(2, std::vector<RegionId>{0, 1});
    const std::array<RegionId, 1> absorbsSelf{0};
    requireThrows<std::invalid_argument>([&] { invalidAssembler.consume(RegionId{0}, absorbsSelf); }, "assembler rejects absorbing its survivor");
}

void testIncrementalBoundaryAgainstRecomputedReference() {
    const auto image = makeImage(3, 4,
                                 {
                                     0, 0, 2, 1,
                                     0, 3, 2, 1,
                                     4, 3, 3, 1,
                                 });
    const RegularGridAdjacency2D adjacency(3, 4, 1.0);
    const auto verifyTree = [&](const WeightedMorphologicalTree<std::uint8_t>& tree, const std::string& label) {
        sdrt::detail::FlatZonePartition<std::uint8_t> partition;
        sdrt::reference::RecomputedProperPartBoundary<std::uint8_t> reference;
        partition.initialize(image, adjacency, NodeId{0});

        for (NodeId node : tree.topology().getAliveNodeIds()) {
            if (tree.topology().getNumProperParts(node) == 0) {
                continue;
            }
            const auto partitionView = partition.viewForNode(node, tree);
            std::vector<NodeId> incrementalSupport(partitionView.supportPixels.begin(), partitionView.supportPixels.end());
            std::vector<NodeId> incrementalBoundary(partitionView.externalPixelsByIncidence.begin(), partitionView.externalPixelsByIncidence.end());
            const auto& recomputed = reference.acquire(node, tree, adjacency);
            auto recomputedSupport = recomputed.supportPixels;
            auto recomputedBoundary = recomputed.externalPixelsByIncidence;
            std::ranges::sort(incrementalSupport);
            std::ranges::sort(incrementalBoundary);
            std::ranges::sort(recomputedSupport);
            std::ranges::sort(recomputedBoundary);
            requireVectorEqual(incrementalSupport, recomputedSupport, label + " support");
            requireVectorEqual(incrementalBoundary, recomputedBoundary, label + " boundary incidences");
        }
    };

    verifyTree(MorphologicalTreeFactory::createMaxTree(image, adjacency), "max-tree incremental boundary");
    verifyTree(MorphologicalTreeFactory::createMinTree(image, adjacency), "min-tree incremental boundary");

    const auto chainImage = makeImage(1, 4, {2, 1, 1, 0});
    const RegularGridAdjacency2D chainAdjacency(1, 4, 1.0);
    sdrt::detail::FlatZonePartition<std::uint8_t> partition;
    partition.initialize(chainImage, chainAdjacency, NodeId{3});
    GenerationStampSet representativeMarks(static_cast<std::size_t>(chainImage->getSize()));
    std::vector<NodeId> mergeRepresentatives;

    const auto requireViewMatchesSupport = [&](const auto& view, std::vector<NodeId> expectedSupport, const std::string& label) {
        std::vector<std::uint8_t> supportMarks(static_cast<std::size_t>(chainImage->getSize()), std::uint8_t{0});
        for (NodeId pixel : expectedSupport) {
            supportMarks[static_cast<std::size_t>(pixel)] = 1;
        }
        std::vector<NodeId> expectedBoundary;
        for (NodeId pixel : expectedSupport) {
            for (NodeId neighbor : chainAdjacency.getNeighborIndices(pixel)) {
                if (supportMarks[static_cast<std::size_t>(neighbor)] == 0) {
                    expectedBoundary.push_back(neighbor);
                }
            }
        }
        std::vector<NodeId> actualSupport(view.supportPixels.begin(), view.supportPixels.end());
        std::vector<NodeId> actualBoundary(view.externalPixelsByIncidence.begin(), view.externalPixelsByIncidence.end());
        std::ranges::sort(expectedSupport);
        std::ranges::sort(expectedBoundary);
        std::ranges::sort(actualSupport);
        std::ranges::sort(actualBoundary);
        requireVectorEqual(actualSupport, expectedSupport, label + " support");
        requireVectorEqual(actualBoundary, expectedBoundary, label + " boundary");
    };

    auto selectedView = partition.viewForPixel(NodeId{0});
    requireViewMatchesSupport(selectedView, {0}, "initial selected flat zone");
    partition.collectAdjacentRepresentativesAtLevel(selectedView.representative, std::uint8_t{1}, selectedView.externalPixelsByIncidence,
                                                    representativeMarks, mergeRepresentatives);
    requireEqual(mergeRepresentatives.size(), std::size_t{1}, "first incremental flat-zone merge count");
    NodeId mergedRepresentative = partition.mergeFlatZonesAtLevel(selectedView.representative, mergeRepresentatives, std::uint8_t{1});
    auto mergedView = partition.viewForPixel(mergedRepresentative);
    requireViewMatchesSupport(mergedView, {0, 1, 2}, "first merged flat zone");
    requireEqual(mergedView.stableSpatialKey, NodeId{0}, "merged flat-zone stable spatial key");
    require(!mergedView.containsExteriorSeed, "first merged flat zone excludes exterior seed");

    partition.collectAdjacentRepresentativesAtLevel(mergedView.representative, std::uint8_t{0}, mergedView.externalPixelsByIncidence,
                                                    representativeMarks, mergeRepresentatives);
    requireEqual(mergeRepresentatives.size(), std::size_t{1}, "second incremental flat-zone merge count");
    mergedRepresentative = partition.mergeFlatZonesAtLevel(mergedView.representative, mergeRepresentatives, std::uint8_t{0});
    const auto terminalView = partition.viewForPixel(mergedRepresentative);
    requireViewMatchesSupport(terminalView, {0, 1, 2, 3}, "terminal merged flat zone");
    require(terminalView.containsExteriorSeed, "terminal merged flat zone retains exterior seed");
    requireEqual(partition.representativeOf(NodeId{3}), mergedRepresentative, "terminal flat-zone representative");
}

template <class T> void testTypedResidualConstruction(const std::string& label) {
    const auto image = makeTypedImage<T>(2, 3, {T{-3}, T{2}, T{1}, T{4}, T{-1}, T{0}});
    const RegularGridAdjacency2D adjacency(2, 3, 1.0);
    const auto unrestricted = MorphologicalTreeFactory::createSelfDualResidualTree(image, adjacency);
    const auto saturated = MorphologicalTreeFactory::createSaturatedSelfDualResidualTree(image, adjacency, NodeId{0});
    requireExactReconstruction(unrestricted, image, label + " unrestricted reconstruction");
    requireExactReconstruction(saturated, image, label + " saturated reconstruction");
}

void testAltitudeTypesAndDegenerateDomain() {
    testTypedResidualConstruction<std::int16_t>("int16");
    testTypedResidualConstruction<std::int32_t>("int32");
    testTypedResidualConstruction<float>("float");
    testTypedResidualConstruction<double>("double");

    const auto unsignedImage = makeTypedImage<std::uint16_t>(2, 3, {0, 5, 2, 7, 1, 3});
    const RegularGridAdjacency2D unsignedAdjacency(2, 3, 1.0);
    requireExactReconstruction(MorphologicalTreeFactory::createSelfDualResidualTree(unsignedImage, unsignedAdjacency), unsignedImage,
                               "uint16 unrestricted reconstruction");
    requireExactReconstruction(MorphologicalTreeFactory::createSaturatedSelfDualResidualTree(unsignedImage, unsignedAdjacency, NodeId{0}), unsignedImage,
                               "uint16 saturated reconstruction");

    const auto singleton = makeTypedImage<double>(1, 1, {2.5});
    const auto unrestrictedSingleton = MorphologicalTreeFactory::createSelfDualResidualTree(singleton, 1.0);
    const auto saturatedSingleton = MorphologicalTreeFactory::createSaturatedSelfDualResidualTree(singleton, NodeId{0}, 1.0);
    requireEqual(unrestrictedSingleton.topology().getNumNodes(), 1, "unrestricted singleton node count");
    requireEqual(saturatedSingleton.topology().getNumNodes(), 1, "saturated singleton node count");
    requireExactReconstruction(unrestrictedSingleton, singleton, "unrestricted singleton reconstruction");
    requireExactReconstruction(saturatedSingleton, singleton, "saturated singleton reconstruction");
}

void testOptionsAndInputContracts() {
    const auto image = makeImage(2, 3, {0, 2, 1, 2, 0, 0});
    const RegularGridAdjacency2D adjacency(2, 3, 1.0);

    sdrt::UnrestrictedResidualTreeOptions unrestrictedOptions;
    unrestrictedOptions.tiePolicy = sdrt::SdrtTiePolicy::MaxBeforeMinThenSpatial;
    const auto unrestrictedWithOptions = MorphologicalTreeFactory::createSelfDualResidualTree(image, adjacency, unrestrictedOptions);
    auto optionMinTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto optionMaxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
    sdrt::UnrestrictedResidualTreeBuilder<std::uint8_t> unrestrictedBuilder(adjacency, unrestrictedOptions);
    unrestrictedBuilder.build(image, std::move(optionMinTree), std::move(optionMaxTree));
    requireVectorEqual(parentBuffer(unrestrictedWithOptions), toVector(unrestrictedBuilder.getNodeParent()), "unrestricted option parents");
    requireVectorEqual(ownerBuffer(unrestrictedWithOptions), toVector(unrestrictedBuilder.getProperPartOwner()), "unrestricted option owners");
    requireVectorEqual(altitudeBuffer(unrestrictedWithOptions), toVector(unrestrictedBuilder.getAltitude()), "unrestricted option altitudes");

    sdrt::SaturatedResidualTreeOptions saturatedOptions;
    saturatedOptions.tiePolicy = sdrt::SdrtTiePolicy::ContrastInvariantSpatial;
    saturatedOptions.lcaPolicy = sdrt::SaturatedMinMaxLcaPolicy::LinkCut;
    saturatedOptions.fallbackPolicy = sdrt::SaturatedMinMaxFallbackPolicy::SingleSourceDepthFirst;
    const auto saturatedWithOptions = MorphologicalTreeFactory::createSaturatedSelfDualResidualTree(image, adjacency, NodeId{0}, saturatedOptions);
    const auto saturatedDefault = MorphologicalTreeFactory::createSaturatedSelfDualResidualTree(image, adjacency, NodeId{0});
    requireVectorEqual(parentBuffer(saturatedWithOptions), parentBuffer(saturatedDefault), "saturated option parents");
    requireVectorEqual(ownerBuffer(saturatedWithOptions), ownerBuffer(saturatedDefault), "saturated option owners");
    requireVectorEqual(altitudeBuffer(saturatedWithOptions), altitudeBuffer(saturatedDefault), "saturated option altitudes");

    sdrt::SaturatedResidualTreeBuilder<std::uint8_t> invalidInfinity(adjacency, NodeId{6});
    auto minTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto maxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
    requireThrows<std::invalid_argument>([&] { invalidInfinity.build(image, std::move(minTree), std::move(maxTree)); },
                                         "saturated builder rejects an out-of-domain infinity pixel");

    sdrt::UnrestrictedResidualTreeBuilder<std::uint8_t> reusable(adjacency);
    auto firstMinTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto secondMinTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    requireThrows<std::invalid_argument>([&] { reusable.build(image, std::move(firstMinTree), std::move(secondMinTree)); },
                                         "unrestricted builder rejects a min-tree in the max-tree slot");
    requireThrows<std::logic_error>([&] { static_cast<void>(reusable.getNodeParent()); }, "failed construction leaves no observable partial result");
    auto validMinTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto validMaxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
    reusable.build(image, std::move(validMinTree), std::move(validMaxTree));
    requireEqual(reusable.getProperPartOwner().size(), static_cast<std::size_t>(image->getSize()), "builder reuse after failure");

    const std::array<GridOffset2D, 1> originOnly{{{0, 0}}};
    const auto disconnected = RegularGridAdjacency2D::fromStructuringElement(2, 3, originOnly);
    sdrt::UnrestrictedResidualTreeBuilder<std::uint8_t> disconnectedBuilder(disconnected);
    auto connectedMinTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto connectedMaxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
    requireThrows<std::invalid_argument>([&] { disconnectedBuilder.build(image, std::move(connectedMinTree), std::move(connectedMaxTree)); },
                                         "residual builder rejects a disconnected adjacency");

    if constexpr (contract::validationsEnabled) {
        auto nonFinite = makeTypedImage<float>(1, 2, {0.0F, std::numeric_limits<float>::quiet_NaN()});
        requireThrows<std::invalid_argument>([&] { static_cast<void>(MorphologicalTreeFactory::createSelfDualResidualTree(nonFinite, 1.0)); },
                                             "unrestricted factory rejects NaN");
        (*nonFinite)[1] = std::numeric_limits<float>::infinity();
        requireThrows<std::invalid_argument>(
            [&] { static_cast<void>(MorphologicalTreeFactory::createSaturatedSelfDualResidualTree(nonFinite, NodeId{0}, 1.0)); },
            "saturated factory rejects infinity");
    }
}

void testExhaustiveDifferentialOracles() {
    using UnrestrictedReference = sdrt::OptimizedUnionFindSelfDualResidualTreeBuilder<std::uint8_t>;
    using SaturatedReference = sdrt::oracle::rag::SingleAdjacencySaturatedResidualTreeOracle<std::uint8_t>;

    const RegularGridAdjacency2D adjacency(2, 3, 1.0);
    for (int code = 0; code < 729; ++code) {
        const auto image = makeRadixImage(code);

        const auto unrestricted = MorphologicalTreeFactory::createSelfDualResidualTree(image, adjacency);
        UnrestrictedReference unrestrictedReference(1.0, sdrt::SdrtTiePolicy::ContrastInvariantSpatial);
        unrestrictedReference.build(image);
        requireVectorEqual(parentBuffer(unrestricted), toVector(unrestrictedReference.getNodeParent()), "exhaustive unrestricted parents");
        requireVectorEqual(ownerBuffer(unrestricted), toVector(unrestrictedReference.getProperPartOwner()), "exhaustive unrestricted owners");
        requireVectorEqual(altitudeBuffer(unrestricted), toVector(unrestrictedReference.getAltitude()), "exhaustive unrestricted altitudes");

        const auto polarityOrdered = MorphologicalTreeFactory::createSelfDualResidualTree(
            image, adjacency, sdrt::UnrestrictedResidualTreeOptions{sdrt::SdrtTiePolicy::MaxBeforeMinThenSpatial});
        UnrestrictedReference polarityOrderedReference(1.0, sdrt::SdrtTiePolicy::MaxBeforeMinThenSpatial);
        polarityOrderedReference.build(image);
        requireVectorEqual(parentBuffer(polarityOrdered), toVector(polarityOrderedReference.getNodeParent()), "exhaustive polarity-ordered parents");
        requireVectorEqual(ownerBuffer(polarityOrdered), toVector(polarityOrderedReference.getProperPartOwner()), "exhaustive polarity-ordered owners");
        requireVectorEqual(altitudeBuffer(polarityOrdered), toVector(polarityOrderedReference.getAltitude()), "exhaustive polarity-ordered altitudes");

        for (NodeId infinityPixel = 0; infinityPixel < 6; ++infinityPixel) {
            const auto saturated = MorphologicalTreeFactory::createSaturatedSelfDualResidualTree(image, adjacency, infinityPixel);
            const auto saturatedReference = SaturatedReference::build(image, adjacency, infinityPixel);
            requireVectorEqual(parentBuffer(saturated), saturatedReference.nodeParent, "exhaustive saturated parents");
            requireVectorEqual(ownerBuffer(saturated), saturatedReference.properPartOwner, "exhaustive saturated owners");
            requireVectorEqual(altitudeBuffer(saturated), saturatedReference.altitude, "exhaustive saturated altitudes");
        }
    }
}

void testSharedArbitrarySymmetricAdjacency() {
    using SaturatedReference = sdrt::oracle::rag::SingleAdjacencySaturatedResidualTreeOracle<std::uint8_t>;
    const std::array<GridOffset2D, 7> offsets{{
        {0, 0},
        {-1, 0},
        {1, 0},
        {0, -1},
        {0, 1},
        {0, -2},
        {0, 2},
    }};
    const RegularGridAdjacency2D adjacency = RegularGridAdjacency2D::fromStructuringElement(4, 5, offsets);
    const auto image = makeImage(4, 5,
                                 {
                                     0, 3, 1, 4, 2, 2, 1, 4, 0, 3, 4, 0, 2, 3, 1, 1, 4, 0, 2, 3,
                                 });
    const auto unrestricted = MorphologicalTreeFactory::createSelfDualResidualTree(image, adjacency);
    const auto saturated = MorphologicalTreeFactory::createSaturatedSelfDualResidualTree(image, adjacency, NodeId{0});
    const auto reference = SaturatedReference::build(image, adjacency, NodeId{0});

    requireExactReconstruction(unrestricted, image, "custom-adjacency unrestricted reconstruction");
    requireVectorEqual(parentBuffer(saturated), reference.nodeParent, "custom-adjacency saturated parents");
    requireVectorEqual(ownerBuffer(saturated), reference.properPartOwner, "custom-adjacency saturated owners");
    requireVectorEqual(altitudeBuffer(saturated), reference.altitude, "custom-adjacency saturated altitudes");
}

void testConfigurablePolicyEquivalence() {
    constexpr std::array lcaPolicies{sdrt::SaturatedMinMaxLcaPolicy::ParentClimb, sdrt::SaturatedMinMaxLcaPolicy::BlockedSnapshot,
                                     sdrt::SaturatedMinMaxLcaPolicy::LinkCut};
    constexpr std::array fallbackPolicies{sdrt::SaturatedMinMaxFallbackPolicy::SingleSourceDepthFirst,
                                          sdrt::SaturatedMinMaxFallbackPolicy::BoundaryMultiSource};
    const RegularGridAdjacency2D adjacency(2, 3, 1.0);

    for (int code = 0; code < 729; ++code) {
        const auto image = makeRadixImage(code);
        const NodeId infinityPixel = static_cast<NodeId>(code % 6);
        const auto expected = MorphologicalTreeFactory::createSaturatedSelfDualResidualTree(image, adjacency, infinityPixel);
        const auto expectedParent = parentBuffer(expected);
        const auto expectedOwner = ownerBuffer(expected);
        const auto expectedAltitude = altitudeBuffer(expected);

        for (const auto lcaPolicy : lcaPolicies) {
            for (const auto fallbackPolicy : fallbackPolicies) {
                const auto actual = buildResidualWithPolicies(image, adjacency, infinityPixel, lcaPolicy, fallbackPolicy);
                requireVectorEqual(actual.parent, expectedParent, "configurable policy parents");
                requireVectorEqual(actual.owner, expectedOwner, "configurable policy owners");
                requireVectorEqual(actual.altitude, expectedAltitude, "configurable policy altitudes");
            }
        }
    }
}

} // namespace

int main() {
    testFactoryModesAndSaliencyIntegration();
    testContrastInversion();
    testDirectBuilderPolicies();
    testResidualTreeEventAssemblerInIsolation();
    testIncrementalBoundaryAgainstRecomputedReference();
    testAltitudeTypesAndDegenerateDomain();
    testOptionsAndInputContracts();
    testExhaustiveDifferentialOracles();
    testSharedArbitrarySymmetricAdjacency();
    testConfigurablePolicyEquivalence();
    return 0;
}
