#include "support/TestSupport.hpp"

#include "mmcfilters/trees/saliency/HierarchySaliencyMap.hpp"
#include "mmcfilters/trees/sdrt/SaturatedResidualTreeBuilder.hpp"
#include "mmcfilters/trees/sdrt/UnrestrictedResidualTreeBuilder.hpp"
#include "mmcfilters/trees/sdrt/detail/FlatZonePartition.hpp"
#include "mmcfilters/trees/sdrt/detail/ResidualTreeCandidatePreparation.hpp"
#include "mmcfilters/trees/sdrt/detail/ResidualTreeEventAssembler.hpp"
#include "mmcfilters/trees/sdrt/detail/ResidualTreeRegionTypes.hpp"
#include "mmcfilters/trees/sdrt/detail/SaturatedResidualEligibility.hpp"

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
concept HasInfinityPixelAccessor = requires(const Engine& engine) { engine.infinityPixel(); };

template <class T>
concept HasComplementarySaturatedResidualFactory = requires(ImagePtr<T> image, ComplementaryAdjacencies adjacencies) {
    MorphologicalTreeFactory::createSaturatedResidualTree(image, adjacencies, PixelId{0});
};

template <class Assembler>
concept LvalueFinalizableAssembler = requires(Assembler& assembler) { assembler.finalize(sdrt::detail::RegionId{0}, 0); };

template <class Assembler>
concept RvalueFinalizableAssembler = requires(Assembler&& assembler) { std::move(assembler).finalize(sdrt::detail::RegionId{0}, 0); };

static_assert(!std::is_same_v<sdrt::detail::SynchronizedResidualTreeEvolution<std::uint8_t, true>,
                              sdrt::detail::SynchronizedResidualTreeEvolution<std::uint8_t, false>>,
              "saturated and unrestricted builders must use distinct engine specializations");
static_assert(HasInfinityPixelAccessor<sdrt::detail::SynchronizedResidualTreeEvolution<std::uint8_t, true>>,
              "saturated construction must expose its configured infinity pixel");
static_assert(!HasInfinityPixelAccessor<sdrt::detail::SynchronizedResidualTreeEvolution<std::uint8_t, false>>,
              "unrestricted construction must not expose infinity-pixel state");
static_assert(!HasComplementarySaturatedResidualFactory<std::uint8_t>,
              "complementary adjacencies must be routed through tree-of-shapes construction");
static_assert(!LvalueFinalizableAssembler<sdrt::detail::ResidualTreeEventAssembler<int>>,
              "residual event assembly must not be finalizable through an lvalue");
static_assert(RvalueFinalizableAssembler<sdrt::detail::ResidualTreeEventAssembler<int>>,
              "residual event assembly must be finalizable through an rvalue");
static_assert(!std::is_constructible_v<sdrt::UnrestrictedResidualTreeBuilder<std::uint8_t>, RegularGridAdjacency2D, NodeId, sdrt::SaturatedResidualTreeOptions>,
              "unrestricted builders must not accept an infinity pixel or saturated policies");
static_assert(!std::is_constructible_v<sdrt::UnrestrictedResidualTreeBuilder<std::uint8_t>, RegularGridAdjacency2D, sdrt::SpatialOrder>,
              "unrestricted builders must require their mode-specific option object");
static_assert(!std::is_constructible_v<sdrt::SaturatedResidualTreeBuilder<std::uint8_t>, RegularGridAdjacency2D, PixelId, sdrt::SpatialOrder>,
              "saturated builders must require their mode-specific option object");

template <class T> std::vector<T> toVector(std::span<const T> values) { return {values.begin(), values.end()}; }

template <class T> std::vector<NodeId> parentBuffer(const ValuedMorphologicalTree<T>& tree) {
    std::vector<NodeId> parent(static_cast<std::size_t>(tree.topology().numInternalNodeSlots()));
    for (NodeId node : tree.topology().aliveNodeIds()) {
        parent[static_cast<std::size_t>(node)] = tree.topology().parent(node);
    }
    return parent;
}

template <class T> std::vector<NodeId> smallestNodeMapBuffer(const ValuedMorphologicalTree<T>& tree) {
    std::vector<NodeId> smallestNodeMap(static_cast<std::size_t>(tree.topology().numPixels()));
    for (PixelId pixel = 0; pixel < tree.topology().numPixels(); ++pixel) {
        smallestNodeMap[static_cast<std::size_t>(pixel)] = tree.topology().smallestNode(pixel);
    }
    return smallestNodeMap;
}

template <class T> std::vector<T> altitudeBuffer(const ValuedMorphologicalTree<T>& tree) {
    std::vector<T> altitude(static_cast<std::size_t>(tree.topology().numInternalNodeSlots()));
    for (NodeId node : tree.topology().aliveNodeIds()) {
        altitude[static_cast<std::size_t>(node)] = tree.nodeAltitude(node);
    }
    return altitude;
}

template <class T> struct ResidualBuffers {
    std::vector<NodeId> parent;
    std::vector<NodeId> smallestNodeMap;
    std::vector<T> altitude;
};

template <class T>
ResidualBuffers<T> buildResidualWithPolicies(const ImagePtr<T>& image, const RegularGridAdjacency2D& adjacency, PixelId infinityPixel,
                                             sdrt::SaturatedMinMaxLcaPolicy lcaPolicy, sdrt::SaturatedMinMaxFallbackPolicy fallbackPolicy) {
    auto minTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto maxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
    sdrt::SaturatedResidualTreeOptions options;
    options.lcaPolicy = lcaPolicy;
    options.fallbackPolicy = fallbackPolicy;
    sdrt::SaturatedResidualTreeBuilder<T> builder(adjacency, infinityPixel, options);
    builder.build(image, std::move(minTree), std::move(maxTree));
    return {toVector(builder.parents()), toVector(builder.smallestNodeMap()), toVector(builder.nodeAltitudes())};
}

ImageUInt8Ptr makeRadixImage(int code) {
    auto image = ImageUInt8::create(2, 3);
    for (PixelId pixel = 0; pixel < 6; ++pixel) {
        (*image)[pixel] = static_cast<std::uint8_t>(code % 3);
        code /= 3;
    }
    return image;
}

template <class T> ImagePtr<T> makeTypedImage(int rows, int columns, std::initializer_list<T> values) {
    requireEqual(static_cast<int>(values.size()), rows * columns, "typed image buffer size");
    auto image = Image<T>::create(rows, columns);
    PixelId pixel = 0;
    for (T value : values) {
        (*image)[pixel++] = value;
    }
    return image;
}

template <class T> void requireExactReconstruction(const ValuedMorphologicalTree<T>& tree, const ImagePtr<T>& image, const std::string& label) {
    for (PixelId pixel = 0; pixel < image->getSize(); ++pixel) {
        const NodeId smallestNode = tree.topology().smallestNode(pixel);
        requireEqual(tree.nodeAltitude(smallestNode), (*image)[pixel], label);
    }
}

void requireResidualSemantics(const MorphologicalTree& tree, const RegularGridAdjacency2D& adjacency, MorphologicalTreeKind expectedKind,
                              const std::string& label) {
    require(tree.kind() == expectedKind, label + ": descriptive kind");
    require(tree.nodeAltitudeOrder() == NodeAltitudeOrder::Unconstrained, label + ": unconstrained altitude order");
    const RegularGridAdjacency2D* stored = nullptr;
    if (const auto* context = tree.sharedAdjacencyContext()) {
        stored = &context->adjacency;
    } else if (const auto* context = tree.saturatedResidualContext()) {
        stored = &context->adjacency;
    }
    require(stored != nullptr, label + ": retained construction adjacency");
    requireEqual(stored->getNumRows(), adjacency.getNumRows(), label + ": adjacency rows");
    requireEqual(stored->getNumColumns(), adjacency.getNumColumns(), label + ": adjacency columns");
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

    const auto unrestricted = MorphologicalTreeFactory::createUnrestrictedResidualTree(image, adjacency);
    const auto saturated = MorphologicalTreeFactory::createSaturatedResidualTree(image, adjacency, PixelId{0});

    requireExactReconstruction(unrestricted, image, "unrestricted residual reconstruction");
    requireExactReconstruction(saturated, image, "saturated residual reconstruction");
    requireResidualSemantics(unrestricted.topology(), adjacency, MorphologicalTreeKind::UnrestrictedResidualTree, "unrestricted residual semantics");
    requireResidualSemantics(saturated.topology(), adjacency, MorphologicalTreeKind::SaturatedResidualTree, "saturated residual semantics");

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

    const auto requireConjugateEvolution = [&](const auto& primary, const auto& conjugate, const std::string& label) {
        const MorphologicalTree& lhs = primary.topology();
        const MorphologicalTree& rhs = conjugate.topology();
        requireEqual(lhs.numNodes(), rhs.numNodes(), label + " node count");
        std::size_t previousSupportCardinality = 0;
        for (NodeId node : lhs.aliveNodeIds()) {
            requireEqual(lhs.parent(node), rhs.parent(node), label + " parent");
            requireEqual(static_cast<int>(primary.nodeAltitude(node)) + static_cast<int>(conjugate.nodeAltitude(node)), 3,
                         label + " complemented altitude");
            const auto lhsSupportRange = lhs.nodeSupport(node);
            const auto rhsSupportRange = rhs.nodeSupport(node);
            const std::vector<NodeId> lhsSupport(lhsSupportRange.begin(), lhsSupportRange.end());
            const std::vector<NodeId> rhsSupport(rhsSupportRange.begin(), rhsSupportRange.end());
            requireVectorEqual(lhsSupport, rhsSupport, label + " event support");
            if (node != lhs.root()) {
                requireEqual(primary.nodeResidue(node), -conjugate.nodeResidue(node), label + " signed residue");
                const auto supportCardinality = lhsSupport.size();
                require(supportCardinality >= previousSupportCardinality, label + " nondecreasing support cardinality");
                previousSupportCardinality = supportCardinality;
            }
        }
        for (PixelId pixel = 0; pixel < image->getSize(); ++pixel) {
            requireEqual(lhs.smallestNode(pixel), rhs.smallestNode(pixel), label + " smallest node");
        }
    };

    requireConjugateEvolution(MorphologicalTreeFactory::createUnrestrictedResidualTree(image, adjacency),
                              MorphologicalTreeFactory::createUnrestrictedResidualTree(inverted, adjacency), "unrestricted contrast");
    requireConjugateEvolution(MorphologicalTreeFactory::createSaturatedResidualTree(image, adjacency, PixelId{0}),
                              MorphologicalTreeFactory::createSaturatedResidualTree(inverted, adjacency, PixelId{0}),
                              "saturated contrast");

    sdrt::UnrestrictedResidualTreeOptions customOrderOptions;
    customOrderOptions.spatialOrder = sdrt::SpatialOrder(std::vector<PixelId>{5, 4, 3, 2, 1, 0});
    requireConjugateEvolution(MorphologicalTreeFactory::createUnrestrictedResidualTree(image, adjacency, customOrderOptions),
                              MorphologicalTreeFactory::createUnrestrictedResidualTree(inverted, adjacency, customOrderOptions),
                              "custom-order unrestricted contrast");
}

void testSelfDualResidualSchedule() {
    const sdrt::SelfDualResidualOrder rowMajorOrder;
    const sdrt::SelfDualResidualKey first{2, 0};
    const sdrt::SelfDualResidualKey second{2, 3};
    require(rowMajorOrder.compareResidualCandidates(first, second), "row-major residual order uses the spatial minimum");
    require(!rowMajorOrder.compareResidualCandidates(first, first), "equal residual keys are equivalent");

    const sdrt::SpatialOrder reverseOrder(std::vector<PixelId>{3, 2, 1, 0});
    const sdrt::SelfDualResidualSchedule schedule(reverseOrder);
    const std::array<sdrt::SelfDualResidualKey, 3> keys{{{2, 1}, {1, 2}, {1, 3}}};
    requireEqual(schedule.selectResidualCandidate(keys), std::size_t{2}, "custom spatial order selects the least equal-cardinality support");
    const std::array<sdrt::SelfDualResidualKey, 2> duplicateKeys{{{1, 2}, {1, 2}}};
    requireThrows<std::invalid_argument>([&] { static_cast<void>(schedule.selectResidualCandidate(duplicateKeys)); },
                                         "self-dual residual schedule rejects duplicate keys");
    requireEqual(reverseOrder.spatialMinimum(std::array<PixelId, 3>{0, 2, 3}), PixelId{3}, "custom spatial minimum");
    requireThrows<std::invalid_argument>([] { static_cast<void>(sdrt::SpatialOrder(std::vector<PixelId>{0, 0})); },
                                         "spatial order rejects duplicate pixels");
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
    sdrt::UnrestrictedResidualTreeBuilder<std::uint8_t> builder(adjacency);
    builder.build(image, std::move(minTree), std::move(maxTree));
    requireEqual(builder.statistics().rejectedExtrema, std::size_t{0}, "unrestricted builder rejects no extrema");
    requireEqual(builder.statistics().complementTraversalCertificates, std::size_t{0}, "unrestricted builder performs no complement traversal");
}

void testResidualTreeEventAssemblerInIsolation() {
    using sdrt::detail::RegionId;
    const std::vector<RegionId> regionByPixel{0, 0, 1, 2};
    sdrt::detail::ResidualTreeEventAssembler<int> assembler(3, regionByPixel);

    const std::array<PixelId, 2> firstSupport{0, 1};
    const sdrt::ResidualCandidate<int> firstCandidate{firstSupport, sdrt::Polarity::Maximum, 5, 3, {2, 0}};
    const auto firstEvent = sdrt::recordResidualEvent(std::size_t{0}, firstCandidate);
    requireEqual(firstEvent.signedResidualValue, std::int64_t{2}, "assembler first signed residual value");
    requireEqual(assembler.emitEvent(RegionId{0}, firstEvent), NodeId{1}, "assembler first event id");
    const std::array<RegionId, 1> firstAbsorption{0};
    assembler.consume(RegionId{1}, firstAbsorption);
    const std::array<PixelId, 3> secondSupport{0, 1, 2};
    const sdrt::ResidualCandidate<int> secondCandidate{secondSupport, sdrt::Polarity::Maximum, 3, 1, {3, 0}};
    const auto secondEvent = sdrt::recordResidualEvent(std::size_t{1}, secondCandidate);
    requireEqual(assembler.emitEvent(RegionId{1}, secondEvent), NodeId{2}, "assembler second event id");
    const std::array<RegionId, 1> secondAbsorption{1};
    assembler.consume(RegionId{2}, secondAbsorption);
    requireEqual(assembler.numEvents(), std::size_t{2}, "assembler event count");

    auto output = std::move(assembler).finalize(RegionId{2}, 1);
    requireVectorEqual(output.parents, std::vector<NodeId>{0, 2, 0}, "assembler parent materialization");
    requireVectorEqual(output.smallestNodeMap, std::vector<NodeId>{1, 1, 2, 0}, "assembler smallest-node-map materialization");
    requireVectorEqual(output.nodeAltitudes, std::vector<int>{1, 5, 3}, "assembler altitude materialization");

    sdrt::detail::ResidualTreeEventAssembler<int> invalidAssembler(2, std::vector<RegionId>{0, 1});
    const std::array<RegionId, 1> absorbsSelf{0};
    requireThrows<std::invalid_argument>([&] { invalidAssembler.consume(RegionId{0}, absorbsSelf); }, "assembler rejects absorbing its survivor");

    const std::array<PixelId, 1> unsignedSupport{0};
    const sdrt::ResidualCandidate<std::uint16_t> minimumCandidate{
        unsignedSupport, sdrt::Polarity::Minimum, std::uint16_t{1}, std::uint16_t{7}, {1, 0}};
    const auto minimumEvent = sdrt::recordResidualEvent(std::size_t{0}, minimumCandidate);
    requireEqual(minimumEvent.signedResidualValue, std::int64_t{-6}, "unsigned altitude produces a signed negative residual value");
}

void testIncrementalFlatZonePartition() {
    const auto chainImage = makeImage(1, 4, {2, 1, 1, 0});
    const RegularGridAdjacency2D chainAdjacency(1, 4, 1.0);
    sdrt::detail::FlatZonePartition<std::uint8_t> partition;
    partition.initialize(chainImage, chainAdjacency, PixelId{3});
    GenerationStampSet representativeMarks(static_cast<std::size_t>(chainImage->getSize()));
    std::vector<PixelId> mergeRepresentatives;

    const auto requireViewMatchesSupport = [&](const auto& view, std::vector<PixelId> expectedSupport, const std::string& label) {
        std::vector<std::uint8_t> supportMarks(static_cast<std::size_t>(chainImage->getSize()), std::uint8_t{0});
        for (PixelId pixel : expectedSupport) {
            supportMarks[static_cast<std::size_t>(pixel)] = 1;
        }
        std::vector<PixelId> expectedBoundary;
        for (PixelId pixel : expectedSupport) {
            for (PixelId neighbor : chainAdjacency.getNeighborIndices(pixel)) {
                if (supportMarks[static_cast<std::size_t>(neighbor)] == 0) {
                    expectedBoundary.push_back(neighbor);
                }
            }
        }
        std::vector<PixelId> actualSupport(view.supportPixels.begin(), view.supportPixels.end());
        std::vector<PixelId> actualBoundary(view.externalPixelsByIncidence.begin(), view.externalPixelsByIncidence.end());
        std::ranges::sort(expectedSupport);
        std::ranges::sort(expectedBoundary);
        std::ranges::sort(actualSupport);
        std::ranges::sort(actualBoundary);
        requireVectorEqual(actualSupport, expectedSupport, label + " support");
        requireVectorEqual(actualBoundary, expectedBoundary, label + " boundary");
    };

    auto selectedView = partition.viewForPixel(PixelId{0});
    requireViewMatchesSupport(selectedView, {0}, "initial selected flat zone");
    partition.collectAdjacentRepresentativesAtLevel(selectedView.representative, std::uint8_t{1}, selectedView.externalPixelsByIncidence,
                                                    representativeMarks, mergeRepresentatives);
    requireEqual(mergeRepresentatives.size(), std::size_t{1}, "first incremental flat-zone merge count");
    NodeId mergedRepresentative = partition.mergeFlatZonesAtLevel(selectedView.representative, mergeRepresentatives, std::uint8_t{1});
    auto mergedView = partition.viewForPixel(mergedRepresentative);
    requireViewMatchesSupport(mergedView, {0, 1, 2}, "first merged flat zone");
    requireEqual(mergedView.spatialMinimum, PixelId{0}, "merged flat-zone spatial minimum");
    require(!mergedView.containsInfinityPixel, "first merged flat zone excludes infinity pixel");

    partition.collectAdjacentRepresentativesAtLevel(mergedView.representative, std::uint8_t{0}, mergedView.externalPixelsByIncidence,
                                                    representativeMarks, mergeRepresentatives);
    requireEqual(mergeRepresentatives.size(), std::size_t{1}, "second incremental flat-zone merge count");
    mergedRepresentative = partition.mergeFlatZonesAtLevel(mergedView.representative, mergeRepresentatives, std::uint8_t{0});
    const auto terminalView = partition.viewForPixel(mergedRepresentative);
    requireViewMatchesSupport(terminalView, {0, 1, 2, 3}, "terminal merged flat zone");
    require(terminalView.containsInfinityPixel, "terminal merged flat zone retains infinity pixel");
    requireEqual(partition.representativeOf(NodeId{3}), mergedRepresentative, "terminal flat-zone representative");
}

template <class T> void testTypedResidualConstruction(const std::string& label) {
    const auto image = makeTypedImage<T>(2, 3, {T{-3}, T{2}, T{1}, T{4}, T{-1}, T{0}});
    const RegularGridAdjacency2D adjacency(2, 3, 1.0);
    const auto unrestricted = MorphologicalTreeFactory::createUnrestrictedResidualTree(image, adjacency);
    const auto saturated = MorphologicalTreeFactory::createSaturatedResidualTree(image, adjacency, PixelId{0});
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
    requireExactReconstruction(MorphologicalTreeFactory::createUnrestrictedResidualTree(unsignedImage, unsignedAdjacency), unsignedImage,
                               "uint16 unrestricted reconstruction");
    requireExactReconstruction(MorphologicalTreeFactory::createSaturatedResidualTree(unsignedImage, unsignedAdjacency, PixelId{0}), unsignedImage,
                               "uint16 saturated reconstruction");

    const auto singleton = makeTypedImage<double>(1, 1, {2.5});
    const auto unrestrictedSingleton = MorphologicalTreeFactory::createUnrestrictedResidualTree(singleton, 1.0);
    const auto saturatedSingleton = MorphologicalTreeFactory::createSaturatedResidualTree(singleton, PixelId{0}, 1.0);
    requireEqual(unrestrictedSingleton.topology().numNodes(), 1, "unrestricted singleton node count");
    requireEqual(saturatedSingleton.topology().numNodes(), 1, "saturated singleton node count");
    requireExactReconstruction(unrestrictedSingleton, singleton, "unrestricted singleton reconstruction");
    requireExactReconstruction(saturatedSingleton, singleton, "saturated singleton reconstruction");
}

void testOptionsAndInputContracts() {
    const auto image = makeImage(2, 3, {0, 2, 1, 2, 0, 0});
    const RegularGridAdjacency2D adjacency(2, 3, 1.0);

    sdrt::UnrestrictedResidualTreeOptions unrestrictedOptions;
    unrestrictedOptions.spatialOrder = sdrt::SpatialOrder(std::vector<PixelId>{5, 4, 3, 2, 1, 0});
    const auto unrestrictedWithOptions = MorphologicalTreeFactory::createUnrestrictedResidualTree(image, adjacency, unrestrictedOptions);
    auto optionMinTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto optionMaxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
    sdrt::UnrestrictedResidualTreeBuilder<std::uint8_t> unrestrictedBuilder(adjacency, unrestrictedOptions);
    unrestrictedBuilder.build(image, std::move(optionMinTree), std::move(optionMaxTree));
    require(!unrestrictedBuilder.spatialOrder().isRowMajor(), "unrestricted builder retains its explicit spatial order");
    requireVectorEqual(parentBuffer(unrestrictedWithOptions), toVector(unrestrictedBuilder.parents()), "unrestricted option parents");
    requireVectorEqual(smallestNodeMapBuffer(unrestrictedWithOptions), toVector(unrestrictedBuilder.smallestNodeMap()),
                       "unrestricted option smallest-node map");
    requireVectorEqual(altitudeBuffer(unrestrictedWithOptions), toVector(unrestrictedBuilder.nodeAltitudes()), "unrestricted option altitudes");

    sdrt::SaturatedResidualTreeOptions saturatedOptions;
    saturatedOptions.lcaPolicy = sdrt::SaturatedMinMaxLcaPolicy::LinkCut;
    saturatedOptions.fallbackPolicy = sdrt::SaturatedMinMaxFallbackPolicy::SingleSourceDepthFirst;
    const auto saturatedWithOptions = MorphologicalTreeFactory::createSaturatedResidualTree(image, adjacency, PixelId{0}, saturatedOptions);
    const auto saturatedDefault = MorphologicalTreeFactory::createSaturatedResidualTree(image, adjacency, PixelId{0});
    requireVectorEqual(parentBuffer(saturatedWithOptions), parentBuffer(saturatedDefault), "saturated option parents");
    requireVectorEqual(smallestNodeMapBuffer(saturatedWithOptions), smallestNodeMapBuffer(saturatedDefault),
                       "saturated option smallest-node maps");
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
    requireThrows<std::logic_error>([&] { static_cast<void>(reusable.parents()); }, "failed construction leaves no observable partial result");
    auto validMinTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto validMaxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
    reusable.build(image, std::move(validMinTree), std::move(validMaxTree));
    requireEqual(reusable.smallestNodeMap().size(), static_cast<std::size_t>(image->getSize()), "builder reuse after failure");

    const std::array<GridOffset2D, 1> originOnly{{{0, 0}}};
    const auto disconnected = RegularGridAdjacency2D::fromStructuringElement(2, 3, originOnly);
    sdrt::UnrestrictedResidualTreeBuilder<std::uint8_t> disconnectedBuilder(disconnected);
    auto connectedMinTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto connectedMaxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
    requireThrows<std::invalid_argument>([&] { disconnectedBuilder.build(image, std::move(connectedMinTree), std::move(connectedMaxTree)); },
                                         "residual builder rejects a disconnected adjacency");

    if constexpr (contract::validationsEnabled) {
        auto nonFinite = makeTypedImage<float>(1, 2, {0.0F, std::numeric_limits<float>::quiet_NaN()});
        requireThrows<std::invalid_argument>([&] { static_cast<void>(MorphologicalTreeFactory::createUnrestrictedResidualTree(nonFinite, 1.0)); },
                                             "unrestricted factory rejects NaN");
        (*nonFinite)[1] = std::numeric_limits<float>::infinity();
        requireThrows<std::invalid_argument>(
            [&] { static_cast<void>(MorphologicalTreeFactory::createSaturatedResidualTree(nonFinite, PixelId{0}, 1.0)); },
            "saturated factory rejects infinity");
    }
}

void testSharedArbitrarySymmetricAdjacency() {
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
    const auto unrestricted = MorphologicalTreeFactory::createUnrestrictedResidualTree(image, adjacency);
    const auto saturated = MorphologicalTreeFactory::createSaturatedResidualTree(image, adjacency, PixelId{0});

    requireExactReconstruction(unrestricted, image, "custom-adjacency unrestricted reconstruction");
    requireExactReconstruction(saturated, image, "custom-adjacency saturated reconstruction");
    requireResidualSemantics(unrestricted.topology(), adjacency, MorphologicalTreeKind::UnrestrictedResidualTree, "custom-adjacency unrestricted semantics");
    requireResidualSemantics(saturated.topology(), adjacency, MorphologicalTreeKind::SaturatedResidualTree, "custom-adjacency saturated semantics");
}

bool evaluateMaximumCandidateEligibility(const ImageUInt8Ptr& image, PixelId candidatePixel, PixelId infinityPixel) {
    const RegularGridAdjacency2D adjacency(image->getNumRows(), image->getNumColumns(), 1.0);
    auto minTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto maxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
    const NodeId candidateNode = maxTree.topology().smallestNode(candidatePixel);
    require(candidateNode != maxTree.topology().root(), "isolated saturated-eligibility candidate must be a non-root maximum");

    sdrt::detail::FlatZonePartition<std::uint8_t> partition;
    partition.initialize(image, adjacency, infinityPixel);
    const auto maxNodeSlots = static_cast<std::size_t>(
        std::max(minTree.topology().numInternalNodeSlots(), maxTree.topology().numInternalNodeSlots()));
    sdrt::detail::ResidualTreeCandidateContext context(static_cast<std::size_t>(image->getSize()), maxNodeSlots);
    const bool containsInfinityPixel =
        sdrt::detail::prepareResidualTreeCandidate(partition, context, candidateNode, maxTree, minTree);

    sdrt::detail::SaturatedResidualEligibility<std::uint8_t> eligibility(
        static_cast<std::size_t>(image->getSize()), adjacency, infinityPixel, sdrt::SaturatedMinMaxLcaPolicy::ParentClimb,
        sdrt::SaturatedMinMaxFallbackPolicy::BoundaryMultiSource, minTree, maxTree);
    sdrt::ResidualTreeBuildStatistics statistics;
    return eligibility.isEligible(maxTree, minTree, candidateNode, true, containsInfinityPixel, context, statistics);
}

void testSaturatedResidualEligibilityInIsolation() {
    require(evaluateMaximumCandidateEligibility(makeImage(1, 3, {0, 0, 1}), PixelId{2}, PixelId{0}),
            "an end maximum has connected complement");
    require(!evaluateMaximumCandidateEligibility(makeImage(1, 3, {0, 1, 0}), PixelId{1}, PixelId{0}),
            "an interior maximum disconnects its complement");
    require(!evaluateMaximumCandidateEligibility(makeImage(1, 3, {0, 0, 1}), PixelId{2}, PixelId{2}),
            "a candidate containing the infinity pixel is ineligible");
}

void testConfigurablePolicyEquivalence() {
    constexpr std::array lcaPolicies{sdrt::SaturatedMinMaxLcaPolicy::ParentClimb, sdrt::SaturatedMinMaxLcaPolicy::BlockedSnapshot,
                                     sdrt::SaturatedMinMaxLcaPolicy::LinkCut};
    constexpr std::array fallbackPolicies{sdrt::SaturatedMinMaxFallbackPolicy::SingleSourceDepthFirst,
                                          sdrt::SaturatedMinMaxFallbackPolicy::BoundaryMultiSource};
    const RegularGridAdjacency2D adjacency(2, 3, 1.0);

    for (int code = 0; code < 729; ++code) {
        const auto image = makeRadixImage(code);
        const PixelId infinityPixel = static_cast<PixelId>(code % 6);
        const auto expected = MorphologicalTreeFactory::createSaturatedResidualTree(image, adjacency, infinityPixel);
        const auto expectedParent = parentBuffer(expected);
        const auto expectedSmallestNodeMap = smallestNodeMapBuffer(expected);
        const auto expectedAltitude = altitudeBuffer(expected);

        for (const auto lcaPolicy : lcaPolicies) {
            for (const auto fallbackPolicy : fallbackPolicies) {
                const auto actual = buildResidualWithPolicies(image, adjacency, infinityPixel, lcaPolicy, fallbackPolicy);
                requireVectorEqual(actual.parent, expectedParent, "configurable policy parents");
                requireVectorEqual(actual.smallestNodeMap, expectedSmallestNodeMap, "configurable policy smallest-node maps");
                requireVectorEqual(actual.altitude, expectedAltitude, "configurable policy altitudes");
            }
        }
    }
}

} // namespace

int main() {
    testFactoryModesAndSaliencyIntegration();
    testContrastInversion();
    testSelfDualResidualSchedule();
    testDirectBuilderPolicies();
    testResidualTreeEventAssemblerInIsolation();
    testIncrementalFlatZonePartition();
    testAltitudeTypesAndDegenerateDomain();
    testOptionsAndInputContracts();
    testSharedArbitrarySymmetricAdjacency();
    testSaturatedResidualEligibilityInIsolation();
    testConfigurablePolicyEquivalence();
    return 0;
}
