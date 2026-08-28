#include "support/TestSupport.hpp"

#include <algorithm>
#include <memory>
#include <type_traits>
#include <vector>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

namespace {

std::vector<std::vector<NodeId>> supportFamily(const MorphologicalTree& tree) {
    std::vector<std::vector<NodeId>> supports;
    supports.reserve(static_cast<std::size_t>(tree.numNodes()));
    for (NodeId nodeId : tree.aliveNodeIds()) {
        std::vector<NodeId> support = collectNodeIds(tree.nodeSupport(nodeId));
        std::sort(support.begin(), support.end());
        supports.push_back(std::move(support));
    }
    std::sort(supports.begin(), supports.end());
    return supports;
}

ImageUInt8Ptr makeImageFromVector(int rows, int columns, const std::vector<uint8_t>& values) {
    requireEqual(static_cast<int>(values.size()), rows * columns, "vector image buffer size");
    auto image = ImageUInt8::create(rows, columns);
    std::copy(values.begin(), values.end(), image->rawData());
    return image;
}

} // namespace

int main() {
    auto image = makeImage(3, 3,
                           {
                               1,
                               2,
                               1,
                               2,
                               3,
                               2,
                               1,
                               2,
                               1,
                           });

    auto tree = makeTreeOfShapes(image, TestTopographicImmersion::Min4Max8);
    auto valuedTree = makeValuedTreeOfShapes(image, TestTopographicImmersion::Min4Max8);

    require(static_cast<bool>(tree), "tree of shapes instance must be created");
    require(tree->kind() == MorphologicalTreeKind::TreeOfShapes, "tree type");
    require(tree->numNodes() > 0, "tree of shapes must contain nodes");
    require(tree->root() != InvalidNode, "tree of shapes root must be valid");
    require(tree->sharedAdjacencyContext() == nullptr, "tree of shapes must not expose a shared-adjacency context");
    const TopographicConvention* convention = tree->topographicConvention();
    require(convention != nullptr, "tree of shapes must retain its topographic convention");
    const auto* complementary = std::get_if<ComplementaryGridImmersion>(&convention->immersion);
    require(complementary != nullptr, "Min4cMax8c must retain a complementary-grid immersion");
    require(tree->nodeAltitudeOrder() == NodeAltitudeOrder::Unconstrained, "tree of shapes altitude order");
    requireNear(complementary->complementaryAdjacencies.minAdjacency.getRadius(), 1.0, 0.0, "Min4cMax8c minimum radius");
    requireNear(complementary->complementaryAdjacencies.maxAdjacency.getRadius(), 1.5, 0.0, "Min4cMax8c maximum radius");
    requireEqual(complementary->complementaryAdjacencies.minAdjacency.getSize(), 5, "Min4cMax8c minimum adjacency size");
    requireEqual(complementary->complementaryAdjacencies.maxAdjacency.getSize(), 9, "Min4cMax8c maximum adjacency size");
    requireEqual(tree->numRows(), 3, "tree of shapes image rows");
    requireEqual(tree->numColumns(), 3, "tree of shapes image columns");
    for (NodeId nodeId : tree->aliveNodeIds()) {
        require(!collectNodeIds(tree->nodeSupport(nodeId)).empty(), "every committed tree-of-shapes node must have non-empty subtree support");
    }

    auto dualTree = makeTreeOfShapes(image, TestTopographicImmersion::SelfDualSpan);
    require(dualTree->kind() == MorphologicalTreeKind::TreeOfShapes, "self-dual construction must retain tree-of-shapes kind");
    require(dualTree->topographicConvention() != nullptr, "self-dual tree of shapes must retain its topographic convention");
    require(std::holds_alternative<SelfDualSpanImmersion>(dualTree->topographicConvention()->immersion),
            "self-dual tree of shapes must retain the span immersion");

    auto min8Max4Tree = makeTreeOfShapes(image, TestTopographicImmersion::Min8Max4);
    require(min8Max4Tree->kind() == MorphologicalTreeKind::TreeOfShapes, "inverse complementary construction must retain tree-of-shapes kind");
    const auto* inverseComplementary =
        std::get_if<ComplementaryGridImmersion>(&min8Max4Tree->topographicConvention()->immersion);
    require(inverseComplementary != nullptr, "Min8cMax4c must retain a complementary-grid immersion");
    requireNear(inverseComplementary->complementaryAdjacencies.minAdjacency.getRadius(), 1.5, 0.0, "Min8cMax4c minimum radius");
    requireNear(inverseComplementary->complementaryAdjacencies.maxAdjacency.getRadius(), 1.0, 0.0, "Min8cMax4c maximum radius");
    requireEqual(inverseComplementary->complementaryAdjacencies.minAdjacency.getSize(), 9, "Min8cMax4c minimum adjacency size");
    requireEqual(inverseComplementary->complementaryAdjacencies.maxAdjacency.getSize(), 5, "Min8cMax4c maximum adjacency size");

    auto reconstruction = valuedTree->reconstructFromNodeAltitudes();
    requireImageShape(reconstruction, 3, 3);
    requireVectorEqual(collectImageValues(reconstruction), std::vector<ToSGrayLevel>{2, 4, 2, 4, 6, 4, 2, 4, 2},
                       "tree of shapes exact doubled-unit reconstruction values");

    auto singlePixel = makeImage(1, 1, {5});
    auto singleTree = makeTreeOfShapes(singlePixel);
    auto singleValuedTree = makeValuedTreeOfShapes(singlePixel);
    require(static_cast<bool>(singleTree), "single-pixel self-dual tree of shapes instance must be created");
    requireEqual(singleTree->numRows(), 1, "single-pixel ToS rows");
    requireEqual(singleTree->numColumns(), 1, "single-pixel ToS columns");
    require(singleTree->root() != InvalidNode, "single-pixel ToS root must be valid");
    requireEqual(singleTree->smallestNode(0), singleTree->root(), "single-pixel ToS smallest node");
    auto singleReconstruction = singleValuedTree->reconstructFromNodeAltitudes();
    requireImageShape(singleReconstruction, 1, 1);
    requireVectorEqual(collectImageValues(singleReconstruction), std::vector<ToSGrayLevel>{10}, "single-pixel exact ToS reconstruction");
    for (TestTopographicImmersion interpolation : {TestTopographicImmersion::SelfDualSpan, TestTopographicImmersion::Min4Max8, TestTopographicImmersion::Min8Max4}) {
        auto unpaddedSingle = MorphologicalTreeFactory::createTreeOfShapes<ToSGrayLevel>(
            singlePixel, makeTopographicConvention(singlePixel, interpolation, TopographicDomainExtension::None));
        requireEqual(unpaddedSingle.topology().numNodes(), 1, "single-pixel unpadded ToS node count");
        requireVectorEqual(collectImageValues(unpaddedSingle.reconstructFromNodeAltitudes()), std::vector<ToSGrayLevel>{10},
                           "single-pixel unpadded exact ToS reconstruction");
    }

    {
        auto virtualRootImage = makeImage(2, 2,
                                          {
                                              1,
                                              1,
                                              0,
                                              0,
                                          });
        auto virtualRootValuedTree = makeValuedTreeOfShapes(virtualRootImage, TestTopographicImmersion::SelfDualSpan);
        const MorphologicalTree& virtualRootTree = virtualRootValuedTree->topology();
        const NodeId virtualRoot = virtualRootTree.root();
        const NodeId upperShape = virtualRootTree.smallestNode(0);
        const NodeId lowerShape = virtualRootTree.smallestNode(2);

        requireEqual(virtualRootTree.numNodes(), 3, "virtual-root ToS node count");
        requireEqual(virtualRootTree.properPartCardinality(virtualRoot), 0, "virtual-root ToS root direct proper part");
        require(virtualRootTree.hasEmptyProperPart(virtualRoot), "virtual-root ToS root must be derived as structural");
        requireEqual(virtualRootTree.numChildren(virtualRoot), 2, "virtual-root ToS root children");
        requireEqual(computeAreaViaAttributeFacade(virtualRootTree, virtualRoot), 4, "virtual-root ToS cumulative area");
        require(upperShape != lowerShape, "virtual-root ToS must preserve both child shapes");
        requireEqual(computeAreaViaAttributeFacade(virtualRootTree, upperShape), 2, "upper shape area");
        requireEqual(computeAreaViaAttributeFacade(virtualRootTree, lowerShape), 2, "lower shape area");
        requireEqual(virtualRootTree.parent(upperShape), virtualRoot, "upper shape parent");
        requireEqual(virtualRootTree.parent(lowerShape), virtualRoot, "lower shape parent");
        requireVectorEqual(collectImageValues(virtualRootValuedTree->reconstructFromNodeAltitudes()), std::vector<ToSGrayLevel>{2, 2, 0, 0},
                           "virtual-root exact ToS reconstruction");
        static_assert(std::is_same_v<std::remove_cvref_t<decltype(*virtualRootValuedTree)>, ValuedMorphologicalTree<ToSGrayLevel>>);
        requireEqual(virtualRootValuedTree->nodeAltitude(virtualRoot), static_cast<ToSGrayLevel>(1), "tree of shapes preserves a half-level virtual root");
        const NodeId lowerOwner = virtualRootTree.smallestNode(2);
        require(lowerOwner != virtualRoot, "exact ToS lower shape must remain below its virtual root");
        requireEqual(virtualRootValuedTree->nodeAltitude(lowerOwner), static_cast<ToSGrayLevel>(0), "exact ToS lower source altitude");
        require(virtualRootValuedTree->nodeAltitude(lowerOwner) != virtualRootValuedTree->nodeAltitude(virtualRoot),
                "exact ToS must preserve the half-level parent-child distinction");

        const auto [exportedParent, exportedAltitude] = virtualRootValuedTree->exportHigraHierarchy();
        requireEqual(static_cast<int>(exportedParent.size()), 7, "virtual-root ToS Higra export vertex count");
        requireEqual(exportedParent.size(), exportedAltitude.size(), "virtual-root ToS Higra export altitude size");

        auto invertedImage = makeImage(2, 2,
                                       {
                                           1,
                                           1,
                                           2,
                                           2,
                                       });
        auto invertedValuedTree = makeValuedTreeOfShapes(invertedImage, TestTopographicImmersion::SelfDualSpan);
        require(supportFamily(virtualRootTree) == supportFamily(invertedValuedTree->topology()),
                "SelfDual ToS support family must be invariant under contrast inversion");
    }

    {
        auto unpaddedImage = makeImage(2, 3,
                                       {
                                           0,
                                           2,
                                           1,
                                           2,
                                           1,
                                           0,
                                       });
        for (TestTopographicImmersion interpolation : {TestTopographicImmersion::SelfDualSpan, TestTopographicImmersion::Min4Max8, TestTopographicImmersion::Min8Max4}) {
            const TopographicConvention convention =
                makeTopographicConvention(unpaddedImage, interpolation, TopographicDomainExtension::None);
            auto unpadded = MorphologicalTreeFactory::createTreeOfShapes<ToSGrayLevel>(unpaddedImage, convention);
            requireEqual(unpadded.topology().numRows(), 2, "unpadded ToS published rows");
            requireEqual(unpadded.topology().numColumns(), 3, "unpadded ToS published columns");
            requireEqual(unpadded.topology().numPixels(), 6, "unpadded ToS pixel domain size");
            requireVectorEqual(collectImageValues(unpadded.reconstructFromNodeAltitudes()), std::vector<ToSGrayLevel>{0, 4, 2, 4, 2, 0},
                               "exact unpadded ToS reconstruction");

            TreeOfShapesProducer producer(convention);
            const TreeOfShapesBuildResult<ToSGrayLevel> result = producer.build<ToSGrayLevel>(unpaddedImage);
            requireEqual(result.smallestNodeMap.size(), std::size_t{6}, "unpadded producer publishes only source-domain smallest nodes");
            requireEqual(result.numRows, 2, "unpadded producer source rows");
            requireEqual(result.numColumns, 3, "unpadded producer source columns");
        }
    }

    // Exhaustive ternary small-domain properties include half-level exterior
    // medians and branching shapes that need virtual projected nodes.
    for (int columns : {2, 3}) {
        const int numPixels = 2 * columns;
        int numConfigurations = 1;
        for (PixelId pixelId = 0; pixelId < numPixels; ++pixelId) {
            numConfigurations *= 3;
        }

        for (int code = 0; code < numConfigurations; ++code) {
            int digits = code;
            std::vector<uint8_t> values(static_cast<std::size_t>(numPixels), uint8_t{0});
            std::vector<uint8_t> inverted(static_cast<std::size_t>(numPixels), uint8_t{0});
            for (PixelId pixelId = 0; pixelId < numPixels; ++pixelId) {
                values[static_cast<std::size_t>(pixelId)] = static_cast<uint8_t>(digits % 3);
                inverted[static_cast<std::size_t>(pixelId)] = static_cast<uint8_t>(2 - values[static_cast<std::size_t>(pixelId)]);
                digits /= 3;
            }
            auto originalTree = makeTreeOfShapes(makeImageFromVector(2, columns, values), TestTopographicImmersion::SelfDualSpan);
            auto invertedTree = makeTreeOfShapes(makeImageFromVector(2, columns, inverted), TestTopographicImmersion::SelfDualSpan);
            require(supportFamily(*originalTree) == supportFamily(*invertedTree),
                    "exhaustive ternary 2x" + std::to_string(columns) + " SelfDual support-family invariance, code " + std::to_string(code));

            auto min4Max8Tree = makeTreeOfShapes(makeImageFromVector(2, columns, values), TestTopographicImmersion::Min4Max8);
            auto dualMin8Max4Tree = makeTreeOfShapes(makeImageFromVector(2, columns, inverted), TestTopographicImmersion::Min8Max4);
            require(supportFamily(*min4Max8Tree) == supportFamily(*dualMin8Max4Tree),
                    "exhaustive ternary 2x" + std::to_string(columns) + " dual-connectivity support-family invariance, code " + std::to_string(code));
        }
    }

    // The default convention uses the unpadded canonical min-4/max-8 immersion,
    // infinity pixel zero, and unchanged 8-bit source levels.
    {
        const TopographicConvention defaultSpecification{};
        const auto* canonicalSpecification = std::get_if<CanonicalComplementaryGridImmersion>(&defaultSpecification.immersion);
        require(canonicalSpecification != nullptr, "default topographic convention must select the canonical complementary grid");
        require(canonicalSpecification->pairing == ComplementaryPairing::Min4Max8,
                "default topographic convention must select minimum-4/maximum-8 connectivity");
        require(defaultSpecification.domainExtension == TopographicDomainExtension::None,
                "default topographic convention must not pad the source domain");
        requireEqual(defaultSpecification.infinityPixel, PixelId{0}, "default topographic convention infinity pixel");
        require(defaultSpecification.altitudeEncoding == TopographicAltitudeEncoding::UInt8,
                "default topographic convention must publish uint8 altitudes");

        auto defaultTree = MorphologicalTreeFactory::createTreeOfShapes(image);
        static_assert(std::is_same_v<decltype(defaultTree), ValuedMorphologicalTree<std::uint8_t>>, "default tree of shapes must publish uint8 altitudes");

        const TopographicConvention* defaultConvention = defaultTree.topology().topographicConvention();
        require(defaultConvention != nullptr, "default tree of shapes must retain its topographic convention");
        require(defaultConvention->domainExtension == TopographicDomainExtension::None,
                "default tree of shapes must not pad the source domain");
        requireEqual(defaultConvention->infinityPixel, PixelId{0}, "default tree of shapes infinity pixel");
        require(defaultConvention->altitudeEncoding == TopographicAltitudeEncoding::UInt8, "default tree of shapes must declare the 8-bit encoding");

        // The canonical immersion is resolved, so the retained convention always
        // exposes explicit adjacencies bound to the source domain.
        const auto* defaultComplementary = std::get_if<ComplementaryGridImmersion>(&defaultConvention->immersion);
        require(defaultComplementary != nullptr, "default convention must retain a resolved complementary-grid immersion");
        require(defaultComplementary->complementaryAdjacencies.minAdjacency.is4connectivity(), "default resolved minimum adjacency must be 4-connected");
        require(defaultComplementary->complementaryAdjacencies.maxAdjacency.is8connectivity(), "default resolved maximum adjacency must be 8-connected");
        requireEqual(defaultComplementary->complementaryAdjacencies.minAdjacency.getNumRows(), image->getNumRows(), "resolved adjacency rows");
        requireEqual(defaultComplementary->complementaryAdjacencies.maxAdjacency.getNumColumns(), image->getNumColumns(), "resolved adjacency columns");
    }

    // The 8-bit encoding is exact: it is the doubled hierarchy with every
    // altitude halved, so no topology and no altitude distinction is lost.
    for (ComplementaryPairing pairing : {ComplementaryPairing::Min4Max8, ComplementaryPairing::Min8Max4}) {
        const TestTopographicImmersion immersion =
            pairing == ComplementaryPairing::Min4Max8 ? TestTopographicImmersion::Min4Max8 : TestTopographicImmersion::Min8Max4;
        const std::string label = pairing == ComplementaryPairing::Min4Max8 ? "min4/max8" : "min8/max4";

        auto eightBit = makeValuedTreeOfShapesUInt8(image, pairing);
        auto doubled = makeValuedTreeOfShapes(image, immersion);

        requireEqual(eightBit->topology().numNodes(), doubled->topology().numNodes(), label + " 8-bit node count matches the doubled hierarchy");
        requireEqual(eightBit->topology().root(), doubled->topology().root(), label + " 8-bit root matches the doubled hierarchy");
        for (NodeId nodeId : eightBit->topology().aliveNodeIds()) {
            requireEqual(eightBit->topology().parent(nodeId), doubled->topology().parent(nodeId), label + " 8-bit parent matches the doubled hierarchy");
            requireEqual(TopographicInterpolationScale * static_cast<int>(eightBit->nodeAltitude(nodeId)), static_cast<int>(doubled->nodeAltitude(nodeId)),
                         label + " 8-bit altitude is the doubled altitude halved");
        }

        // Reconstruction now returns the source image itself instead of 2 * image.
        requireVectorEqual(collectImageValues(eightBit->reconstructFromNodeAltitudes()), collectImageValues(image),
                           label + " 8-bit reconstruction reproduces the source image");
    }

    // Without an exterior ring the self-dual span never reaches a half level,
    // because the boundary reference level is cropped away and no interior cell
    // reads it. The 8-bit encoding is then exact there too.
    {
        TopographicConvention unpaddedUInt8{SelfDualSpanImmersion{}, TopographicDomainExtension::None, PixelId{0}, TopographicAltitudeEncoding::UInt8};
        auto eightBitSelfDual = MorphologicalTreeFactory::createTreeOfShapes<std::uint8_t>(image, unpaddedUInt8);
        auto doubledSelfDual = MorphologicalTreeFactory::createTreeOfShapes<ToSGrayLevel>(
            image, makeTopographicConvention(image, TestTopographicImmersion::SelfDualSpan, TopographicDomainExtension::None));

        requireEqual(eightBitSelfDual.topology().numNodes(), doubledSelfDual.topology().numNodes(),
                     "unpadded 8-bit self-dual node count matches the doubled hierarchy");
        requireEqual(eightBitSelfDual.topology().root(), doubledSelfDual.topology().root(),
                     "unpadded 8-bit self-dual root matches the doubled hierarchy");
        for (NodeId nodeId : eightBitSelfDual.topology().aliveNodeIds()) {
            requireEqual(eightBitSelfDual.topology().parent(nodeId), doubledSelfDual.topology().parent(nodeId),
                         "unpadded 8-bit self-dual parent matches the doubled hierarchy");
            requireEqual(TopographicInterpolationScale * static_cast<int>(eightBitSelfDual.nodeAltitude(nodeId)),
                         static_cast<int>(doubledSelfDual.nodeAltitude(nodeId)), "unpadded 8-bit self-dual altitude is the doubled altitude halved");
        }
        requireVectorEqual(collectImageValues(eightBitSelfDual.reconstructFromNodeAltitudes()), collectImageValues(image),
                           "unpadded 8-bit self-dual reconstruction reproduces the source image");
    }

    if constexpr (contract::validationsEnabled) {
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(ImageUInt8::create(0, 0)); }, "empty image creation must throw");

        // The exterior ring carries the boundary reference level, which may fall
        // on a half level, so it admits no 8-bit encoding under a self-dual span.
        requireThrows<std::invalid_argument>(
            [&]() {
                TopographicConvention rejected{SelfDualSpanImmersion{}, TopographicDomainExtension::ExteriorRing, PixelId{0},
                                               TopographicAltitudeEncoding::UInt8};
                static_cast<void>(MorphologicalTreeFactory::createTreeOfShapes(image, std::move(rejected)));
            },
            "self-dual span immersion over an exterior ring must reject the 8-bit altitude encoding");

        // The published altitude type and the declared encoding must agree.
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(MorphologicalTreeFactory::createTreeOfShapes<ToSGrayLevel>(image)); },
                                             "doubled altitudes must be rejected under an 8-bit declaration");
        requireThrows<std::invalid_argument>(
            [&]() { static_cast<void>(MorphologicalTreeFactory::createTreeOfShapes<std::uint8_t>(image, selfDualSpanConvention())); },
            "8-bit altitudes must be rejected under a doubled declaration");
    }

    return 0;
}
