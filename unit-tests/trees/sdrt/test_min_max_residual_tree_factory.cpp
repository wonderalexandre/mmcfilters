#include "support/TestSupport.hpp"

#include "mmcfilters/trees/sdrt/MinMaxResidualTreeBuilder.hpp"
#include "sdrt_reference/OptimizedUnionFindSelfDualResidualTreeBuilder.hpp"
#include "sdrt_reference/oracle/rag/SingleAdjacencySaturatedResidualTreeOracle.hpp"

#include <array>
#include <cstdint>
#include <span>
#include <vector>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

namespace {

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
                                             sdrt::SaturatedMinMaxLcaPolicy lcaPolicy, sdrt::SaturatedMinMaxFallbackPolicy fallbackPolicy,
                                             sdrt::SaturatedMinMaxBoundaryPolicy boundaryPolicy) {
    auto minTree = MorphologicalTreeFactory::createMinTree(image, adjacency);
    auto maxTree = MorphologicalTreeFactory::createMaxTree(image, adjacency);
    sdrt::MinMaxResidualTreeBuilder<T> builder(adjacency, infinityPixel, sdrt::SdrtTiePolicy::ContrastInvariantSpatial, lcaPolicy, fallbackPolicy,
                                               boundaryPolicy, sdrt::MinMaxResidualEligibilityPolicy::SaturatedOnly);
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

void testFactoryModesAndSemantics() {
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
    sdrt::MinMaxResidualTreeBuilder<std::uint8_t> builder(adjacency, NodeId{0}, sdrt::SdrtTiePolicy::ContrastInvariantSpatial,
                                                          sdrt::SaturatedMinMaxLcaPolicy::ParentClimb, sdrt::SaturatedMinMaxFallbackPolicy::BoundaryMultiSource,
                                                          sdrt::SaturatedMinMaxBoundaryPolicy::IncrementalSmallToLarge,
                                                          sdrt::MinMaxResidualEligibilityPolicy::AllRegionalExtrema);
    builder.build(image, std::move(minTree), std::move(maxTree));
    requireEqual(builder.getStatistics().rejectedExtrema, std::size_t{0}, "unrestricted builder rejects no extrema");
    requireEqual(builder.getStatistics().complementTraversalCertificates, std::size_t{0}, "unrestricted builder performs no complement traversal");
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

        const auto polarityOrdered = MorphologicalTreeFactory::createSelfDualResidualTree(image, adjacency, sdrt::SdrtTiePolicy::MaxBeforeMinThenSpatial);
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
    constexpr std::array boundaryPolicies{sdrt::SaturatedMinMaxBoundaryPolicy::RecomputeFromSupport,
                                          sdrt::SaturatedMinMaxBoundaryPolicy::IncrementalSmallToLarge};
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
                for (const auto boundaryPolicy : boundaryPolicies) {
                    const auto actual = buildResidualWithPolicies(image, adjacency, infinityPixel, lcaPolicy, fallbackPolicy, boundaryPolicy);
                    requireVectorEqual(actual.parent, expectedParent, "configurable policy parents");
                    requireVectorEqual(actual.owner, expectedOwner, "configurable policy owners");
                    requireVectorEqual(actual.altitude, expectedAltitude, "configurable policy altitudes");
                }
            }
        }
    }
}

} // namespace

int main() {
    testFactoryModesAndSemantics();
    testContrastInversion();
    testDirectBuilderPolicies();
    testExhaustiveDifferentialOracles();
    testSharedArbitrarySymmetricAdjacency();
    testConfigurablePolicyEquivalence();
    return 0;
}
