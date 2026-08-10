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
    supports.reserve(static_cast<std::size_t>(tree.getNumNodes()));
    for (NodeId nodeId : tree.getAliveNodeIds()) {
        std::vector<NodeId> support = collectNodeIds(tree.getConnectedComponent(nodeId));
        std::sort(support.begin(), support.end());
        supports.push_back(std::move(support));
    }
    std::sort(supports.begin(), supports.end());
    return supports;
}

ImageUInt8Ptr makeImageFromVector(int rows, int cols, const std::vector<uint8_t>& values) {
    requireEqual(static_cast<int>(values.size()), rows * cols, "vector image buffer size");
    auto image = ImageUInt8::create(rows, cols);
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

    auto tree = makeTreeOfShapes(image, ToSInterpolation::Min4cMax8c);
    auto weighted = makeWeightedTreeOfShapes(image, ToSInterpolation::Min4cMax8c);

    require(static_cast<bool>(tree), "tree of shapes instance must be created");
    require(tree->getDescriptiveKind() == MorphologicalTreeKind::TREE_OF_SHAPES, "tree type");
    require(tree->getNumNodes() > 0, "tree of shapes must contain nodes");
    require(tree->getRoot() != InvalidNode, "tree of shapes root must be valid");
    require(!tree->hasUniformGridAdjacency2D(), "tree of shapes must not expose adjacency relation by default");
    require(tree->getUniformGridAdjacency2D() == nullptr, "tree of shapes adjacency pointer must be null");
    require(tree->hasDirectionalGridAdjacency2D(), "tree of shapes must expose generic dual adjacency context");
    require(tree->hasDirectionalGridAdjacency2D(), "tree of shapes must expose directional adjacency capability");
    require(tree->getAdjacencyMode() == AdjacencyMode::DIRECTIONAL, "tree of shapes adjacency mode");
    require(tree->getAltitudeOrder() == AltitudeOrder::UNCONSTRAINED, "tree of shapes altitude order");
    requireNear(tree->getDecreasingGridAdjacency2D()->getRadius(), 1.0, 0.0, "Min4cMax8c decreasing radius");
    requireNear(tree->getIncreasingGridAdjacency2D()->getRadius(), 1.5, 0.0, "Min4cMax8c increasing radius");
    requireEqual(tree->getDecreasingGridAdjacency2D()->getSize(), 5, "Min4cMax8c decreasing adjacency size");
    requireEqual(tree->getIncreasingGridAdjacency2D()->getSize(), 9, "Min4cMax8c increasing adjacency size");
    requireEqual(tree->getNumRowsOfGridDomain2D(), 3, "tree of shapes image rows");
    requireEqual(tree->getNumColsOfGridDomain2D(), 3, "tree of shapes image cols");
    for (NodeId nodeId : tree->getAliveNodeIds()) {
        require(!collectNodeIds(tree->getConnectedComponent(nodeId)).empty(), "every committed tree-of-shapes node must have non-empty subtree support");
    }

    auto dualTree = makeTreeOfShapes(image, ToSInterpolation::SelfDual);
    require(dualTree->hasDirectionalGridAdjacency2D(), "self-dual tree of shapes must expose generic dual adjacency context");
    requireNear(dualTree->getDecreasingGridAdjacency2D()->getRadius(), 1.0, 0.0, "SelfDual decreasing radius");
    requireNear(dualTree->getIncreasingGridAdjacency2D()->getRadius(), 1.0, 0.0, "SelfDual increasing radius");

    auto min8Max4Tree = makeTreeOfShapes(image, ToSInterpolation::Min8cMax4c);
    require(min8Max4Tree->hasDirectionalGridAdjacency2D(), "Min8cMax4c tree of shapes must expose generic dual adjacency context");
    requireNear(min8Max4Tree->getDecreasingGridAdjacency2D()->getRadius(), 1.5, 0.0, "Min8cMax4c decreasing radius");
    requireNear(min8Max4Tree->getIncreasingGridAdjacency2D()->getRadius(), 1.0, 0.0, "Min8cMax4c increasing radius");
    requireEqual(min8Max4Tree->getDecreasingGridAdjacency2D()->getSize(), 9, "Min8cMax4c decreasing adjacency size");
    requireEqual(min8Max4Tree->getIncreasingGridAdjacency2D()->getSize(), 5, "Min8cMax4c increasing adjacency size");

    auto reconstruction = weighted->reconstructionImage();
    requireImageShape(reconstruction, 3, 3);
    requireVectorEqual(collectImageValues(reconstruction), std::vector<uint8_t>{1, 2, 1, 2, 3, 2, 1, 2, 1}, "tree of shapes reconstruction values");

    auto singlePixel = makeImage(1, 1, {5});
    auto singleTree = makeTreeOfShapes(singlePixel);
    auto singleWeighted = makeWeightedTreeOfShapes(singlePixel);
    require(static_cast<bool>(singleTree), "single-pixel self-dual tree of shapes instance must be created");
    requireEqual(singleTree->getNumRowsOfGridDomain2D(), 1, "single-pixel ToS rows");
    requireEqual(singleTree->getNumColsOfGridDomain2D(), 1, "single-pixel ToS cols");
    require(singleTree->getRoot() != InvalidNode, "single-pixel ToS root must be valid");
    requireEqual(singleTree->getProperPartOwner(0), singleTree->getRoot(), "single-pixel ToS owner");
    auto singleReconstruction = singleWeighted->reconstructionImage();
    requireImageShape(singleReconstruction, 1, 1);
    requireVectorEqual(collectImageValues(singleReconstruction), std::vector<uint8_t>{5}, "single-pixel ToS reconstruction");
    for (ToSInterpolation interpolation : {ToSInterpolation::SelfDual, ToSInterpolation::Min4cMax8c, ToSInterpolation::Min8cMax4c}) {
        auto unpaddedSingle =
            MorphologicalTreeFactory::createTreeOfShapes(singlePixel, TreeOfShapesProducerOptions{interpolation, ToSPaddingPolicy::NoPadding, 0, 0});
        requireEqual(unpaddedSingle.topology().getNumNodes(), 1, "single-pixel unpadded ToS node count");
        requireVectorEqual(collectImageValues(unpaddedSingle.reconstructionImage()), std::vector<std::uint8_t>{5}, "single-pixel unpadded ToS reconstruction");
    }

    {
        auto virtualRootImage = makeImage(2, 2,
                                          {
                                              1,
                                              1,
                                              0,
                                              0,
                                          });
        auto virtualRootWeighted = makeWeightedTreeOfShapes(virtualRootImage, ToSInterpolation::SelfDual);
        const MorphologicalTree& virtualRootTree = virtualRootWeighted->topology();
        const NodeId virtualRoot = virtualRootTree.getRoot();
        const NodeId upperShape = virtualRootTree.getProperPartOwner(0);
        const NodeId lowerShape = virtualRootTree.getProperPartOwner(2);

        requireEqual(virtualRootTree.getNumNodes(), 3, "virtual-root ToS node count");
        requireEqual(virtualRootTree.getNumProperParts(virtualRoot), 0, "virtual-root ToS root direct proper part");
        require(virtualRootTree.isStructuralNode(virtualRoot), "virtual-root ToS root must be derived as structural");
        requireEqual(virtualRootTree.getNumChildren(virtualRoot), 2, "virtual-root ToS root children");
        requireEqual(computeAreaViaAttributeFacade(virtualRootTree, virtualRoot), 4, "virtual-root ToS cumulative area");
        require(upperShape != lowerShape, "virtual-root ToS must preserve both child shapes");
        requireEqual(computeAreaViaAttributeFacade(virtualRootTree, upperShape), 2, "upper shape area");
        requireEqual(computeAreaViaAttributeFacade(virtualRootTree, lowerShape), 2, "lower shape area");
        requireEqual(virtualRootTree.getNodeParent(upperShape), virtualRoot, "upper shape parent");
        requireEqual(virtualRootTree.getNodeParent(lowerShape), virtualRoot, "lower shape parent");
        requireVectorEqual(collectImageValues(virtualRootWeighted->reconstructionImage()), std::vector<uint8_t>{1, 1, 0, 0}, "virtual-root ToS reconstruction");

        auto exactVirtualRoot = MorphologicalTreeFactory::createTreeOfShapesExact(virtualRootImage);
        static_assert(std::is_same_v<decltype(exactVirtualRoot), WeightedMorphologicalTree<ToSGrayLevel>>);
        const MorphologicalTree& exactTopology = exactVirtualRoot.topology();
        requireEqual(exactTopology.getNumNodes(), virtualRootTree.getNumNodes(), "exact ToS topology node count");
        requireEqual(exactTopology.getRoot(), virtualRoot, "exact ToS root");
        for (NodeId nodeId : virtualRootTree.getAliveNodeIds()) {
            requireEqual(exactTopology.getNodeParent(nodeId), virtualRootTree.getNodeParent(nodeId), "exact ToS parent topology");
        }
        for (NodeId elementId = 0; elementId < virtualRootTree.getNumTotalProperParts(); ++elementId) {
            requireEqual(exactTopology.getProperPartOwner(elementId), virtualRootTree.getProperPartOwner(elementId), "exact ToS proper-part owner");
        }
        requireEqual(exactVirtualRoot.getAltitude(virtualRoot), static_cast<ToSGrayLevel>(1), "exact ToS preserves a half-level virtual root");
        requireEqual(virtualRootWeighted->getAltitude(virtualRoot), static_cast<std::uint8_t>(0), "compatibility ToS quantizes a half-level virtual root");
        const NodeId lowerOwner = exactTopology.getProperPartOwner(2);
        require(lowerOwner != virtualRoot, "exact ToS lower shape must remain below its virtual root");
        requireEqual(exactVirtualRoot.getAltitude(lowerOwner), static_cast<ToSGrayLevel>(0), "exact ToS lower source altitude");
        requireEqual(exactVirtualRoot.getAltitude(virtualRoot), static_cast<ToSGrayLevel>(1), "exact ToS distinguishes the quantized parent arc");
        requireEqual(virtualRootWeighted->getAltitude(lowerOwner), virtualRootWeighted->getAltitude(virtualRoot),
                     "compatibility ToS collapses only the half-level altitude arc");
        requireVectorEqual(collectImageValues(exactVirtualRoot.reconstructionImage()), std::vector<ToSGrayLevel>{2, 2, 0, 0},
                           "exact ToS reconstruction uses doubled gray units");

        const auto [exportedParent, exportedAltitude] = virtualRootWeighted->exportHigraHierarchy();
        requireEqual(static_cast<int>(exportedParent.size()), 7, "virtual-root ToS Higra export vertex count");
        requireEqual(exportedParent.size(), exportedAltitude.size(), "virtual-root ToS Higra export altitude size");

        auto invertedImage = makeImage(2, 2,
                                       {
                                           1,
                                           1,
                                           2,
                                           2,
                                       });
        auto invertedWeighted = makeWeightedTreeOfShapes(invertedImage, ToSInterpolation::SelfDual);
        require(supportFamily(virtualRootTree) == supportFamily(invertedWeighted->topology()),
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
        for (ToSInterpolation interpolation : {ToSInterpolation::SelfDual, ToSInterpolation::Min4cMax8c, ToSInterpolation::Min8cMax4c}) {
            const TreeOfShapesProducerOptions options{interpolation, ToSPaddingPolicy::NoPadding, 0, 0};
            auto unpadded = MorphologicalTreeFactory::createTreeOfShapes(unpaddedImage, options);
            requireEqual(unpadded.topology().getNumRowsOfGridDomain2D(), 2, "unpadded ToS published rows");
            requireEqual(unpadded.topology().getNumColsOfGridDomain2D(), 3, "unpadded ToS published columns");
            requireEqual(unpadded.topology().getNumTotalProperParts(), 6, "unpadded ToS proper-part domain size");
            requireVectorEqual(collectImageValues(unpadded.reconstructionImage()), std::vector<std::uint8_t>{0, 2, 1, 2, 1, 0}, "unpadded ToS reconstruction");
            auto exactUnpadded = MorphologicalTreeFactory::createTreeOfShapesExact(unpaddedImage, options);
            requireVectorEqual(collectImageValues(exactUnpadded.reconstructionImage()), std::vector<ToSGrayLevel>{0, 4, 2, 4, 2, 0},
                               "exact unpadded ToS reconstruction");
            require(supportFamily(exactUnpadded.topology()) == supportFamily(unpadded.topology()),
                    "exact altitude encoding must not alter unpadded ToS topology");

            TreeOfShapesProducer producer(options);
            const TreeOfShapesBuildResult result = producer.build(unpaddedImage);
            requireEqual(result.properPartOwner.size(), std::size_t{6}, "unpadded producer publishes only source-domain owners");
            requireEqual(result.numRows, 2, "unpadded producer source rows");
            requireEqual(result.numCols, 3, "unpadded producer source columns");
        }
    }

    // Exhaustive ternary small-domain properties include half-level exterior
    // medians and branching shapes that need virtual projected nodes.
    for (int cols : {2, 3}) {
        const int numPixels = 2 * cols;
        int numConfigurations = 1;
        for (int pixelId = 0; pixelId < numPixels; ++pixelId) {
            numConfigurations *= 3;
        }

        for (int code = 0; code < numConfigurations; ++code) {
            int digits = code;
            std::vector<uint8_t> values(static_cast<std::size_t>(numPixels), uint8_t{0});
            std::vector<uint8_t> inverted(static_cast<std::size_t>(numPixels), uint8_t{0});
            for (int pixelId = 0; pixelId < numPixels; ++pixelId) {
                values[static_cast<std::size_t>(pixelId)] = static_cast<uint8_t>(digits % 3);
                inverted[static_cast<std::size_t>(pixelId)] = static_cast<uint8_t>(2 - values[static_cast<std::size_t>(pixelId)]);
                digits /= 3;
            }
            auto originalTree = makeTreeOfShapes(makeImageFromVector(2, cols, values), ToSInterpolation::SelfDual);
            auto invertedTree = makeTreeOfShapes(makeImageFromVector(2, cols, inverted), ToSInterpolation::SelfDual);
            require(supportFamily(*originalTree) == supportFamily(*invertedTree),
                    "exhaustive ternary 2x" + std::to_string(cols) + " SelfDual support-family invariance, code " + std::to_string(code));

            auto min4Max8Tree = makeTreeOfShapes(makeImageFromVector(2, cols, values), ToSInterpolation::Min4cMax8c);
            auto dualMin8Max4Tree = makeTreeOfShapes(makeImageFromVector(2, cols, inverted), ToSInterpolation::Min8cMax4c);
            require(supportFamily(*min4Max8Tree) == supportFamily(*dualMin8Max4Tree),
                    "exhaustive ternary 2x" + std::to_string(cols) + " dual-connectivity support-family invariance, code " + std::to_string(code));
        }
    }

    if constexpr (contract::validationsEnabled) {
        requireThrows<std::invalid_argument>([&]() { static_cast<void>(ImageUInt8::create(0, 0)); }, "empty image creation must throw");
    }

    return 0;
}
