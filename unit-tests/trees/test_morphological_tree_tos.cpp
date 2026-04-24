#include "support/TestSupport.hpp"

#include <memory>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

int main() {
    auto image = makeImage(
        3,
        3,
        {
            1, 2, 1,
            2, 3, 2,
            1, 2, 1,
        }
    );

    auto tree = makeTreeOfShapes(image, ToSInterpolation::Min4cMax8c);
    auto weighted = makeWeightedTreeOfShapes(image, ToSInterpolation::Min4cMax8c);

    require(static_cast<bool>(tree), "tree of shapes instance must be created");
    requireEqual(tree->getTreeType(), 2, "tree type");
    require(tree->getNumNodes() > 0, "tree of shapes must contain nodes");
    require(tree->getRoot() != InvalidNode, "tree of shapes root must be valid");
    require(!tree->hasAdjacencyRelation(), "tree of shapes must not expose adjacency relation by default");
    require(tree->getAdjacencyRelation() == nullptr, "tree of shapes adjacency pointer must be null");
    requireEqual(tree->getNumRowsOfImage(), 3, "tree of shapes image rows");
    requireEqual(tree->getNumColsOfImage(), 3, "tree of shapes image cols");

    auto reconstruction = weighted->reconstructionImage();
    requireImageShape(reconstruction, 3, 3);
    requireVectorEqual(
        collectImageValues(reconstruction),
        std::vector<uint8_t>{1, 2, 1, 2, 3, 2, 1, 2, 1},
        "tree of shapes reconstruction values"
    );

    auto singlePixel = makeImage(1, 1, {5});
    auto singleTree = makeTreeOfShapes(singlePixel);
    auto singleWeighted = makeWeightedTreeOfShapes(singlePixel);
    require(static_cast<bool>(singleTree), "single-pixel self-dual tree of shapes instance must be created");
    requireEqual(singleTree->getNumRowsOfImage(), 1, "single-pixel ToS rows");
    requireEqual(singleTree->getNumColsOfImage(), 1, "single-pixel ToS cols");
    require(singleTree->getRoot() != InvalidNode, "single-pixel ToS root must be valid");
    requireEqual(singleTree->getSmallestComponent(0), singleTree->getRoot(), "single-pixel ToS owner");
    auto singleReconstruction = singleWeighted->reconstructionImage();
    requireImageShape(singleReconstruction, 1, 1);
    requireVectorEqual(collectImageValues(singleReconstruction), std::vector<uint8_t>{5}, "single-pixel ToS reconstruction");

    auto emptyImage = ImageUInt8::create(0, 0);
    requireThrows<std::invalid_argument>(
        [&]() { static_cast<void>(MorphologicalTree::createComponentTree(emptyImage, true)); },
        "empty component tree must throw");
    requireThrows<std::invalid_argument>(
        [&]() { static_cast<void>(MorphologicalTree::createTreeOfShapes(emptyImage)); },
        "empty tree of shapes must throw");
    requireThrows<std::invalid_argument>(
        [&]() { static_cast<void>(WeightedMorphologicalTree::createComponentTree(emptyImage, true)); },
        "empty weighted component tree must throw");
    requireThrows<std::invalid_argument>(
        [&]() { static_cast<void>(WeightedMorphologicalTree::createTreeOfShapes(emptyImage)); },
        "empty weighted tree of shapes must throw");

    return 0;
}
