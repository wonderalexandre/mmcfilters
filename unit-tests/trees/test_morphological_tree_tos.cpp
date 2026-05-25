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
    require(tree->getTreeType() == MorphologicalTreeKind::TREE_OF_SHAPES, "tree type");
    require(tree->getNumNodes() > 0, "tree of shapes must contain nodes");
    require(tree->getRoot() != InvalidNode, "tree of shapes root must be valid");
    require(!tree->hasAdjacencyRelation(), "tree of shapes must not expose adjacency relation by default");
    require(tree->getAdjacencyRelation() == nullptr, "tree of shapes adjacency pointer must be null");
    require(tree->hasTreeOfShapesAdjacencyPolicy(), "tree of shapes must expose auxiliary min/max adjacency policy");
    requireNear(tree->getTreeOfShapesMinTreeAdjacencyRadius(), 1.0, 0.0, "Min4cMax8c min-tree radius");
    requireNear(tree->getTreeOfShapesMaxTreeAdjacencyRadius(), 1.5, 0.0, "Min4cMax8c max-tree radius");
    require(tree->getTreeOfShapesMinTreeAdjacencyRelation() != nullptr, "Min4cMax8c min-tree adjacency relation");
    require(tree->getTreeOfShapesMaxTreeAdjacencyRelation() != nullptr, "Min4cMax8c max-tree adjacency relation");
    requireEqual(tree->getTreeOfShapesMinTreeAdjacencyRelation()->getSize(), 5, "Min4cMax8c min-tree adjacency relation size");
    requireEqual(tree->getTreeOfShapesMaxTreeAdjacencyRelation()->getSize(), 9, "Min4cMax8c max-tree adjacency relation size");
    requireEqual(tree->getNumRowsOfImage(), 3, "tree of shapes image rows");
    requireEqual(tree->getNumColsOfImage(), 3, "tree of shapes image cols");

    auto dualTree = makeTreeOfShapes(image, ToSInterpolation::SelfDual);
    require(dualTree->hasTreeOfShapesAdjacencyPolicy(), "self-dual tree of shapes must expose auxiliary adjacency policy");
    requireNear(dualTree->getTreeOfShapesMinTreeAdjacencyRadius(), 1.0, 0.0, "SelfDual min-tree radius");
    requireNear(dualTree->getTreeOfShapesMaxTreeAdjacencyRadius(), 1.0, 0.0, "SelfDual max-tree radius");

    auto min8Max4Tree = makeTreeOfShapes(image, ToSInterpolation::Min8cMax4c);
    require(min8Max4Tree->hasTreeOfShapesAdjacencyPolicy(), "Min8cMax4c tree of shapes must expose auxiliary adjacency policy");
    requireNear(min8Max4Tree->getTreeOfShapesMinTreeAdjacencyRadius(), 1.5, 0.0, "Min8cMax4c min-tree radius");
    requireNear(min8Max4Tree->getTreeOfShapesMaxTreeAdjacencyRadius(), 1.0, 0.0, "Min8cMax4c max-tree radius");
    requireEqual(min8Max4Tree->getTreeOfShapesMinTreeAdjacencyRelation()->getSize(), 9, "Min8cMax4c min-tree adjacency relation size");
    requireEqual(min8Max4Tree->getTreeOfShapesMaxTreeAdjacencyRelation()->getSize(), 5, "Min8cMax4c max-tree adjacency relation size");

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
    requireEqual(singleTree->getProperPartOwner(0), singleTree->getRoot(), "single-pixel ToS owner");
    auto singleReconstruction = singleWeighted->reconstructionImage();
    requireImageShape(singleReconstruction, 1, 1);
    requireVectorEqual(collectImageValues(singleReconstruction), std::vector<uint8_t>{5}, "single-pixel ToS reconstruction");

    requireThrows<std::invalid_argument>(
        [&]() { static_cast<void>(ImageUInt8::create(0, 0)); },
        "empty image creation must throw");

    return 0;
}
