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

    auto tree = std::make_shared<MorphologicalTree>(image, ToSInterpolation::Min4cMax8c);
    auto weighted = std::make_shared<WeightedMorphologicalTree>(image, ToSInterpolation::Min4cMax8c);

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

    return 0;
}
