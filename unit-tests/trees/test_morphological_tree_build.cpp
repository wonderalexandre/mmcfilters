#include "support/TestSupport.hpp"
#include "mmcfilters/attributes/AttributeComputedIncrementally.hpp"

#include <memory>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

int main() {
    auto computeArea = [](MorphologicalTree& tree, NodeId nodeId) {
        auto [attrNames, buffer] = AttributeComputedIncrementally::computeSingleAttribute(tree, AREA);
        return static_cast<int>(buffer[attrNames.linearIndex(nodeId, AREA)]);
    };
    auto requireThrows = [](auto&& fn, const std::string& label) {
        bool threw = false;
        try {
            fn();
        } catch (const std::exception&) {
            threw = true;
        }
        require(threw, label);
    };

    auto image = makeComponentTreeFixture();
    auto maxTree = makeWeightedComponentTree(image, true);
    auto minTree = makeWeightedComponentTree(image, false);
    auto [higraParent, higraAltitude] = exportHigraHierarchy(*maxTree);

    require(static_cast<bool>(maxTree), "max-tree instance must be created");
    require(static_cast<bool>(minTree), "min-tree instance must be created");

    require(maxTree->tree.hasAdjacencyRelation(), "max-tree must expose adjacency relation");
    requireEqual(maxTree->tree.getNumRowsOfImage(), 4, "max-tree image rows");
    requireEqual(maxTree->tree.getNumColsOfImage(), 4, "max-tree image cols");
    require(maxTree->tree.isMaxtree(), "tree must be a max-tree");
    requireEqual(maxTree->tree.getNumNodes(), 6, "max-tree node count");
    requireEqual(static_cast<int>(maxTree->tree.getLeaves().size()), 1, "max-tree leaves");
    requireVectorEqual(maxTree->tree.getLeaves(), std::vector<NodeId>{5}, "max-tree leaves should be dense NodeIds");

    require(minTree->tree.hasAdjacencyRelation(), "min-tree must expose adjacency relation");
    requireEqual(minTree->tree.getNumRowsOfImage(), 4, "min-tree image rows");
    requireEqual(minTree->tree.getNumColsOfImage(), 4, "min-tree image cols");
    require(!minTree->tree.isMaxtree(), "tree must be a min-tree");
    requireEqual(minTree->tree.getNumNodes(), 6, "min-tree node count");
    requireEqual(static_cast<int>(minTree->tree.getLeaves().size()), 2, "min-tree leaves");
    requireVectorEqual(minTree->tree.getLeaves(), std::vector<NodeId>({4, 5}), "min-tree leaves should be dense NodeIds");

    requireEqual(maxTree->tree.getRoot(), 0, "max-tree root dense node id");
    requireEqual(maxTree->getAltitude(maxTree->tree.getRoot()), 0, "max-tree root altitude");
    requireEqual(computeArea(maxTree->tree, maxTree->tree.getRoot()), 16, "max-tree root area");

    requireEqual(minTree->tree.getRoot(), 0, "min-tree root dense node id");
    requireEqual(minTree->getAltitude(minTree->tree.getRoot()), 5, "min-tree root altitude");
    requireEqual(computeArea(minTree->tree, minTree->tree.getRoot()), 16, "min-tree root area");

    auto maxReconstruction = maxTree->reconstructionImage();
    auto minReconstruction = minTree->reconstructionImage();
    requireImageShape(maxReconstruction, 4, 4);
    requireImageShape(minReconstruction, 4, 4);
    requireVectorEqual(collectImageValues(maxReconstruction), collectImageValues(image), "max-tree reconstruction values");
    requireVectorEqual(collectImageValues(minReconstruction), collectImageValues(image), "min-tree reconstruction values");

    auto importedFromHigra = WeightedMorphologicalTree::createFromHigraParent(higraParent, higraAltitude, 4, 4, MorphologicalTree::MAX_TREE, AdjacencyRelation(4, 4, 1.5));
    requireEqual(importedFromHigra.tree.getNumHigraNodes(), 22, "Higra import node-id domain size");
    requireEqual(importedFromHigra.tree.getHigraNodeId(3), 19, "Higra import slot->node-id mapping");
    requireEqual(computeArea(importedFromHigra.tree, importedFromHigra.tree.getRoot()), 16, "Higra import area");
    requireImageShape(importedFromHigra.reconstructionImage(), 4, 4);
    requireVectorEqual(collectImageValues(importedFromHigra.reconstructionImage()), collectImageValues(image), "Higra import reconstruction values");
    auto [reexportedHigraParent, reexportedHigraAltitude] = importedFromHigra.exportHigraHierarchy();
    requireVectorEqual(reexportedHigraParent, higraParent, "Higra export/import parent round-trip");
    requireVectorEqual(reexportedHigraAltitude, higraAltitude, "Higra export/import altitude round-trip");

    auto imageBuiltTopology = MorphologicalTree::createComponentTree(image, true);
    requireThrows([&]() { imageBuiltTopology.getNumHigraNodes(); }, "image-built topology must not expose a Higra node-id domain");

    auto mutatedImport = WeightedMorphologicalTree::createFromHigraParent(higraParent, higraAltitude, 4, 4, MorphologicalTree::MAX_TREE, AdjacencyRelation(4, 4, 1.5));
    mutatedImport.mergeNodeIntoParent(0);
    requireThrows([&]() { mutatedImport.tree.getNumHigraNodes(); }, "mutated Higra import must drop the original Higra node-id domain");
    requireEqual(mutatedImport.tree.getHigraNodeId(3), InvalidNode, "mutated Higra import must not expose stale slot->node-id mapping");

    {
        auto sparseTree = makeWeightedComponentTree(image, true);
        sparseTree->mergeNodeIntoParent(5);
        auto [sparseHigraParent, sparseHigraAltitude] = sparseTree->exportHigraHierarchy();
        requireEqual(static_cast<int>(sparseHigraParent.size()), 16 + 5, "Higra export must compact dead node slots");
        auto sparseImported = WeightedMorphologicalTree::createFromHigraParent(
            sparseHigraParent,
            sparseHigraAltitude,
            4,
            4,
            MorphologicalTree::MAX_TREE,
            AdjacencyRelation(4, 4, 1.5));
        requireEqual(sparseImported.tree.getNumNodes(), 5, "compact Higra round-trip node count");
        requireEqual(sparseImported.tree.getNumInternalNodeSlots(), 5, "compact Higra round-trip slot count");
        requireEqual(computeArea(sparseImported.tree, sparseImported.tree.getRoot()), 16, "compact Higra round-trip area");
        auto [sparseReexportedParent, sparseReexportedAltitude] = sparseImported.exportHigraHierarchy();
        requireVectorEqual(sparseReexportedParent, sparseHigraParent, "compact Higra parent round-trip");
        requireVectorEqual(sparseReexportedAltitude, sparseHigraAltitude, "compact Higra altitude round-trip");
    }

    return 0;
}
