#include "support/TestSupport.hpp"
#include "mmcfilters/attributes/AttributeComputation.hpp"

#include <limits>
#include <memory>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

int main() {
    auto computeArea = [](const MorphologicalTree& tree, NodeId nodeId) {
        auto [attrNames, buffer] = AttributeComputation::computeSingleTopologyAttribute(tree, AREA);
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
    mmcfilters::unit_tests::requireThrows<std::invalid_argument>([]() { (void)ImageUInt8::create(0, 1); }, "zero-row image creation must throw");
    mmcfilters::unit_tests::requireThrows<std::invalid_argument>([]() { (void)ImageUInt8::create(1, 0); }, "zero-column image creation must throw");
    mmcfilters::unit_tests::requireThrows<std::invalid_argument>([]() { (void)ImageUInt8::create(-1, 1); }, "negative-row image creation must throw");
    mmcfilters::unit_tests::requireThrows<std::overflow_error>(
        []() { (void)ImageUInt8::create(std::numeric_limits<int>::max(), 2); },
        "overflowing image creation must throw before allocation");

    auto maxTree = makeWeightedComponentTree(image, true);
    auto minTree = makeWeightedComponentTree(image, false);
    auto [higraParent, higraAltitude] = exportHigraHierarchy(*maxTree);
    auto maxTreeShortcut = MorphologicalTreeFactory::createMaxTree(image);
    auto minTreeShortcut = MorphologicalTreeFactory::createMinTree(image);

    require(static_cast<bool>(maxTree), "max-tree instance must be created");
    require(static_cast<bool>(minTree), "min-tree instance must be created");

    require(maxTree->topology().hasAdjacencyRelation(), "max-tree must expose adjacency relation");
    requireEqual(maxTree->topology().getNumRowsOfImage(), 4, "max-tree image rows");
    requireEqual(maxTree->topology().getNumColsOfImage(), 4, "max-tree image cols");
    require(maxTree->topology().getTreeType() == MorphologicalTreeKind::MAX_TREE, "tree must be a max-tree");
    requireEqual(maxTree->topology().getNumNodes(), 6, "max-tree node count");
    requireEqual(static_cast<int>(maxTree->topology().getLeaves().size()), 1, "max-tree leaves");
    requireVectorEqual(maxTree->topology().getLeaves(), std::vector<NodeId>{5}, "max-tree leaves should be dense NodeIds");

    require(minTree->topology().hasAdjacencyRelation(), "min-tree must expose adjacency relation");
    requireEqual(minTree->topology().getNumRowsOfImage(), 4, "min-tree image rows");
    requireEqual(minTree->topology().getNumColsOfImage(), 4, "min-tree image cols");
    require(minTree->topology().getTreeType() == MorphologicalTreeKind::MIN_TREE, "tree must be a min-tree");
    requireEqual(minTree->topology().getNumNodes(), 6, "min-tree node count");
    requireEqual(static_cast<int>(minTree->topology().getLeaves().size()), 2, "min-tree leaves");
    requireVectorEqual(minTree->topology().getLeaves(), std::vector<NodeId>({4, 5}), "min-tree leaves should be dense NodeIds");

    requireEqual(maxTree->topology().getRoot(), 0, "max-tree root dense node id");
    requireEqual(maxTree->getAltitude(maxTree->topology().getRoot()), 0, "max-tree root altitude");
    requireEqual(computeArea(maxTree->topology(), maxTree->topology().getRoot()), 16, "max-tree root area");
    require(maxTreeShortcut.topology().getTreeType() == MorphologicalTreeKind::MAX_TREE, "createMaxTree shortcut must create a max-tree");
    requireEqual(maxTreeShortcut.topology().getRoot(), maxTree->topology().getRoot(), "createMaxTree shortcut root");
    requireEqual(computeArea(maxTreeShortcut.topology(), maxTreeShortcut.topology().getRoot()), 16, "createMaxTree shortcut root area");

    requireEqual(minTree->topology().getRoot(), 0, "min-tree root dense node id");
    requireEqual(minTree->getAltitude(minTree->topology().getRoot()), 5, "min-tree root altitude");
    requireEqual(computeArea(minTree->topology(), minTree->topology().getRoot()), 16, "min-tree root area");
    require(minTreeShortcut.topology().getTreeType() == MorphologicalTreeKind::MIN_TREE, "createMinTree shortcut must create a min-tree");
    requireEqual(minTreeShortcut.topology().getRoot(), minTree->topology().getRoot(), "createMinTree shortcut root");
    requireEqual(computeArea(minTreeShortcut.topology(), minTreeShortcut.topology().getRoot()), 16, "createMinTree shortcut root area");

    auto maxReconstruction = maxTree->reconstructionImage();
    auto minReconstruction = minTree->reconstructionImage();
    requireImageShape(maxReconstruction, 4, 4);
    requireImageShape(minReconstruction, 4, 4);
    requireVectorEqual(collectImageValues(maxReconstruction), collectImageValues(image), "max-tree reconstruction values");
    requireVectorEqual(collectImageValues(minReconstruction), collectImageValues(image), "min-tree reconstruction values");

    auto importedFromHigra = MorphologicalTreeFactory::createFromHigraParent(
        std::span<const NodeId>(higraParent),
        std::span<const std::uint8_t>(higraAltitude),
        4,
        4,
        MorphologicalTreeKind::MAX_TREE,
        AdjacencyRelation(4, 4, 1.5));
    requireEqual(importedFromHigra.topology().getNumHigraNodes(), 22, "Higra import node-id domain size");
    requireEqual(importedFromHigra.topology().getHigraNodeId(3), 19, "Higra import slot->node-id mapping");
    requireEqual(computeArea(importedFromHigra.topology(), importedFromHigra.topology().getRoot()), 16, "Higra import area");
    requireImageShape(importedFromHigra.reconstructionImage(), 4, 4);
    requireVectorEqual(collectImageValues(importedFromHigra.reconstructionImage()), collectImageValues(image), "Higra import reconstruction values");
    auto [reexportedHigraParent, reexportedHigraAltitude] = importedFromHigra.exportHigraHierarchy();
    requireVectorEqual(reexportedHigraParent, higraParent, "Higra export/import parent round-trip");
    requireVectorEqual(reexportedHigraAltitude, higraAltitude, "Higra export/import altitude round-trip");

    auto imageBuiltTopology = MorphologicalTreeFactory::createMaxTree(image);
    requireThrows([&]() { imageBuiltTopology.topology().getNumHigraNodes(); }, "image-built topology must not expose a Higra node-id domain");

    auto mutatedImport = MorphologicalTreeFactory::createFromHigraParent(
        std::span<const NodeId>(higraParent),
        std::span<const std::uint8_t>(higraAltitude),
        4,
        4,
        MorphologicalTreeKind::MAX_TREE,
        AdjacencyRelation(4, 4, 1.5));
    mutatedImport.mergeNodeIntoParent(0);
    requireThrows([&]() { mutatedImport.topology().getNumHigraNodes(); }, "mutated Higra import must drop the original Higra node-id domain");
    requireEqual(mutatedImport.topology().getHigraNodeId(3), InvalidNode, "mutated Higra import must not expose stale slot->node-id mapping");

    {
        auto sparseTree = makeWeightedComponentTree(image, true);
        sparseTree->mergeNodeIntoParent(5);
        auto [sparseHigraParent, sparseHigraAltitude] = sparseTree->exportHigraHierarchy();
        requireEqual(static_cast<int>(sparseHigraParent.size()), 16 + 5, "Higra export must compact dead node slots");
        auto sparseImported = MorphologicalTreeFactory::createFromHigraParent(
            std::span<const NodeId>(sparseHigraParent),
            std::span<const std::uint8_t>(sparseHigraAltitude),
            4,
            4,
            MorphologicalTreeKind::MAX_TREE,
            AdjacencyRelation(4, 4, 1.5));
        requireEqual(sparseImported.topology().getNumNodes(), 5, "compact Higra round-trip node count");
        requireEqual(sparseImported.topology().getNumInternalNodeSlots(), 5, "compact Higra round-trip slot count");
        requireEqual(computeArea(sparseImported.topology(), sparseImported.topology().getRoot()), 16, "compact Higra round-trip area");
        auto [sparseReexportedParent, sparseReexportedAltitude] = sparseImported.exportHigraHierarchy();
        requireVectorEqual(sparseReexportedParent, sparseHigraParent, "compact Higra parent round-trip");
        requireVectorEqual(sparseReexportedAltitude, sparseHigraAltitude, "compact Higra altitude round-trip");
    }

    return 0;
}
