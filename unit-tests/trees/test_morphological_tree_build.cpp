#include "support/TestSupport.hpp"
#include "mmcfilters/attributes/AttributeComputation.hpp"

#include <limits>
#include <memory>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

int main() {
    auto computeArea = [](const MorphologicalTree& tree, NodeId nodeId) {
        auto [attrNames, buffer] = AttributeComputation::computeSingleTopologyAttribute(tree, Area);
        return static_cast<int>(buffer[attrNames.linearIndex(nodeId, Area)]);
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
    if constexpr (contract::validationsEnabled) {
        mmcfilters::unit_tests::requireThrows<std::invalid_argument>([]() { (void)ImageUInt8::create(0, 1); }, "zero-row image creation must throw");
        mmcfilters::unit_tests::requireThrows<std::invalid_argument>([]() { (void)ImageUInt8::create(1, 0); }, "zero-column image creation must throw");
        mmcfilters::unit_tests::requireThrows<std::invalid_argument>([]() { (void)ImageUInt8::create(-1, 1); }, "negative-row image creation must throw");
        mmcfilters::unit_tests::requireThrows<std::overflow_error>([]() { (void)ImageUInt8::create(std::numeric_limits<int>::max(), 2); },
                                                                   "overflowing image creation must throw before allocation");
        mmcfilters::unit_tests::requireThrows<std::invalid_argument>([]() { (void)RegularGridAdjacency2D(4, 4, -1.0); },
                                                                     "negative adjacency radius must throw before stencil construction");
        mmcfilters::unit_tests::requireThrows<std::invalid_argument>(
            []() { (void)RegularGridAdjacency2D(4, 4, std::numeric_limits<double>::quiet_NaN()); },
            "NaN adjacency radius must throw before stencil construction");
        mmcfilters::unit_tests::requireThrows<std::invalid_argument>(
            []() { (void)RegularGridAdjacency2D(4, 4, std::numeric_limits<double>::infinity()); },
            "infinite adjacency radius must throw before stencil construction");
        mmcfilters::unit_tests::requireThrows<std::invalid_argument>(
            []() { (void)RegularGridAdjacency2D(4, 4, std::numeric_limits<double>::max()); },
            "unrepresentable finite adjacency radius must throw before stencil construction");
        mmcfilters::unit_tests::requireThrows<std::invalid_argument>([]() { (void)RegularGridAdjacency2D(-1, 4, 1.0); },
                                                                     "negative adjacency dimensions must throw before stencil construction");
    }

    auto maxTree = makeValuedComponentTree(image, true);
    auto minTree = makeValuedComponentTree(image, false);
    auto [higraParent, higraAltitude] = exportHigraHierarchy(*maxTree);
    auto maxTreeShortcut = MorphologicalTreeFactory::createMaxTree(image);
    auto minTreeShortcut = MorphologicalTreeFactory::createMinTree(image);
    const auto collectParents = [](const MorphologicalTree& tree) {
        std::vector<NodeId> parent;
        parent.reserve(static_cast<std::size_t>(tree.numInternalNodeSlots()));
        for (NodeId node = 0; node < tree.numInternalNodeSlots(); ++node) {
            parent.push_back(tree.parent(node));
        }
        return parent;
    };
    const auto collectOwners = [](const MorphologicalTree& tree) {
        std::vector<NodeId> smallestNodeId;
        smallestNodeId.reserve(static_cast<std::size_t>(tree.numPixels()));
        for (PixelId pixel = 0; pixel < tree.numPixels(); ++pixel) {
            smallestNodeId.push_back(tree.smallestNode(pixel));
        }
        return smallestNodeId;
    };

    require(static_cast<bool>(maxTree), "max-tree instance must be created");
    require(static_cast<bool>(minTree), "min-tree instance must be created");

    require(maxTree->topology().sharedAdjacencyContext() != nullptr, "max-tree must retain its shared-adjacency context");
    requireEqual(maxTree->topology().numRows(), 4, "max-tree image rows");
    requireEqual(maxTree->topology().numColumns(), 4, "max-tree image columns");
    require(maxTree->topology().kind() == MorphologicalTreeKind::MaxTree, "tree must be a max-tree");
    require(maxTree->topology().nodeAltitudeOrder() == NodeAltitudeOrder::Increasing, "max-tree altitude order capability");
    require(std::holds_alternative<SharedAdjacencyContext>(maxTree->topology().constructionContext()), "max-tree shared-adjacency context");
    requireEqual(maxTree->topology().numNodes(), 6, "max-tree node count");
    requireEqual(maxTree->topology().getMutationVersion(), std::size_t{7}, "max-tree construction mutation version");
    requireVectorEqual(collectParents(maxTree->topology()), std::vector<NodeId>{0, 0, 1, 2, 3, 4}, "max-tree native parent buffer");
    requireVectorEqual(collectOwners(maxTree->topology()), std::vector<NodeId>{3, 3, 2, 2, 3, 4, 4, 2, 1, 4, 5, 2, 1, 1, 5, 0},
                       "max-tree native smallest node buffer");
    requireVectorEqual(maxTree->nodeAltitudes(), std::vector<std::uint8_t>{0, 1, 2, 3, 4, 5}, "max-tree native altitude buffer");
    requireEqual(static_cast<int>(maxTree->topology().leaves().size()), 1, "max-tree leaves");
    requireVectorEqual(maxTree->topology().leaves(), std::vector<NodeId>{5}, "max-tree leaves should be dense NodeIds");

    require(minTree->topology().sharedAdjacencyContext() != nullptr, "min-tree must retain its shared-adjacency context");
    requireEqual(minTree->topology().numRows(), 4, "min-tree image rows");
    requireEqual(minTree->topology().numColumns(), 4, "min-tree image columns");
    require(minTree->topology().kind() == MorphologicalTreeKind::MinTree, "tree must be a min-tree");
    require(minTree->topology().nodeAltitudeOrder() == NodeAltitudeOrder::Decreasing, "min-tree altitude order capability");
    requireEqual(minTree->topology().numNodes(), 6, "min-tree node count");
    requireEqual(minTree->topology().getMutationVersion(), std::size_t{7}, "min-tree construction mutation version");
    requireVectorEqual(collectParents(minTree->topology()), std::vector<NodeId>{0, 0, 1, 2, 2, 3}, "min-tree native parent buffer");
    requireVectorEqual(collectOwners(minTree->topology()), std::vector<NodeId>{2, 2, 3, 3, 2, 1, 1, 3, 4, 1, 0, 3, 4, 4, 0, 5},
                       "min-tree native smallest node buffer");
    requireVectorEqual(minTree->nodeAltitudes(), std::vector<std::uint8_t>{5, 4, 3, 2, 1, 0}, "min-tree native altitude buffer");
    requireEqual(static_cast<int>(minTree->topology().leaves().size()), 2, "min-tree leaves");
    requireVectorEqual(minTree->topology().leaves(), std::vector<NodeId>({4, 5}), "min-tree leaves should be dense NodeIds");

    requireEqual(maxTree->topology().root(), 0, "max-tree root dense node id");
    requireEqual(maxTree->nodeAltitude(maxTree->topology().root()), 0, "max-tree root altitude");
    requireEqual(computeArea(maxTree->topology(), maxTree->topology().root()), 16, "max-tree root area");
    require(maxTreeShortcut.topology().kind() == MorphologicalTreeKind::MaxTree, "createMaxTree shortcut must create a max-tree");
    requireEqual(maxTreeShortcut.topology().root(), maxTree->topology().root(), "createMaxTree shortcut root");
    requireEqual(computeArea(maxTreeShortcut.topology(), maxTreeShortcut.topology().root()), 16, "createMaxTree shortcut root area");

    requireEqual(minTree->topology().root(), 0, "min-tree root dense node id");
    requireEqual(minTree->nodeAltitude(minTree->topology().root()), 5, "min-tree root altitude");
    requireEqual(computeArea(minTree->topology(), minTree->topology().root()), 16, "min-tree root area");
    require(minTreeShortcut.topology().kind() == MorphologicalTreeKind::MinTree, "createMinTree shortcut must create a min-tree");
    requireEqual(minTreeShortcut.topology().root(), minTree->topology().root(), "createMinTree shortcut root");
    requireEqual(computeArea(minTreeShortcut.topology(), minTreeShortcut.topology().root()), 16, "createMinTree shortcut root area");

    auto maxReconstruction = maxTree->reconstructFromNodeAltitudes();
    auto minReconstruction = minTree->reconstructFromNodeAltitudes();
    requireImageShape(maxReconstruction, 4, 4);
    requireImageShape(minReconstruction, 4, 4);
    requireVectorEqual(collectImageValues(maxReconstruction), collectImageValues(image), "max-tree reconstruction values");
    requireVectorEqual(collectImageValues(minReconstruction), collectImageValues(image), "min-tree reconstruction values");

    auto importedFromHigra = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(higraParent), std::span<const std::uint8_t>(higraAltitude),
                                                                             4, 4, MorphologicalTreeKind::MaxTree, RegularGridAdjacency2D(4, 4, 1.5));
    requireEqual(importedFromHigra.topology().getNumHigraNodes(), 22, "Higra import node-id domain size");
    requireEqual(importedFromHigra.topology().getMutationVersion(), maxTree->topology().getMutationVersion(), "Higra import construction mutation version");
    requireEqual(importedFromHigra.topology().getHigraNodeId(3), 19, "Higra import slot->node-id mapping");
    requireEqual(computeArea(importedFromHigra.topology(), importedFromHigra.topology().root()), 16, "Higra import area");
    requireImageShape(importedFromHigra.reconstructFromNodeAltitudes(), 4, 4);
    requireVectorEqual(collectImageValues(importedFromHigra.reconstructFromNodeAltitudes()), collectImageValues(image), "Higra import reconstruction values");
    auto [reexportedHigraParent, reexportedHigraAltitude] = importedFromHigra.exportHigraHierarchy();
    requireVectorEqual(reexportedHigraParent, higraParent, "Higra export/import parent round-trip");
    requireVectorEqual(reexportedHigraAltitude, higraAltitude, "Higra export/import altitude round-trip");

    auto imageBuiltTopology = MorphologicalTreeFactory::createMaxTree(image);
    requireThrows([&]() { imageBuiltTopology.topology().getNumHigraNodes(); }, "image-built topology must not expose a Higra node-id domain");

    auto mutatedImport = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(higraParent), std::span<const std::uint8_t>(higraAltitude), 4,
                                                                         4, MorphologicalTreeKind::MaxTree, RegularGridAdjacency2D(4, 4, 1.5));
    mutatedImport.mergeNodeIntoParent(0);
    requireThrows([&]() { mutatedImport.topology().getNumHigraNodes(); }, "mutated Higra import must drop the original Higra node-id domain");
    requireEqual(mutatedImport.topology().getHigraNodeId(3), InvalidNode, "mutated Higra import must not expose stale slot->node-id mapping");

    {
        auto sparseTree = makeValuedComponentTree(image, true);
        sparseTree->mergeNodeIntoParent(5);
        auto [sparseHigraParent, sparseHigraAltitude] = sparseTree->exportHigraHierarchy();
        requireEqual(static_cast<int>(sparseHigraParent.size()), 16 + 5, "Higra export must compact dead node slots");
        auto sparseImported =
            MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(sparseHigraParent), std::span<const std::uint8_t>(sparseHigraAltitude), 4,
                                                            4, MorphologicalTreeKind::MaxTree, RegularGridAdjacency2D(4, 4, 1.5));
        requireEqual(sparseImported.topology().numNodes(), 5, "compact Higra round-trip node count");
        requireEqual(sparseImported.topology().numInternalNodeSlots(), 5, "compact Higra round-trip slot count");
        requireEqual(computeArea(sparseImported.topology(), sparseImported.topology().root()), 16, "compact Higra round-trip area");
        auto [sparseReexportedParent, sparseReexportedAltitude] = sparseImported.exportHigraHierarchy();
        requireVectorEqual(sparseReexportedParent, sparseHigraParent, "compact Higra parent round-trip");
        requireVectorEqual(sparseReexportedAltitude, sparseHigraAltitude, "compact Higra altitude round-trip");
    }

    return 0;
}
