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

    auto maxTree = makeWeightedComponentTree(image, true);
    auto minTree = makeWeightedComponentTree(image, false);
    auto [higraParent, higraAltitude] = exportHigraHierarchy(*maxTree);
    auto maxTreeShortcut = MorphologicalTreeFactory::createMaxTree(image);
    auto minTreeShortcut = MorphologicalTreeFactory::createMinTree(image);
    const auto collectParents = [](const MorphologicalTree& tree) {
        std::vector<NodeId> parent;
        parent.reserve(static_cast<std::size_t>(tree.getNumInternalNodeSlots()));
        for (NodeId node = 0; node < tree.getNumInternalNodeSlots(); ++node) {
            parent.push_back(tree.getNodeParent(node));
        }
        return parent;
    };
    const auto collectOwners = [](const MorphologicalTree& tree) {
        std::vector<NodeId> owner;
        owner.reserve(static_cast<std::size_t>(tree.getNumTotalProperParts()));
        for (NodeId properPart = 0; properPart < tree.getNumTotalProperParts(); ++properPart) {
            owner.push_back(tree.getProperPartOwner(properPart));
        }
        return owner;
    };

    require(static_cast<bool>(maxTree), "max-tree instance must be created");
    require(static_cast<bool>(minTree), "min-tree instance must be created");

    require(maxTree->topology().hasUniformGridAdjacency2D(), "max-tree must expose adjacency relation");
    requireEqual(maxTree->topology().getNumRowsOfGridDomain2D(), 4, "max-tree image rows");
    requireEqual(maxTree->topology().getNumColsOfGridDomain2D(), 4, "max-tree image cols");
    require(maxTree->topology().getDescriptiveKind() == MorphologicalTreeKind::MAX_TREE, "tree must be a max-tree");
    require(maxTree->topology().getAltitudeOrder() == AltitudeOrder::INCREASING_FROM_ROOT, "max-tree altitude order capability");
    require(maxTree->topology().getAdjacencyMode() == AdjacencyMode::UNIFORM, "max-tree uniform adjacency capability");
    require(maxTree->topology().getUniformGridAdjacency2D() != nullptr, "max-tree uniform adjacency relation");
    requireEqual(maxTree->topology().getNumNodes(), 6, "max-tree node count");
    requireEqual(maxTree->topology().getMutationVersion(), std::size_t{7}, "max-tree construction mutation version");
    requireVectorEqual(collectParents(maxTree->topology()), std::vector<NodeId>{0, 0, 1, 2, 3, 4}, "max-tree native parent buffer");
    requireVectorEqual(collectOwners(maxTree->topology()), std::vector<NodeId>{3, 3, 2, 2, 3, 4, 4, 2, 1, 4, 5, 2, 1, 1, 5, 0},
                       "max-tree native proper-part owner buffer");
    requireVectorEqual(maxTree->getAltitudeBuffer(), std::vector<std::uint8_t>{0, 1, 2, 3, 4, 5}, "max-tree native altitude buffer");
    requireEqual(static_cast<int>(maxTree->topology().getLeaves().size()), 1, "max-tree leaves");
    requireVectorEqual(maxTree->topology().getLeaves(), std::vector<NodeId>{5}, "max-tree leaves should be dense NodeIds");

    require(minTree->topology().hasUniformGridAdjacency2D(), "min-tree must expose adjacency relation");
    requireEqual(minTree->topology().getNumRowsOfGridDomain2D(), 4, "min-tree image rows");
    requireEqual(minTree->topology().getNumColsOfGridDomain2D(), 4, "min-tree image cols");
    require(minTree->topology().getDescriptiveKind() == MorphologicalTreeKind::MIN_TREE, "tree must be a min-tree");
    require(minTree->topology().getAltitudeOrder() == AltitudeOrder::DECREASING_FROM_ROOT, "min-tree altitude order capability");
    requireEqual(minTree->topology().getNumNodes(), 6, "min-tree node count");
    requireEqual(minTree->topology().getMutationVersion(), std::size_t{7}, "min-tree construction mutation version");
    requireVectorEqual(collectParents(minTree->topology()), std::vector<NodeId>{0, 0, 1, 2, 2, 3}, "min-tree native parent buffer");
    requireVectorEqual(collectOwners(minTree->topology()), std::vector<NodeId>{2, 2, 3, 3, 2, 1, 1, 3, 4, 1, 0, 3, 4, 4, 0, 5},
                       "min-tree native proper-part owner buffer");
    requireVectorEqual(minTree->getAltitudeBuffer(), std::vector<std::uint8_t>{5, 4, 3, 2, 1, 0}, "min-tree native altitude buffer");
    requireEqual(static_cast<int>(minTree->topology().getLeaves().size()), 2, "min-tree leaves");
    requireVectorEqual(minTree->topology().getLeaves(), std::vector<NodeId>({4, 5}), "min-tree leaves should be dense NodeIds");

    requireEqual(maxTree->topology().getRoot(), 0, "max-tree root dense node id");
    requireEqual(maxTree->getAltitude(maxTree->topology().getRoot()), 0, "max-tree root altitude");
    requireEqual(computeArea(maxTree->topology(), maxTree->topology().getRoot()), 16, "max-tree root area");
    require(maxTreeShortcut.topology().getDescriptiveKind() == MorphologicalTreeKind::MAX_TREE, "createMaxTree shortcut must create a max-tree");
    requireEqual(maxTreeShortcut.topology().getRoot(), maxTree->topology().getRoot(), "createMaxTree shortcut root");
    requireEqual(computeArea(maxTreeShortcut.topology(), maxTreeShortcut.topology().getRoot()), 16, "createMaxTree shortcut root area");

    requireEqual(minTree->topology().getRoot(), 0, "min-tree root dense node id");
    requireEqual(minTree->getAltitude(minTree->topology().getRoot()), 5, "min-tree root altitude");
    requireEqual(computeArea(minTree->topology(), minTree->topology().getRoot()), 16, "min-tree root area");
    require(minTreeShortcut.topology().getDescriptiveKind() == MorphologicalTreeKind::MIN_TREE, "createMinTree shortcut must create a min-tree");
    requireEqual(minTreeShortcut.topology().getRoot(), minTree->topology().getRoot(), "createMinTree shortcut root");
    requireEqual(computeArea(minTreeShortcut.topology(), minTreeShortcut.topology().getRoot()), 16, "createMinTree shortcut root area");

    auto maxReconstruction = maxTree->reconstructionImage();
    auto minReconstruction = minTree->reconstructionImage();
    requireImageShape(maxReconstruction, 4, 4);
    requireImageShape(minReconstruction, 4, 4);
    requireVectorEqual(collectImageValues(maxReconstruction), collectImageValues(image), "max-tree reconstruction values");
    requireVectorEqual(collectImageValues(minReconstruction), collectImageValues(image), "min-tree reconstruction values");

    auto importedFromHigra = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(higraParent), std::span<const std::uint8_t>(higraAltitude),
                                                                             4, 4, MorphologicalTreeKind::MAX_TREE, RegularGridAdjacency2D(4, 4, 1.5));
    requireEqual(importedFromHigra.topology().getNumHigraNodes(), 22, "Higra import node-id domain size");
    requireEqual(importedFromHigra.topology().getMutationVersion(), maxTree->topology().getMutationVersion(), "Higra import construction mutation version");
    requireEqual(importedFromHigra.topology().getHigraNodeId(3), 19, "Higra import slot->node-id mapping");
    requireEqual(computeArea(importedFromHigra.topology(), importedFromHigra.topology().getRoot()), 16, "Higra import area");
    requireImageShape(importedFromHigra.reconstructionImage(), 4, 4);
    requireVectorEqual(collectImageValues(importedFromHigra.reconstructionImage()), collectImageValues(image), "Higra import reconstruction values");
    auto [reexportedHigraParent, reexportedHigraAltitude] = importedFromHigra.exportHigraHierarchy();
    requireVectorEqual(reexportedHigraParent, higraParent, "Higra export/import parent round-trip");
    requireVectorEqual(reexportedHigraAltitude, higraAltitude, "Higra export/import altitude round-trip");

    auto imageBuiltTopology = MorphologicalTreeFactory::createMaxTree(image);
    requireThrows([&]() { imageBuiltTopology.topology().getNumHigraNodes(); }, "image-built topology must not expose a Higra node-id domain");

    auto mutatedImport = MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(higraParent), std::span<const std::uint8_t>(higraAltitude), 4,
                                                                         4, MorphologicalTreeKind::MAX_TREE, RegularGridAdjacency2D(4, 4, 1.5));
    mutatedImport.mergeNodeIntoParent(0);
    requireThrows([&]() { mutatedImport.topology().getNumHigraNodes(); }, "mutated Higra import must drop the original Higra node-id domain");
    requireEqual(mutatedImport.topology().getHigraNodeId(3), InvalidNode, "mutated Higra import must not expose stale slot->node-id mapping");

    {
        auto sparseTree = makeWeightedComponentTree(image, true);
        sparseTree->mergeNodeIntoParent(5);
        auto [sparseHigraParent, sparseHigraAltitude] = sparseTree->exportHigraHierarchy();
        requireEqual(static_cast<int>(sparseHigraParent.size()), 16 + 5, "Higra export must compact dead node slots");
        auto sparseImported =
            MorphologicalTreeFactory::createFromHigraParent(std::span<const NodeId>(sparseHigraParent), std::span<const std::uint8_t>(sparseHigraAltitude), 4,
                                                            4, MorphologicalTreeKind::MAX_TREE, RegularGridAdjacency2D(4, 4, 1.5));
        requireEqual(sparseImported.topology().getNumNodes(), 5, "compact Higra round-trip node count");
        requireEqual(sparseImported.topology().getNumInternalNodeSlots(), 5, "compact Higra round-trip slot count");
        requireEqual(computeArea(sparseImported.topology(), sparseImported.topology().getRoot()), 16, "compact Higra round-trip area");
        auto [sparseReexportedParent, sparseReexportedAltitude] = sparseImported.exportHigraHierarchy();
        requireVectorEqual(sparseReexportedParent, sparseHigraParent, "compact Higra parent round-trip");
        requireVectorEqual(sparseReexportedAltitude, sparseHigraAltitude, "compact Higra altitude round-trip");
    }

    return 0;
}
