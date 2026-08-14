#include "support/TestSupport.hpp"
#include "mmcfilters/trees/TreeEditor.hpp"

#include <memory>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

int main() {
    {
        auto image = makeComponentTreeFixture();
        auto tree = makeComponentTree(image, true);

        if constexpr (contract::validationsEnabled) {
            requireThrows<std::invalid_argument>([&]() { tree->mergeNodeIntoParent(InvalidNode); }, "mergeNodeIntoParent must reject invalid NodeId");
            requireThrows<std::invalid_argument>([&]() { tree->mergeNodeIntoParent(tree->root()); }, "mergeNodeIntoParent must reject root");
            requireThrows<std::invalid_argument>([&]() { tree->pruneNode(InvalidNode); }, "pruneNode must reject invalid NodeId");
            requireThrows<std::invalid_argument>([&]() { tree->pruneNode(tree->root()); }, "pruneNode must reject root");
        }

        const std::size_t versionBeforeMerge = tree->getMutationVersion();
        tree->mergeNodeIntoParent(5);
        require(tree->getMutationVersion() > versionBeforeMerge, "mergeNodeIntoParent must advance the mutation version");

        requireVectorEqual(collectNodeIds(tree->aliveNodeIds()), {0, 1, 2, 3, 4}, "alive ids after merge");
        require(!tree->isAlive(5), "node 5 must be released after merge");
        requireEqual(tree->getNumFreeNodeSlots(), 1, "free slots after merge");
        requireVectorEqual(collectNodeIds(tree->children(4)), {}, "node 4 children after merge");
        requireVectorEqual(collectPixelIds(tree->properPart(4)), {5, 6, 9, 10, 14}, "node 4 proper parts after merge");
        if constexpr (contract::validationsEnabled) {
            requireThrows<std::invalid_argument>([&]() { tree->mergeNodeIntoParent(5); }, "mergeNodeIntoParent must reject dead slot");
            requireThrows<std::invalid_argument>([&]() { tree->pruneNode(5); }, "pruneNode must reject dead slot");
        }

        auto editor = tree->edit();
        const NodeId reusedNode = editor.createDetachedNode();
        requireEqual(reusedNode, 5, "createDetachedNode must reuse released slot");
        require(tree->isAlive(reusedNode), "reused node must be alive");
        requireEqual(tree->parent(reusedNode), reusedNode, "reused node must start detached");
        requireVectorEqual(collectPixelIds(tree->properPart(reusedNode)), {}, "reused node must start empty");

        editor.attach(4, reusedNode);
        requireEqual(tree->parent(reusedNode), 4, "attached node parent");
        requireVectorEqual(collectNodeIds(tree->children(4)), {5}, "children after attach");

        editor.movePixelToProperPart(5, 4, 5);
        requireEqual(tree->smallestNode(5), 5, "pixel 5 smallest node after movePixelToProperPart");
        requireVectorEqual(collectPixelIds(tree->properPart(4)), {6, 9, 10, 14}, "source proper parts after movePixelToProperPart");
        requireVectorEqual(collectPixelIds(tree->properPart(5)), {5}, "target proper parts after movePixelToProperPart");
        requireEqual(static_cast<int>(collectNodeIds(tree->nodeSupport(4)).size()), 5, "source subtree area after movePixelToProperPart");
        requireEqual(static_cast<int>(collectNodeIds(tree->nodeSupport(5)).size()), 1, "target subtree area after movePixelToProperPart");

        editor.mergeProperParts(5, 4);
        requireVectorEqual(collectPixelIds(tree->properPart(4)), {}, "source proper parts after mergeProperParts");
        requireVectorEqual(collectPixelIds(tree->properPart(5)), {5, 6, 9, 10, 14}, "target proper parts after mergeProperParts");
        requireEqual(static_cast<int>(collectNodeIds(tree->nodeSupport(4)).size()), 5, "source subtree area after mergeProperParts");
        requireEqual(static_cast<int>(collectNodeIds(tree->nodeSupport(5)).size()), 5, "target subtree area after mergeProperParts");

        editor.detach(5);
        requireEqual(tree->parent(5), 5, "detached node parent");
        requireVectorEqual(collectNodeIds(tree->children(4)), {}, "children after detach");
        requireEqual(static_cast<int>(collectNodeIds(tree->nodeSupport(4)).size()), 0, "source subtree area after detach");

        editor.attach(3, 5);
        requireEqual(tree->parent(5), 3, "re-attached node parent");
        editor.reparent(5, 2);
        requireEqual(tree->parent(5), 2, "node parent after moveNode");
        requireVectorEqual(collectNodeIds(tree->children(3)), {4}, "node 3 children after moveNode");
        requireVectorEqual(collectNodeIds(tree->children(2)), {3, 5}, "node 2 children after moveNode");

        editor.setRoot(2);
        requireEqual(tree->root(), 2, "new root after setRoot");
        requireEqual(tree->parent(2), 2, "new root parent must point to itself");
        requireEqual(tree->parent(0), 0, "previous root must become detached after setRoot");
    }

    {
        auto image = makeComponentTreeFixture();
        auto tree = makeComponentTree(image, true);
        auto editor = tree->edit();

        editor.movePixelToProperPart(5, 4, 5);
        editor.mergeProperParts(5, 4);
        editor.commit();
        auto parent = exportFlatHigraHierarchy(*tree).first;
        auto rebuilt = makeTreeFromHigraParent(parent, 4, 4, true);

        requireVectorEqual(collectNodeIds(rebuilt->aliveNodeIds()), {0, 1, 2, 3, 4, 5}, "alive ids after rebuilding ownership-moved tree");
        requireEqual(rebuilt->parent(rebuilt->root()), rebuilt->root(), "rebuilt ownership-moved tree root parent must point to itself");
        const NodeId targetNodeId = rebuilt->smallestNode(10);
        const NodeId sourceNodeId = rebuilt->parent(targetNodeId);
        require(sourceNodeId != InvalidNode && sourceNodeId != targetNodeId, "rebuilt ownership-moved tree must keep the source above the target");
        requireVectorEqual(collectNodeIds(rebuilt->children(sourceNodeId)), {targetNodeId}, "children after rebuilding ownership-moved tree");
        requireVectorEqual(collectPixelIds(rebuilt->properPart(sourceNodeId)), {}, "source proper parts after rebuilding ownership-moved tree");
        requireVectorEqual(collectPixelIds(rebuilt->properPart(targetNodeId)), {5, 6, 9, 10, 14},
                           "target proper parts after rebuilding ownership-moved tree");
        requireEqual(rebuilt->smallestNode(5), targetNodeId, "moved proper part must share the rebuilt target component");
        requireEqual(computeAreaViaAttributeFacade(*rebuilt, sourceNodeId), 5, "source subtree area after rebuilding ownership-moved tree");
        requireEqual(computeAreaViaAttributeFacade(*rebuilt, targetNodeId), 5, "target subtree area after rebuilding ownership-moved tree");
    }

    {
        auto image = makeComponentTreeFixture();
        auto tree = makeComponentTree(image, false);
        auto editor = tree->edit();

        editor.moveChildren(4, 3);
        requireVectorEqual(collectNodeIds(tree->children(3)), {}, "source children after moveChildren");
        requireVectorEqual(collectNodeIds(tree->children(4)), {5}, "target children after moveChildren");
        requireEqual(tree->parent(5), 4, "moved child parent");
        editor.commit();
    }

    {
        auto image = makeComponentTreeFixture();
        auto tree = makeComponentTree(image, true);
        auto editor = tree->edit();

        editor.detach(5);
        requireEqual(tree->parent(5), 5, "detached node parent before rejected merge");
        requireThrows<std::invalid_argument>([&]() { editor.mergeNodeIntoParent(5); }, "mergeNodeIntoParent must reject detached nodes");
        require(tree->isAlive(5), "rejected detached merge must keep node alive");
        requireEqual(tree->parent(5), 5, "rejected detached merge must keep node detached");
    }

    {
        auto image = makeComponentTreeFixture();
        auto tree = makeComponentTree(image, true);
        auto editor = tree->edit();

        editor.detach(5);
        requireEqual(tree->parent(5), 5, "detached node parent before rejected prune");
        requireThrows<std::invalid_argument>([&]() { editor.pruneNode(5); }, "pruneNode must reject detached nodes");
        require(tree->isAlive(5), "rejected detached prune must keep node alive");
        requireEqual(tree->parent(5), 5, "rejected detached prune must keep node detached");
    }

    {
        auto image = makeComponentTreeFixture();
        auto tree = makeComponentTree(image, false);

        requireVectorEqual(collectNodeIds(tree->children(2)), {3, 4}, "fixture children before positional merge");
        requireVectorEqual(collectNodeIds(tree->children(3)), {5}, "fixture merged-node children before positional merge");

        tree->mergeNodeIntoParent(3);
        require(!tree->isAlive(3), "merged middle child slot must be released");
        requireEqual(tree->parent(5), 2, "merged child parent after positional merge");
        requireVectorEqual(collectNodeIds(tree->children(2)), {5, 4}, "mergeNodeIntoParent must preserve sibling position");
    }

    {
        auto image = makeComponentTreeFixture();
        auto tree = makeComponentTree(image, true);

        const std::size_t versionBeforePrune = tree->getMutationVersion();
        tree->pruneNode(4);
        require(tree->getMutationVersion() > versionBeforePrune, "pruneNode must advance the mutation version");

        requireVectorEqual(collectNodeIds(tree->aliveNodeIds()), {0, 1, 2, 3}, "alive ids after prune");
        require(!tree->isAlive(4), "node 4 must be released after prune");
        require(!tree->isAlive(5), "node 5 must be released after prune");
        requireEqual(tree->getNumFreeNodeSlots(), 2, "free slots after prune");
        requireVectorEqual(collectNodeIds(tree->children(3)), {}, "node 3 children after prune");
        requireVectorEqual(collectPixelIds(tree->properPart(3)), {0, 1, 4, 10, 14, 5, 6, 9}, "node 3 proper parts after prune");
        requireEqual(computeAreaViaAttributeFacade(*tree, 3), 8, "pruned node area via attribute computer");
    }

    return 0;
}
