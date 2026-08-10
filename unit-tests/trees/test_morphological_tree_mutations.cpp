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
            requireThrows<std::invalid_argument>([&]() { tree->mergeNodeIntoParent(tree->getRoot()); }, "mergeNodeIntoParent must reject root");
            requireThrows<std::invalid_argument>([&]() { tree->pruneNode(InvalidNode); }, "pruneNode must reject invalid NodeId");
            requireThrows<std::invalid_argument>([&]() { tree->pruneNode(tree->getRoot()); }, "pruneNode must reject root");
        }

        const std::size_t versionBeforeMerge = tree->getMutationVersion();
        tree->mergeNodeIntoParent(5);
        require(tree->getMutationVersion() > versionBeforeMerge, "mergeNodeIntoParent must advance the mutation version");

        requireVectorEqual(collectNodeIds(tree->getAliveNodeIds()), {0, 1, 2, 3, 4}, "alive ids after merge");
        require(!tree->isAlive(5), "node 5 must be released after merge");
        requireEqual(tree->getNumFreeNodeSlots(), 1, "free slots after merge");
        requireVectorEqual(collectNodeIds(tree->getChildren(4)), {}, "node 4 children after merge");
        requireVectorEqual(collectNodeIds(tree->getProperParts(4)), {5, 6, 9, 10, 14}, "node 4 proper parts after merge");
        if constexpr (contract::validationsEnabled) {
            requireThrows<std::invalid_argument>([&]() { tree->mergeNodeIntoParent(5); }, "mergeNodeIntoParent must reject dead slot");
            requireThrows<std::invalid_argument>([&]() { tree->pruneNode(5); }, "pruneNode must reject dead slot");
        }

        auto editor = tree->edit();
        const NodeId reusedNode = editor.createDetachedNode();
        requireEqual(reusedNode, 5, "createDetachedNode must reuse released slot");
        require(tree->isAlive(reusedNode), "reused node must be alive");
        requireEqual(tree->getNodeParent(reusedNode), reusedNode, "reused node must start detached");
        requireVectorEqual(collectNodeIds(tree->getProperParts(reusedNode)), {}, "reused node must start empty");

        editor.attach(4, reusedNode);
        requireEqual(tree->getNodeParent(reusedNode), 4, "attached node parent");
        requireVectorEqual(collectNodeIds(tree->getChildren(4)), {5}, "children after attach");

        editor.moveProperPart(5, 4, 5);
        requireEqual(tree->getProperPartOwner(5), 5, "pixel 5 owner after moveProperPart");
        requireVectorEqual(collectNodeIds(tree->getProperParts(4)), {6, 9, 10, 14}, "source proper parts after moveProperPart");
        requireVectorEqual(collectNodeIds(tree->getProperParts(5)), {5}, "target proper parts after moveProperPart");
        requireEqual(static_cast<int>(collectNodeIds(tree->getConnectedComponent(4)).size()), 5, "source subtree area after moveProperPart");
        requireEqual(static_cast<int>(collectNodeIds(tree->getConnectedComponent(5)).size()), 1, "target subtree area after moveProperPart");

        editor.moveProperParts(5, 4);
        requireVectorEqual(collectNodeIds(tree->getProperParts(4)), {}, "source proper parts after moveProperParts");
        requireVectorEqual(collectNodeIds(tree->getProperParts(5)), {5, 6, 9, 10, 14}, "target proper parts after moveProperParts");
        requireEqual(static_cast<int>(collectNodeIds(tree->getConnectedComponent(4)).size()), 5, "source subtree area after moveProperParts");
        requireEqual(static_cast<int>(collectNodeIds(tree->getConnectedComponent(5)).size()), 5, "target subtree area after moveProperParts");

        editor.detach(5);
        requireEqual(tree->getNodeParent(5), 5, "detached node parent");
        requireVectorEqual(collectNodeIds(tree->getChildren(4)), {}, "children after detach");
        requireEqual(static_cast<int>(collectNodeIds(tree->getConnectedComponent(4)).size()), 0, "source subtree area after detach");

        editor.attach(3, 5);
        requireEqual(tree->getNodeParent(5), 3, "re-attached node parent");
        editor.reparent(5, 2);
        requireEqual(tree->getNodeParent(5), 2, "node parent after moveNode");
        requireVectorEqual(collectNodeIds(tree->getChildren(3)), {4}, "node 3 children after moveNode");
        requireVectorEqual(collectNodeIds(tree->getChildren(2)), {3, 5}, "node 2 children after moveNode");

        editor.setRoot(2);
        requireEqual(tree->getRoot(), 2, "new root after setRoot");
        requireEqual(tree->getNodeParent(2), 2, "new root parent must point to itself");
        requireEqual(tree->getNodeParent(0), 0, "previous root must become detached after setRoot");
    }

    {
        auto image = makeComponentTreeFixture();
        auto tree = makeComponentTree(image, true);
        auto editor = tree->edit();

        editor.moveProperPart(5, 4, 5);
        editor.moveProperParts(5, 4);
        editor.commit();
        auto parent = exportFlatHigraHierarchy(*tree).first;
        auto rebuilt = makeTreeFromHigraParent(parent, 4, 4, true);

        requireVectorEqual(collectNodeIds(rebuilt->getAliveNodeIds()), {0, 1, 2, 3, 4, 5}, "alive ids after rebuilding ownership-moved tree");
        requireEqual(rebuilt->getNodeParent(rebuilt->getRoot()), rebuilt->getRoot(), "rebuilt ownership-moved tree root parent must point to itself");
        const NodeId targetNodeId = rebuilt->getProperPartOwner(10);
        const NodeId sourceNodeId = rebuilt->getNodeParent(targetNodeId);
        require(sourceNodeId != InvalidNode && sourceNodeId != targetNodeId, "rebuilt ownership-moved tree must keep the source above the target");
        requireVectorEqual(collectNodeIds(rebuilt->getChildren(sourceNodeId)), {targetNodeId}, "children after rebuilding ownership-moved tree");
        requireVectorEqual(collectNodeIds(rebuilt->getProperParts(sourceNodeId)), {}, "source proper parts after rebuilding ownership-moved tree");
        requireVectorEqual(collectNodeIds(rebuilt->getProperParts(targetNodeId)), {5, 6, 9, 10, 14},
                           "target proper parts after rebuilding ownership-moved tree");
        requireEqual(rebuilt->getProperPartOwner(5), targetNodeId, "moved proper part must share the rebuilt target component");
        requireEqual(computeAreaViaAttributeFacade(*rebuilt, sourceNodeId), 5, "source subtree area after rebuilding ownership-moved tree");
        requireEqual(computeAreaViaAttributeFacade(*rebuilt, targetNodeId), 5, "target subtree area after rebuilding ownership-moved tree");
    }

    {
        auto image = makeComponentTreeFixture();
        auto tree = makeComponentTree(image, false);
        auto editor = tree->edit();

        editor.moveChildren(4, 3);
        requireVectorEqual(collectNodeIds(tree->getChildren(3)), {}, "source children after moveChildren");
        requireVectorEqual(collectNodeIds(tree->getChildren(4)), {5}, "target children after moveChildren");
        requireEqual(tree->getNodeParent(5), 4, "moved child parent");
        editor.commit();
    }

    {
        auto image = makeComponentTreeFixture();
        auto tree = makeComponentTree(image, true);
        auto editor = tree->edit();

        editor.detach(5);
        requireEqual(tree->getNodeParent(5), 5, "detached node parent before rejected merge");
        requireThrows<std::invalid_argument>([&]() { editor.mergeNodeIntoParent(5); }, "mergeNodeIntoParent must reject detached nodes");
        require(tree->isAlive(5), "rejected detached merge must keep node alive");
        requireEqual(tree->getNodeParent(5), 5, "rejected detached merge must keep node detached");
    }

    {
        auto image = makeComponentTreeFixture();
        auto tree = makeComponentTree(image, true);
        auto editor = tree->edit();

        editor.detach(5);
        requireEqual(tree->getNodeParent(5), 5, "detached node parent before rejected prune");
        requireThrows<std::invalid_argument>([&]() { editor.pruneNode(5); }, "pruneNode must reject detached nodes");
        require(tree->isAlive(5), "rejected detached prune must keep node alive");
        requireEqual(tree->getNodeParent(5), 5, "rejected detached prune must keep node detached");
    }

    {
        auto image = makeComponentTreeFixture();
        auto tree = makeComponentTree(image, false);

        requireVectorEqual(collectNodeIds(tree->getChildren(2)), {3, 4}, "fixture children before positional merge");
        requireVectorEqual(collectNodeIds(tree->getChildren(3)), {5}, "fixture merged-node children before positional merge");

        tree->mergeNodeIntoParent(3);
        require(!tree->isAlive(3), "merged middle child slot must be released");
        requireEqual(tree->getNodeParent(5), 2, "merged child parent after positional merge");
        requireVectorEqual(collectNodeIds(tree->getChildren(2)), {5, 4}, "mergeNodeIntoParent must preserve sibling position");
    }

    {
        auto image = makeComponentTreeFixture();
        auto tree = makeComponentTree(image, true);

        const std::size_t versionBeforePrune = tree->getMutationVersion();
        tree->pruneNode(4);
        require(tree->getMutationVersion() > versionBeforePrune, "pruneNode must advance the mutation version");

        requireVectorEqual(collectNodeIds(tree->getAliveNodeIds()), {0, 1, 2, 3}, "alive ids after prune");
        require(!tree->isAlive(4), "node 4 must be released after prune");
        require(!tree->isAlive(5), "node 5 must be released after prune");
        requireEqual(tree->getNumFreeNodeSlots(), 2, "free slots after prune");
        requireVectorEqual(collectNodeIds(tree->getChildren(3)), {}, "node 3 children after prune");
        requireVectorEqual(collectNodeIds(tree->getProperParts(3)), {0, 1, 4, 10, 14, 5, 6, 9}, "node 3 proper parts after prune");
        requireEqual(computeAreaViaAttributeFacade(*tree, 3), 8, "pruned node area via attribute computer");
    }

    return 0;
}
