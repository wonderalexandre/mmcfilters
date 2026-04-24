#include "support/TestSupport.hpp"
#include "mmcfilters/trees/TreeEditor.hpp"

#include <memory>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

int main() {
    {
        auto image = makeComponentTreeFixture();
        auto tree = std::make_shared<MorphologicalTree>(image, true);
        TreeEditor editor(*tree);

        tree->mergeNodeIntoParent(5);

        requireVectorEqual(collectNodeIds(tree->getAliveNodeIds()), {0, 1, 2, 3, 4}, "alive ids after merge");
        require(!tree->isAlive(5), "node 5 must be released after merge");
        requireEqual(tree->getNumFreeNodeSlots(), 1, "free slots after merge");
        requireVectorEqual(collectNodeIds(tree->getChildren(4)), {}, "node 4 children after merge");
        requireVectorEqual(collectNodeIds(tree->getProperParts(4)), {5, 6, 9, 10, 14}, "node 4 proper parts after merge");

        const NodeId reusedNode = editor.createDetachedNode();
        requireEqual(reusedNode, 5, "createDetachedNode must reuse released slot");
        require(tree->isAlive(reusedNode), "reused node must be alive");
        requireEqual(tree->getNodeParent(reusedNode), reusedNode, "reused node must start detached");
        requireVectorEqual(collectNodeIds(tree->getProperParts(reusedNode)), {}, "reused node must start empty");

        editor.attach(4, reusedNode);
        requireEqual(tree->getNodeParent(reusedNode), 4, "attached node parent");
        requireVectorEqual(collectNodeIds(tree->getChildren(4)), {5}, "children after attach");

        editor.moveProperPart(5, 4, 5);
        requireEqual(tree->getSmallestComponent(5), 5, "pixel 5 owner after moveProperPart");
        requireVectorEqual(collectNodeIds(tree->getProperParts(4)), {6, 9, 10, 14}, "source proper parts after moveProperPart");
        requireVectorEqual(collectNodeIds(tree->getProperParts(5)), {5}, "target proper parts after moveProperPart");
        requireEqual(computeAreaAttribute(*tree, 4), 5, "source subtree area after moveProperPart");
        requireEqual(computeAreaAttribute(*tree, 5), 1, "target subtree area after moveProperPart");

        editor.moveProperParts(5, 4);
        requireVectorEqual(collectNodeIds(tree->getProperParts(4)), {}, "source proper parts after moveProperParts");
        requireVectorEqual(collectNodeIds(tree->getProperParts(5)), {5, 6, 9, 10, 14}, "target proper parts after moveProperParts");
        requireEqual(computeAreaAttribute(*tree, 4), 5, "source subtree area after moveProperParts");
        requireEqual(computeAreaAttribute(*tree, 5), 5, "target subtree area after moveProperParts");

        editor.detach(5);
        requireEqual(tree->getNodeParent(5), 5, "detached node parent");
        requireVectorEqual(collectNodeIds(tree->getChildren(4)), {}, "children after detach");
        requireEqual(computeAreaAttribute(*tree, 4), 0, "source subtree area after detach");

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
        auto tree = std::make_shared<MorphologicalTree>(image, true);
        TreeEditor editor(*tree);

        editor.moveProperPart(5, 4, 5);
        editor.moveProperParts(5, 4);
        auto parent = exportParentArray(*tree);
        auto rebuilt = std::make_shared<MorphologicalTree>(std::span<const NodeId>(parent), 4, 4, true);

        requireVectorEqual(collectNodeIds(rebuilt->getAliveNodeIds()), {0, 1, 2, 3, 4, 5}, "alive ids after rebuilding ownership-moved tree");
        requireEqual(rebuilt->getNodeParent(rebuilt->getRoot()), rebuilt->getRoot(), "rebuilt ownership-moved tree root parent must point to itself");
        requireVectorEqual(collectNodeIds(rebuilt->getChildren(4)), {5}, "children after rebuilding ownership-moved tree");
        requireVectorEqual(collectNodeIds(rebuilt->getProperParts(4)), {}, "source proper parts after rebuilding ownership-moved tree");
        requireVectorEqual(collectNodeIds(rebuilt->getProperParts(5)), {5, 6, 9, 10, 14}, "target proper parts after rebuilding ownership-moved tree");
        requireEqual(rebuilt->getSmallestComponent(10), 5, "smallest component after rebuilding merged tree");
        requireEqual(computeAreaAttribute(*rebuilt, 4), 5, "source subtree area after rebuilding ownership-moved tree");
        requireEqual(computeAreaAttribute(*rebuilt, 5), 5, "target subtree area after rebuilding ownership-moved tree");
    }

    {
        auto image = makeComponentTreeFixture();
        auto tree = std::make_shared<MorphologicalTree>(image, false);
        TreeEditor editor(*tree);

        editor.moveChildren(4, 3);
        requireVectorEqual(collectNodeIds(tree->getChildren(3)), {}, "source children after moveChildren");
        requireVectorEqual(collectNodeIds(tree->getChildren(4)), {5}, "target children after moveChildren");
        requireEqual(tree->getNodeParent(5), 4, "moved child parent");
    }

    {
        auto image = makeComponentTreeFixture();
        auto tree = std::make_shared<MorphologicalTree>(image, true);

        tree->pruneNode(4);

        requireVectorEqual(collectNodeIds(tree->getAliveNodeIds()), {0, 1, 2, 3}, "alive ids after prune");
        require(!tree->isAlive(4), "node 4 must be released after prune");
        require(!tree->isAlive(5), "node 5 must be released after prune");
        requireEqual(tree->getNumFreeNodeSlots(), 2, "free slots after prune");
        requireVectorEqual(collectNodeIds(tree->getChildren(3)), {}, "node 3 children after prune");
        requireVectorEqual(collectNodeIds(tree->getProperParts(3)), {0, 1, 4, 5, 6, 9, 10, 14}, "node 3 proper parts after prune");
        requireEqual(computeAreaAttribute(*tree, 3), 8, "pruned node area via attribute computer");
    }

    return 0;
}
