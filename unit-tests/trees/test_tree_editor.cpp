#include "support/TestSupport.hpp"
#include "mmcfilters/trees/TreeEditor.hpp"

#include <memory>
#include <stdexcept>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

int main() {
    auto requireThrows = [](auto&& fn, const std::string& label) {
        bool threw = false;
        try {
            fn();
        } catch (const std::exception&) {
            threw = true;
        }
        require(threw, label);
    };

    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), false);
        TreeEditor editor(*tree);

        requireEqual(tree->getNumInternalNodeSlots(), 6, "initial internal slot count");
        requireEqual(tree->getNumFreeNodeSlots(), 0, "initial free slot count");

        const int expectedInsertedArea = computeAreaAttribute(*tree, 3) + computeAreaAttribute(*tree, 4);
        const NodeId insertedNode = editor.createDetachedNode();

        requireEqual(insertedNode, 6, "createDetachedNode must append a fresh slot when none is free");
        require(tree->isAlive(insertedNode), "inserted node must be alive");
        requireEqual(tree->getNodeParent(insertedNode), insertedNode, "inserted node must start detached");
        require(editor.hasDetachedAliveNodes(), "editor must report detached nodes during the staged edit");

        editor.reparent(3, insertedNode);
        editor.reparent(4, insertedNode);
        requireVectorEqual(collectNodeIds(tree->getChildren(insertedNode)), {3, 4}, "new node children before attach");

        editor.attach(2, insertedNode);
        require(!editor.hasDetachedAliveNodes(), "editor must stop reporting detached nodes after re-attachment");

        editor.commit();

        requireEqual(tree->getNodeParent(insertedNode), 2, "inserted node parent after commit");
        requireVectorEqual(collectNodeIds(tree->getChildren(2)), {6}, "parent children after commit");
        requireVectorEqual(collectNodeIds(tree->getChildren(insertedNode)), {3, 4}, "inserted node children after commit");
        requireVectorEqual(collectNodeIds(tree->getPathToRootNodes(5)), {5, 3, 6, 2, 1, 0}, "path to root after staged insertion");
        requireEqual(computeAreaAttribute(*tree, insertedNode), expectedInsertedArea, "inserted node area after commit");

        const auto exportedHigra = exportFlatHigraHierarchy(*tree);
        requireEqual(static_cast<int>(exportedHigra.first.size()), tree->getNumTotalProperParts() + tree->getNumNodes(), "exported Higra parent size after commit");
    }

    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), true);
        TreeEditor editor(*tree);

        const NodeId detachedNode = editor.createDetachedNode();
        require(editor.hasDetachedAliveNodes(), "detached edit must be visible before commit");
        requireThrows([&]() { (void)exportFlatHigraHierarchy(*tree); }, "Higra export must reject detached alive nodes");

        bool threw = false;
        try {
            editor.commit();
        } catch (const std::runtime_error&) {
            threw = true;
        }
        require(threw, "commit must reject a tree that still has detached alive nodes");

        editor.attach(4, detachedNode);
        editor.commit();
        requireEqual(tree->getNodeParent(detachedNode), 4, "editor must remain usable after a failed commit");
        require(!editor.hasDetachedAliveNodes(), "successful re-attachment must clear the detached state");
    }

    return 0;
}
