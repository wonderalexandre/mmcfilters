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
        const int expectedInsertedArea = computeAreaViaAttributeFacade(*tree, 3) + computeAreaViaAttributeFacade(*tree, 4);
        const std::size_t versionBeforeEdit = tree->getMutationVersion();
        auto editor = tree->edit();

        requireEqual(tree->getNumInternalNodeSlots(), 6, "initial internal slot count");
        requireEqual(tree->getNumFreeNodeSlots(), 0, "initial free slot count");

        const NodeId insertedNode = editor.createDetachedNode();
        require(tree->isEditing(), "tree must report an open edit session");

        requireEqual(insertedNode, 6, "createDetachedNode must append a fresh slot when none is free");
        require(tree->isAlive(insertedNode), "inserted node must be alive");
        requireEqual(tree->getNodeParent(insertedNode), insertedNode, "inserted node must start detached");
        require(editor.hasDetachedAliveNodes(), "editor must report detached nodes during the staged edit");

        editor.reparent(3, insertedNode);
        editor.reparent(4, insertedNode);
        requireVectorEqual(collectNodeIds(tree->getChildren(insertedNode)), {3, 4}, "new node children before attach");

        editor.attach(2, insertedNode);
        require(!editor.hasDetachedAliveNodes(), "editor must stop reporting detached nodes after re-attachment");

        const TreeValidationResult commitResult = editor.validateAndCommit();
        require(commitResult.ok, "validateAndCommit must accept a repaired staged insertion");
        require(!tree->isEditing(), "successful validateAndCommit must close the edit session");
        require(tree->getMutationVersion() > versionBeforeEdit, "committed staged edit must advance the mutation version");

        requireEqual(tree->getNodeParent(insertedNode), 2, "inserted node parent after commit");
        requireVectorEqual(collectNodeIds(tree->getChildren(2)), {6}, "parent children after commit");
        requireVectorEqual(collectNodeIds(tree->getChildren(insertedNode)), {3, 4}, "inserted node children after commit");
        requireVectorEqual(collectNodeIds(tree->getPathToRootNodes(5)), {5, 3, 6, 2, 1, 0}, "path to root after staged insertion");
        requireEqual(computeAreaViaAttributeFacade(*tree, insertedNode), expectedInsertedArea, "inserted node area after commit");

        const auto exportedHigra = exportFlatHigraHierarchy(*tree);
        requireEqual(static_cast<int>(exportedHigra.first.size()), tree->getNumTotalProperParts() + tree->getNumNodes(), "exported Higra parent size after commit");
    }

    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), true);
        auto editor = tree->edit();

        const NodeId detachedNode = editor.createDetachedNode();
        require(editor.hasDetachedAliveNodes(), "detached edit must be visible before commit");
        requireThrows([&]() { (void)exportFlatHigraHierarchy(*tree); }, "Higra export must reject detached alive nodes");
        requireThrows([&]() { (void)computeAreaViaAttributeFacade(*tree, tree->getRoot()); }, "attribute computation must reject an open edit session");

        bool threw = false;
        try {
            editor.commit();
        } catch (const std::runtime_error&) {
            threw = true;
        }
        require(threw, "commit must reject a tree that still has detached alive nodes");
        require(tree->isEditing(), "failed commit must leave the edit session open");

        editor.attach(4, detachedNode);
        const TreeValidationResult repairedCommit = editor.validateAndCommit();
        require(repairedCommit.ok, "validateAndCommit must succeed after repairing the edit");
        requireEqual(tree->getNodeParent(detachedNode), 4, "editor must remain usable after a failed commit");
        require(!editor.hasDetachedAliveNodes(), "successful re-attachment must clear the detached state");
        require(!tree->isEditing(), "successful validateAndCommit must close the repaired session");
    }

    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), true);
        auto editor = tree->edit();

        editor.reparent(3, 4);
        requireThrows([&]() { editor.commit(); }, "commit must reject a cycle introduced by staged reparent");
        require(tree->isEditing(), "failed cyclic commit must leave the session open");
    }

    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), true);
        auto editor = tree->edit();

        editor.detach(3);
        editor.attach(4, 3);
        requireThrows([&]() { editor.commit(); }, "commit must reject a cycle introduced by staged attach");
        require(tree->isEditing(), "failed cyclic attach commit must leave the session open");
    }

    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), true);
        auto editor = tree->edit();
        requireThrows([&]() { (void)tree->edit(); }, "nested edit sessions must be rejected");
        editor.commitUnchecked();
        require(!tree->isEditing(), "commitUnchecked must close an internal hot-path edit session");
    }

    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), true);
        {
            auto editor = tree->edit();
            static_cast<void>(editor.createDetachedNode());
            require(tree->isEditing(), "uncommitted editor must keep the tree in editing mode");
        }
        require(tree->isEditing(), "destroying an uncommitted editor must leave the tree in editing mode");
        requireThrows([&]() { (void)tree->edit(); }, "tree with abandoned edit session must reject a new editor");
        requireThrows([&]() { (void)computeAreaViaAttributeFacade(*tree, tree->getRoot()); }, "tree with abandoned edit session must reject committed-tree APIs");
    }

    return 0;
}
