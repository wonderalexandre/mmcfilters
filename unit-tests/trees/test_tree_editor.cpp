#include "support/TestSupport.hpp"
#include "mmcfilters/trees/TreeEditor.hpp"

#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

static_assert(std::is_move_constructible_v<TreeEditor::IncrementalProof>);
static_assert(!std::is_copy_constructible_v<TreeEditor::IncrementalProof>);
static_assert(std::is_move_constructible_v<MorphologicalTree>);
static_assert(std::is_move_assignable_v<MorphologicalTree>);
static_assert(!std::is_nothrow_move_constructible_v<MorphologicalTree>);

namespace {

struct PublicTreeState {
    NodeId root = InvalidNode;
    int numNodes = 0;
    int numSlots = 0;
    int numFreeSlots = 0;
    std::size_t mutationVersion = 0;
    std::vector<std::uint8_t> alive;
    std::vector<NodeId> parents;
    std::vector<NodeId> orderedChildren;
    std::vector<NodeId> orderedProperParts;
    std::vector<NodeId> properPartOwners;
};

PublicTreeState capturePublicTreeState(const MorphologicalTree& tree) {
    PublicTreeState state;
    state.root = tree.getRoot();
    state.numNodes = tree.getNumNodes();
    state.numSlots = tree.getNumInternalNodeSlots();
    state.numFreeSlots = tree.getNumFreeNodeSlots();
    state.mutationVersion = tree.getMutationVersion();
    state.alive.reserve(static_cast<std::size_t>(state.numSlots));
    state.parents.reserve(static_cast<std::size_t>(state.numSlots));

    for (NodeId node = 0; node < state.numSlots; ++node) {
        const bool alive = tree.isAlive(node);
        state.alive.push_back(static_cast<std::uint8_t>(alive));
        state.parents.push_back(alive ? tree.getNodeParent(node) : InvalidNode);
        if (alive) {
            for (NodeId child : tree.getChildren(node)) {
                state.orderedChildren.push_back(child);
            }
            for (NodeId properPart : tree.getProperParts(node)) {
                state.orderedProperParts.push_back(properPart);
            }
        }
        state.orderedChildren.push_back(InvalidNode);
        state.orderedProperParts.push_back(InvalidNode);
    }

    for (NodeId properPart = 0; properPart < tree.getNumTotalProperParts(); ++properPart) {
        state.properPartOwners.push_back(tree.getProperPartOwner(properPart));
    }
    return state;
}

void requirePublicTreeStateEqual(const MorphologicalTree& tree, const PublicTreeState& expected, const std::string& label) {
    const PublicTreeState actual = capturePublicTreeState(tree);
    requireEqual(actual.root, expected.root, label + " root");
    requireEqual(actual.numNodes, expected.numNodes, label + " live nodes");
    requireEqual(actual.numSlots, expected.numSlots, label + " node slots");
    requireEqual(actual.numFreeSlots, expected.numFreeSlots, label + " free slots");
    requireEqual(actual.mutationVersion, expected.mutationVersion, label + " mutation version");
    requireVectorEqual(actual.alive, expected.alive, label + " alive slots");
    requireVectorEqual(actual.parents, expected.parents, label + " parents");
    requireVectorEqual(actual.orderedChildren, expected.orderedChildren, label + " ordered children");
    requireVectorEqual(actual.orderedProperParts, expected.orderedProperParts, label + " ordered proper parts");
    requireVectorEqual(actual.properPartOwners, expected.properPartOwners, label + " proper-part owners");
}

} // namespace

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
        requireEqual(static_cast<int>(exportedHigra.first.size()), tree->getNumTotalProperParts() + tree->getNumNodes(),
                     "exported Higra parent size after commit");
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

        const auto donorProperParts = collectNodeIds(tree->getProperParts(4));
        require(!donorProperParts.empty(), "repair fixture must expose a direct proper part");
        editor.moveProperPart(detachedNode, 4, donorProperParts.front());
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

        static_cast<void>(editor.createDetachedNode());
        requireThrows([&]() { editor.commit(); }, "commit must reject a detached node");
        require(tree->isEditing(), "failed commit must leave the edit session open");
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
        const NodeId node = 3;
        const NodeId parent = tree->getNodeParent(node);
        auto editor = tree->edit();
        editor.detach(node);
        editor.attach(parent, node);
        auto proof = editor.proveIncremental();
        require(!proof.usedCompleteValidation(), "supported topology edits must produce an incremental proof");
        editor.commit(std::move(proof));

        const auto& statistics = tree->getEditValidationStatistics();
        requireEqual(statistics.incrementalValidationCommits, std::size_t{1}, "incremental topology proof commit count");
        requireEqual(statistics.completeValidationCommits, std::size_t{0}, "incremental topology edit must not run complete commit validation");
    }

    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), true);
        const NodeId node = 3;
        const NodeId parent = tree->getNodeParent(node);
        auto editor = tree->edit();
        auto staleProof = editor.proveIncremental();
        editor.detach(node);
        requireThrows([&] { editor.commit(std::move(staleProof)); }, "an incremental proof must not commit a later edit revision");
        editor.attach(parent, node);
        editor.commit();
        requireEqual(tree->getEditValidationStatistics().completeValidationCommits, std::size_t{1}, "checked recovery after stale proof");
    }

    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), true);
        auto editor = tree->edit();
        editor.reparent(3, 4);
        requireThrows([&] { static_cast<void>(editor.proveIncremental()); }, "incremental proof must reject a touched parent cycle");
    }

    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), true);
        auto editor = tree->edit();
        requireThrows([&]() { (void)tree->edit(); }, "nested edit sessions must be rejected");
        editor.commit();
        require(!tree->isEditing(), "commit must close a valid edit session");
    }

    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), true);
        const PublicTreeState original = capturePublicTreeState(*tree);
        auto editor = tree->edit();
        editor.detach(3);

        requireThrows(
            [&] {
                MorphologicalTree moved(std::move(*tree));
                static_cast<void>(moved);
            },
            "move construction must reject a tree with an active editor");
        requireThrows([&] { static_cast<void>(tree->clone()); }, "clone must reject a staged forest");
        require(tree->isEditing(), "failed move/clone must leave the editor bound to its original owner");

        editor.rollback();
        requirePublicTreeStateEqual(*tree, original, "failed move/clone followed by rollback");

        MorphologicalTree moved(std::move(*tree));
        moved.validateConnectedRootedTree();
        requireEqual(moved.getNumTotalProperParts(), static_cast<int>(original.properPartOwners.size()), "committed topology must remain movable");
    }

    {
        auto destination = makeComponentTree(makeComponentTreeFixture(), true);
        auto source = makeComponentTree(makeComponentTreeFixture(), false);
        const PublicTreeState destinationBefore = capturePublicTreeState(*destination);
        const PublicTreeState sourceBefore = capturePublicTreeState(*source);

        auto sourceEditor = source->edit();
        sourceEditor.detach(3);
        requireThrows([&] { *destination = std::move(*source); }, "move assignment must reject an actively edited source");
        requirePublicTreeStateEqual(*destination, destinationBefore, "rejected move assignment must preserve its destination");
        require(source->isEditing(), "rejected move assignment must preserve its source editor");
        sourceEditor.rollback();
        requirePublicTreeStateEqual(*source, sourceBefore, "rejected move assignment source rollback");

        auto destinationEditor = destination->edit();

        requireThrows([&] { *destination = std::move(*source); }, "move assignment must reject an actively edited destination");
        requirePublicTreeStateEqual(*source, sourceBefore, "rejected move assignment must preserve its source");
        destinationEditor.commit();

        const NodeId expectedLca = source->getLowestCommonAncestor(3, 4);
        static_cast<void>(destination->getLowestCommonAncestor(3, 4));
        *destination = std::move(*source);
        destination->validateConnectedRootedTree();
        requireEqual(destination->getNumTotalProperParts(), static_cast<int>(sourceBefore.properPartOwners.size()),
                     "committed topology must remain move-assignable");
        requireEqual(destination->getLowestCommonAncestor(3, 4), expectedLca, "move assignment must rebuild an owner-bound LCA cache");
    }

    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), true);
        const NodeId originalParent = tree->getNodeParent(3);
        const int originalSlots = tree->getNumInternalNodeSlots();
        const std::size_t originalVersion = tree->getMutationVersion();
        {
            auto editor = tree->edit();
            require(editor.canRollback(), "ordinary public editor must advertise delta rollback");
            const NodeId inserted = editor.createDetachedNode();
            editor.reparent(3, inserted);
            require(tree->isEditing(), "uncommitted editor must keep the tree in editing mode");
        }
        require(!tree->isEditing(), "destroying an ordinary editor must roll back and close the session");
        requireEqual(tree->getNodeParent(3), originalParent, "ordinary editor rollback must restore parent links");
        requireEqual(tree->getNumInternalNodeSlots(), originalSlots, "ordinary editor rollback must remove appended slots");
        requireEqual(tree->getMutationVersion(), originalVersion, "ordinary editor rollback must restore the committed version");
        requireEqual(computeAreaViaAttributeFacade(*tree, tree->getRoot()), tree->getNumTotalProperParts(),
                     "rolled-back ordinary edit must leave the tree usable");
        auto nextEditor = tree->edit();
        nextEditor.commit();
    }

    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), false);
        const PublicTreeState original = capturePublicTreeState(*tree);
        {
            auto editor = tree->edit();
            editor.moveChildren(4, 3);
            editor.moveProperParts(4, 3);
        }
        requirePublicTreeStateEqual(*tree, original, "abandoned child/proper-part splice");
    }

    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), false);
        const PublicTreeState original = capturePublicTreeState(*tree);
        {
            auto editor = tree->edit();
            const NodeId appended = editor.createDetachedNode();
            const NodeId root = tree->getRoot();
            editor.attach(root, appended);
            editor.removeChild(root, appended, true);
            const NodeId reused = editor.createDetachedNode();
            requireEqual(reused, appended, "edit-local free slot reuse");
            editor.releaseNode(reused);
        }
        requirePublicTreeStateEqual(*tree, original, "abandoned append/release/reuse sequence");
    }

    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), false);
        const PublicTreeState original = capturePublicTreeState(*tree);
        {
            auto editor = tree->edit();
            editor.setRoot(2);
        }
        requirePublicTreeStateEqual(*tree, original, "abandoned root promotion");
    }

    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), false);
        const PublicTreeState original = capturePublicTreeState(*tree);
        try {
            auto editor = tree->edit();
            editor.pruneNode(3);
            throw std::runtime_error("synthetic failure after prune");
        } catch (const std::runtime_error&) {
        }
        requirePublicTreeStateEqual(*tree, original, "exception rollback after subtree prune");
    }

    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), false);
        const PublicTreeState original = capturePublicTreeState(*tree);
        {
            auto editor = tree->edit();
            editor.mergeNodeIntoParent(3);
        }
        requirePublicTreeStateEqual(*tree, original, "abandoned parent merge");
    }

    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), true);
        const NodeId originalParent = tree->getNodeParent(3);
        const int originalSlots = tree->getNumInternalNodeSlots();
        const std::size_t originalVersion = tree->getMutationVersion();

        auto editor = tree->edit();
        require(editor.canRollback(), "ordinary editor must support explicit rollback");
        editor.detach(3);
        static_cast<void>(editor.createDetachedNode());
        editor.rollback();

        require(!tree->isEditing(), "explicit rollback must close the edit session");
        requireEqual(tree->getNodeParent(3), originalParent, "explicit rollback must restore parent links");
        requireEqual(tree->getNumInternalNodeSlots(), originalSlots, "explicit rollback must restore dense node slots");
        requireEqual(tree->getMutationVersion(), originalVersion, "explicit rollback must restore the committed version");
    }

    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), true);
        const NodeId originalParent = tree->getNodeParent(3);
        {
            auto editor = tree->edit();
            editor.reparent(3, 4);
            require(tree->isEditing(), "abandoned transactional edit must be visible while active");
        }
        require(!tree->isEditing(), "transactional editor destructor must close through rollback");
        requireEqual(tree->getNodeParent(3), originalParent, "transactional editor destructor must restore topology");
        requireEqual(computeAreaViaAttributeFacade(*tree, tree->getRoot()), tree->getNumTotalProperParts(), "rolled-back hierarchy must remain usable");
    }

    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), true);
        requireEqual(tree->getProperPartOwner(5), 4, "transaction commit fixture ownership");
        {
            auto editor = tree->edit();
            editor.moveProperPart(5, 4, 5);
            editor.commit();
        }
        requireEqual(tree->getProperPartOwner(5), 5, "successful transactional commit must discard rollback snapshot");
        requireEqual(computeAreaViaAttributeFacade(*tree, tree->getRoot()), tree->getNumTotalProperParts(),
                     "successful transactional edit must publish a valid root support");
    }

    return 0;
}
