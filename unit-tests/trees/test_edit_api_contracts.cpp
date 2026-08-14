#include "support/TestSupport.hpp"
#include "mmcfilters/trees/TreeEditor.hpp"
#include "mmcfilters/trees/ValuedMorphologicalTree.hpp"

#include <stdexcept>
#include <type_traits>
#include <utility>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

namespace {

template <class T>
concept HasPublicTree = requires(T& value) { value.tree; };

template <class T>
concept HasPublicTreeUnderscore = requires(T& value) { value.tree_; };

template <class T>
concept HasPublicAltitudeUnderscore = requires(T& value) { value.altitude_; };

template <class T>
concept HasPublicCreateDetachedNode = requires(T& value) { value.createDetachedNode(); };

template <class T>
concept HasPublicAttachNode = requires(T& value, NodeId parentId, NodeId nodeId) { value.attachNode(parentId, nodeId); };

template <class T>
concept HasPublicDetachNode = requires(T& value, NodeId nodeId) { value.detachNode(nodeId); };

template <class T>
concept HasPublicMoveNode = requires(T& value, NodeId nodeId, NodeId parentId) { value.moveNode(nodeId, parentId); };

template <class T>
concept HasPublicMoveChildren = requires(T& value, NodeId parentId, NodeId sourceId) { value.moveChildren(parentId, sourceId); };

template <class T>
concept HasPublicMoveProperPart =
    requires(T& value, NodeId targetId, NodeId sourceId, PixelId pixel) { value.movePixelToProperPart(targetId, sourceId, pixel); };

template <class T>
concept HasPublicMoveProperParts = requires(T& value, NodeId targetId, NodeId sourceId) { value.mergeProperParts(targetId, sourceId); };

template <class T>
concept ValuedTopologyAllowsEdit = requires(T& valuedTree) { valuedTree.topology().edit(); };

template <class T>
concept ValuedTopologyAllowsPrune = requires(T& valuedTree, NodeId nodeId) { valuedTree.topology().pruneNode(nodeId); };

template <class T>
concept HasPublicUncheckedAltitude = requires(T& valuedTree, NodeId nodeId, std::uint8_t altitude) { valuedTree.setAltitudeUnchecked(nodeId, altitude); };

template <class T>
concept HasPublicUncheckedAltitudeBuffer =
    requires(T& valuedTree, NodeAltitudeBuffer<std::uint8_t> altitude) { valuedTree.setAltitudeBufferUnchecked(std::move(altitude)); };

static_assert(!std::is_constructible_v<TreeEditor, MorphologicalTree&>);
static_assert(!std::is_constructible_v<ValuedMorphologicalTreeEditor<std::uint8_t>, ValuedMorphologicalTree<std::uint8_t>&>);

static_assert(std::is_same_v<decltype(std::declval<MorphologicalTree&>().edit()), TreeEditor>);
static_assert(std::is_same_v<decltype(std::declval<ValuedMorphologicalTree<std::uint8_t>&>().edit()), ValuedMorphologicalTreeEditor<std::uint8_t>>);

static_assert(!std::is_base_of_v<MorphologicalTree, ValuedMorphologicalTree<std::uint8_t>>);
static_assert(!std::is_convertible_v<ValuedMorphologicalTree<std::uint8_t>&, MorphologicalTree&>);
static_assert(std::is_same_v<decltype(std::declval<ValuedMorphologicalTree<std::uint8_t>&>().topology()), const MorphologicalTree&>);
static_assert(std::is_same_v<decltype(std::declval<const ValuedMorphologicalTree<std::uint8_t>&>().topology()), const MorphologicalTree&>);
static_assert(
    std::is_same_v<decltype(std::declval<const ValuedMorphologicalTree<std::uint8_t>&>().nodeAltitudes()), const NodeAltitudeBuffer<std::uint8_t>&>);

static_assert(std::is_invocable_r_v<void, decltype(&MorphologicalTree::pruneNode), MorphologicalTree&, NodeId>);
static_assert(std::is_invocable_r_v<void, decltype(&MorphologicalTree::mergeNodeIntoParent), MorphologicalTree&, NodeId>);
static_assert(std::is_invocable_r_v<void, decltype(&ValuedMorphologicalTree<std::uint8_t>::pruneNode), ValuedMorphologicalTree<std::uint8_t>&, NodeId>);
static_assert(
    std::is_invocable_r_v<void, decltype(&ValuedMorphologicalTree<std::uint8_t>::mergeNodeIntoParent), ValuedMorphologicalTree<std::uint8_t>&, NodeId>);
static_assert(std::is_invocable_r_v<void, decltype(&ValuedMorphologicalTreeEditor<std::uint8_t>::pruneNode), ValuedMorphologicalTreeEditor<std::uint8_t>&, NodeId>);
static_assert(std::is_invocable_r_v<void, decltype(&ValuedMorphologicalTreeEditor<std::uint8_t>::mergeNodeIntoParent), ValuedMorphologicalTreeEditor<std::uint8_t>&, NodeId>);

static_assert(!HasPublicCreateDetachedNode<MorphologicalTree>);
static_assert(!HasPublicAttachNode<MorphologicalTree>);
static_assert(!HasPublicDetachNode<MorphologicalTree>);
static_assert(!HasPublicMoveNode<MorphologicalTree>);
static_assert(!HasPublicMoveChildren<MorphologicalTree>);
static_assert(!HasPublicMoveProperPart<MorphologicalTree>);
static_assert(!HasPublicMoveProperParts<MorphologicalTree>);

static_assert(!HasPublicTree<ValuedMorphologicalTree<std::uint8_t>>);
static_assert(!HasPublicTreeUnderscore<ValuedMorphologicalTree<std::uint8_t>>);
static_assert(!HasPublicAltitudeUnderscore<ValuedMorphologicalTree<std::uint8_t>>);
static_assert(!HasPublicAttachNode<ValuedMorphologicalTree<std::uint8_t>>);
static_assert(!HasPublicDetachNode<ValuedMorphologicalTree<std::uint8_t>>);
static_assert(!HasPublicMoveNode<ValuedMorphologicalTree<std::uint8_t>>);
static_assert(!HasPublicMoveChildren<ValuedMorphologicalTree<std::uint8_t>>);
static_assert(!HasPublicMoveProperPart<ValuedMorphologicalTree<std::uint8_t>>);
static_assert(!HasPublicMoveProperParts<ValuedMorphologicalTree<std::uint8_t>>);
static_assert(!ValuedTopologyAllowsEdit<ValuedMorphologicalTree<std::uint8_t>>);
static_assert(!ValuedTopologyAllowsPrune<ValuedMorphologicalTree<std::uint8_t>>);
static_assert(!HasPublicUncheckedAltitude<ValuedMorphologicalTree<std::uint8_t>>);
static_assert(!HasPublicUncheckedAltitudeBuffer<ValuedMorphologicalTree<std::uint8_t>>);

} // namespace

int main() {
    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), true);
        auto editor = tree->edit();

        const NodeId root = tree->root();
        const NodeId formerChild = tree->getFirstChild(root);
        const NodeId insertedNode = editor.createDetachedNode();
        editor.attach(root, insertedNode);
        requireThrows<std::runtime_error>([&]() { editor.commit(); }, "TreeEditor commit must reject an attached node with empty subtree support");
        editor.reparent(formerChild, insertedNode);
        editor.commit();

        requireEqual(tree->parent(insertedNode), root, "TreeEditor must be the topology edit entrypoint");
        require(tree->hasEmptyProperPart(insertedNode), "inserted hierarchy-only node must be derived as structural");
    }

    {
        auto valuedTree = makeValuedComponentTree(makeComponentTreeFixture(), false);
        auto spacedAltitude = valuedTree->nodeAltitudes();
        for (auto& value : spacedAltitude) {
            value = static_cast<std::uint8_t>(value * 2);
        }
        valuedTree->setNodeAltitudes(std::move(spacedAltitude));
        auto editor = valuedTree->edit();

        const NodeId root = valuedTree->topology().root();
        const NodeId formerChild = valuedTree->topology().getFirstChild(root);
        const NodeId insertedNode = editor.createDetachedNode(static_cast<std::uint8_t>(valuedTree->nodeAltitude(root) + 1));
        editor.reparent(formerChild, insertedNode);
        editor.attach(root, insertedNode);

        requireThrows<std::runtime_error>([&]() { editor.commit(); },
                                          "ValuedMorphologicalTreeEditor<std::uint8_t> commit must reject non-monotone altitude after topology validation");
        require(valuedTree->topology().isEditing(), "failed valuedTree commit must leave the edit session open");

        editor.setNodeAltitude(insertedNode, static_cast<std::uint8_t>(valuedTree->nodeAltitude(root) - 1));
        const TreeValidationResult repaired = editor.validateAndCommit();
        require(repaired.ok, "ValuedMorphologicalTreeEditor<std::uint8_t> validateAndCommit must accept repaired altitude monotonicity");
        require(!valuedTree->topology().isEditing(), "successful valuedTree validateAndCommit must close the edit session");
    }

    {
        auto valuedTree = makeValuedComponentTree(makeComponentTreeFixture(), false);
        auto spacedAltitude = valuedTree->nodeAltitudes();
        for (auto& value : spacedAltitude) {
            value = static_cast<std::uint8_t>(value * 2);
        }
        valuedTree->setNodeAltitudes(std::move(spacedAltitude));
        auto editor = valuedTree->edit();

        const NodeId root = valuedTree->topology().root();
        const NodeId formerChild = valuedTree->topology().getFirstChild(root);
        const NodeId insertedNode = editor.createDetachedNode(static_cast<std::uint8_t>(valuedTree->nodeAltitude(root) + 1));
        editor.reparent(formerChild, insertedNode);
        editor.attach(root, insertedNode);

        requireThrows<std::runtime_error>([&] { static_cast<void>(editor.proveIncremental()); },
                                          "incremental valuedTree proof must reject non-monotone touched arcs");

        editor.setNodeAltitude(insertedNode, static_cast<std::uint8_t>(valuedTree->nodeAltitude(root) - 1));
        auto proof = editor.proveIncremental();
        require(!proof.usedCompleteValidation(), "supported valuedTree edit must use delta validation");
        editor.commit(std::move(proof));

        const auto& statistics = valuedTree->topology().getEditValidationStatistics();
        requireEqual(statistics.incrementalValidationCommits, std::size_t{1}, "valuedTree incremental proof commit count");
        requireEqual(statistics.completeValidationCommits, std::size_t{0}, "valuedTree incremental proof must avoid complete commit validation");
    }

    return 0;
}
