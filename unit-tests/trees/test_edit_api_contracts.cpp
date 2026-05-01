#include "support/TestSupport.hpp"
#include "mmcfilters/trees/TreeEditor.hpp"
#include "mmcfilters/trees/WeightedMorphologicalTree.hpp"

#include <stdexcept>
#include <type_traits>
#include <utility>

using namespace mmcfilters;
using namespace mmcfilters::unit_tests;

namespace {

template <class T>
concept HasPublicTree = requires(T& value) {
    value.tree;
};

template <class T>
concept HasPublicTreeUnderscore = requires(T& value) {
    value.tree_;
};

template <class T>
concept HasPublicAltitudeUnderscore = requires(T& value) {
    value.altitude_;
};

template <class T>
concept HasPublicCreateDetachedNode = requires(T& value) {
    value.createDetachedNode();
};

template <class T>
concept HasPublicAttachNode = requires(T& value, NodeId parentId, NodeId nodeId) {
    value.attachNode(parentId, nodeId);
};

template <class T>
concept HasPublicDetachNode = requires(T& value, NodeId nodeId) {
    value.detachNode(nodeId);
};

template <class T>
concept HasPublicMoveNode = requires(T& value, NodeId nodeId, NodeId parentId) {
    value.moveNode(nodeId, parentId);
};

template <class T>
concept HasPublicMoveChildren = requires(T& value, NodeId parentId, NodeId sourceId) {
    value.moveChildren(parentId, sourceId);
};

template <class T>
concept HasPublicMoveProperPart = requires(T& value, NodeId targetId, NodeId sourceId, NodeId properPartId) {
    value.moveProperPart(targetId, sourceId, properPartId);
};

template <class T>
concept HasPublicMoveProperParts = requires(T& value, NodeId targetId, NodeId sourceId) {
    value.moveProperParts(targetId, sourceId);
};

template <class T>
concept WeightedTopologyAllowsEdit = requires(T& weighted) {
    weighted.topology().edit();
};

template <class T>
concept WeightedTopologyAllowsPrune = requires(T& weighted, NodeId nodeId) {
    weighted.topology().pruneNode(nodeId);
};

static_assert(!std::is_constructible_v<TreeEditor, MorphologicalTree&>);
static_assert(!std::is_constructible_v<WeightedTreeEditor, WeightedMorphologicalTree&>);

static_assert(std::is_same_v<decltype(std::declval<MorphologicalTree&>().edit()), TreeEditor>);
static_assert(std::is_same_v<decltype(std::declval<WeightedMorphologicalTree&>().edit()), WeightedTreeEditor>);

static_assert(!std::is_base_of_v<MorphologicalTree, WeightedMorphologicalTree>);
static_assert(!std::is_convertible_v<WeightedMorphologicalTree&, MorphologicalTree&>);
static_assert(std::is_same_v<decltype(std::declval<WeightedMorphologicalTree&>().topology()), const MorphologicalTree&>);
static_assert(std::is_same_v<decltype(std::declval<const WeightedMorphologicalTree&>().topology()), const MorphologicalTree&>);
static_assert(std::is_same_v<decltype(std::declval<const WeightedMorphologicalTree&>().getAltitudeBuffer()), const AltitudeBuffer&>);

static_assert(std::is_invocable_r_v<void, decltype(&MorphologicalTree::pruneNode), MorphologicalTree&, NodeId>);
static_assert(std::is_invocable_r_v<void, decltype(&MorphologicalTree::mergeNodeIntoParent), MorphologicalTree&, NodeId>);
static_assert(std::is_invocable_r_v<void, decltype(&WeightedMorphologicalTree::pruneNode), WeightedMorphologicalTree&, NodeId>);
static_assert(std::is_invocable_r_v<void, decltype(&WeightedMorphologicalTree::mergeNodeIntoParent), WeightedMorphologicalTree&, NodeId>);

static_assert(!HasPublicCreateDetachedNode<MorphologicalTree>);
static_assert(!HasPublicAttachNode<MorphologicalTree>);
static_assert(!HasPublicDetachNode<MorphologicalTree>);
static_assert(!HasPublicMoveNode<MorphologicalTree>);
static_assert(!HasPublicMoveChildren<MorphologicalTree>);
static_assert(!HasPublicMoveProperPart<MorphologicalTree>);
static_assert(!HasPublicMoveProperParts<MorphologicalTree>);

static_assert(!HasPublicTree<WeightedMorphologicalTree>);
static_assert(!HasPublicTreeUnderscore<WeightedMorphologicalTree>);
static_assert(!HasPublicAltitudeUnderscore<WeightedMorphologicalTree>);
static_assert(!HasPublicAttachNode<WeightedMorphologicalTree>);
static_assert(!HasPublicDetachNode<WeightedMorphologicalTree>);
static_assert(!HasPublicMoveNode<WeightedMorphologicalTree>);
static_assert(!HasPublicMoveChildren<WeightedMorphologicalTree>);
static_assert(!HasPublicMoveProperPart<WeightedMorphologicalTree>);
static_assert(!HasPublicMoveProperParts<WeightedMorphologicalTree>);
static_assert(!WeightedTopologyAllowsEdit<WeightedMorphologicalTree>);
static_assert(!WeightedTopologyAllowsPrune<WeightedMorphologicalTree>);

} // namespace

int main() {
    {
        auto tree = makeComponentTree(makeComponentTreeFixture(), true);
        auto editor = tree->edit();

        const NodeId insertedNode = editor.createDetachedNode();
        editor.attach(tree->getRoot(), insertedNode);
        editor.commit();

        requireEqual(tree->getNodeParent(insertedNode), tree->getRoot(), "TreeEditor must be the topology edit entrypoint");
    }

    {
        auto weighted = makeWeightedComponentTree(makeComponentTreeFixture(), true);
        auto editor = weighted->edit();

        const NodeId root = weighted->topology().getRoot();
        const NodeId insertedNode = editor.createDetachedNode(weighted->getAltitude(root) - 1);
        editor.attach(root, insertedNode);

        requireThrows<std::runtime_error>(
            [&]() { editor.commit(); },
            "WeightedTreeEditor commit must reject non-monotone altitude after topology validation");
    }

    return 0;
}
