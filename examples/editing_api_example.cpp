/**
 * Demonstrate the public editing boundary for topology-only and valuedTree
 * morphological trees.
 *
 * Build with `-DMMCFILTERS_BUILD_EXAMPLES=ON` and run:
 * `./build/examples/mmcfilters_example_editing_api`.
 *
 * The example covers safe public mutators, staged `TreeEditor` commits,
 * valuedTree altitude insertion, and rejection of a non-monotone valuedTree edit.
 */
#include "../mmcfilters/trees/TreeEditor.hpp"
#include "../mmcfilters/trees/ValuedMorphologicalTree.hpp"
#include "../mmcfilters/trees/MorphologicalTreeFactory.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {

mmcfilters::ImageUInt8Ptr makeFixtureImage() {
    const std::vector<uint8_t> values = {
        3, 3, 2, 2, 3, 4, 4, 2, 1, 4, 5, 2, 1, 1, 5, 0,
    };

    auto image = mmcfilters::ImageUInt8::create(4, 4);
    for (std::size_t i = 0; i < values.size(); ++i) {
        (*image)[static_cast<int>(i)] = values[i];
    }
    return image;
}

void require(bool condition, const char* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

mmcfilters::NodeId firstNonRootLeaf(const mmcfilters::MorphologicalTree& tree) {
    for (mmcfilters::NodeId nodeId : tree.leaves()) {
        if (!tree.isRoot(nodeId)) {
            return nodeId;
        }
    }
    throw std::runtime_error("fixture must expose at least one non-root leaf");
}

void safeMutatorExample() {
    auto topologySource = mmcfilters::MorphologicalTreeFactory::createMaxTree(makeFixtureImage());
    auto tree = topologySource.topology().clone();
    const int nodesBeforePrune = tree.numNodes();

    tree.pruneNode(firstNonRootLeaf(tree));
    require(tree.numNodes() < nodesBeforePrune, "pruneNode should remove a live subtree");

    auto valuedTree = mmcfilters::MorphologicalTreeFactory::createMaxTree(makeFixtureImage());
    const int nodesBeforeMerge = valuedTree.topology().numNodes();

    valuedTree.mergeNodeIntoParent(firstNonRootLeaf(valuedTree.topology()));
    require(valuedTree.topology().numNodes() < nodesBeforeMerge, "valuedTree merge should update the owned topology");
}

void treeEditorExample() {
    auto topologySource = mmcfilters::MorphologicalTreeFactory::createMinTree(makeFixtureImage());
    auto tree = topologySource.topology().clone();
    auto editor = tree.edit();

    const mmcfilters::NodeId insertedNode = editor.createDetachedNode();
    editor.reparent(3, insertedNode);
    editor.reparent(4, insertedNode);
    editor.attach(2, insertedNode);
    editor.commit();

    require(tree.parent(insertedNode) == 2, "TreeEditor should validate and commit the staged topology");
}

void valuedTreeEditorExample() {
    auto valuedTree = mmcfilters::MorphologicalTreeFactory::createMinTree(makeFixtureImage());
    auto editor = valuedTree.edit();

    const std::uint8_t insertedAltitude = std::min(valuedTree.nodeAltitude(2), std::max(valuedTree.nodeAltitude(3), valuedTree.nodeAltitude(4)));
    const mmcfilters::NodeId insertedNode = editor.createDetachedNode(insertedAltitude);

    editor.reparent(3, insertedNode);
    editor.reparent(4, insertedNode);
    editor.attach(2, insertedNode);
    editor.commit();

    require(valuedTree.nodeAltitude(insertedNode) == insertedAltitude, "ValuedMorphologicalTreeEditor<std::uint8_t> should preserve inserted altitude");
}

void rejectedValuedTreeCommitExample() {
    auto valuedTree = mmcfilters::MorphologicalTreeFactory::createMinTree(makeFixtureImage());
    auto editor = valuedTree.edit();

    const mmcfilters::NodeId root = valuedTree.topology().root();
    const auto invalidAltitude = static_cast<std::uint8_t>(valuedTree.nodeAltitude(root) + 1);
    const mmcfilters::NodeId insertedNode = editor.createDetachedNode(invalidAltitude);
    editor.attach(root, insertedNode);

    const auto result = editor.validateAndCommit();
    require(!result.ok, "ValuedMorphologicalTreeEditor<std::uint8_t>::validateAndCommit should reject non-monotone altitude");
    require(valuedTree.topology().isEditing(), "failed validateAndCommit should keep the edit session open");
}

} // namespace

int main() {
    safeMutatorExample();
    treeEditorExample();
    valuedTreeEditorExample();
    rejectedValuedTreeCommitExample();

    std::cout << "Editing API example completed.\n";
    return 0;
}
